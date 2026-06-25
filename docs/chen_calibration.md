# Chen calibration fields — what they mean and how they drive the carrier call

The columns `chen_evidence`, `chen_calibration_approach`, `chen_domain`,
`chen_vep_score`, and `in_chen_table` in `A_variant_level_per_person.tsv`
(and downstream tables) all come from one external table:

> **Chen et al., 2026** (the "Chen/Pejaver paper"). They re-calibrated
> AlphaMissense and other in-silico predictors against the ACMG/AMP 2015
> evidence framework, on a per-gene basis, using ClinVar-validated
> pathogenic/benign labels. The output is a per-variant table that says:
> "given this variant's AM score, here is the evidence level you may
> claim under PP3/BP4, for THIS gene specifically."

That per-variant lookup is the input to our carrier call. The fields
below explain what each Chen column tells you, and how the pipeline
converts them into the boolean `is_AM_carrier_primary`.

---

## `in_chen_table` — was this variant in Chen's calibration set?

Two values:

- **TRUE** — `(CHROM, POS, REF, ALT)` matched a row in
  `refs/chen_calibration.target_genes.tsv.gz`. Chen has a specific
  evidence label and points score for this exact variant.
- **FALSE** — no match. The pipeline still records the AlphaMissense
  score (`am_pathogenicity`, `am_class`) but cannot say what evidence
  level Chen's calibration assigns. Such variants are **not** counted
  as AM carriers under the primary policy — they will have
  `is_AM_carrier_primary = FALSE` regardless of how high the AM score is.

Why "not counted": the entire point of Chen's calibration is to
replace the original single-cutoff (0.864) treatment of AM with a
*per-variant, ACMG/AMP-anchored* evidence label. A variant Chen
didn't score has no calibrated label, so its AM score has no
calibration anchor to stand on.

---

## `chen_evidence` — the ACMG/AMP evidence label Chen assigned

The full enum (case-preserving):

| value | meaning (ACMG/AMP semantics) |
|---|---|
| `PP3_VeryStrong` | strongest "computational evidence supports pathogenic" — large effect |
| `PP3_Strong` | strong PP3 evidence |
| `PP3_Moderate` | moderate PP3 evidence — the **threshold our primary policy uses** |
| `PP3_Supporting` | weak PP3 evidence — does **not** qualify for the primary policy |
| `PP3_<level>+` | same as the `<level>` above, but Chen's "+" marker means the variant is *at or above* the cutoff with margin to spare. Still counted under the corresponding non-+ tier. |
| `BP4_Strong` / `BP4_Moderate` / `BP4_Supporting` | "computational evidence supports benign" — opposite sign. Never counted as a carrier. |
| `` (blank / empty) | Chen has the variant but it falls in the neutral middle band (neither PP3 nor BP4 was claimed). Not counted as a carrier. |

The mapping from `chen_evidence` to the primary carrier flag is
**explicit and conservative**, defined in `config/config.yaml` as
`calibration.min_evidence_strength: Moderate`:

```
is_AM_carrier_primary = TRUE
  iff  chen_evidence ∈ {PP3_Moderate, PP3_Moderate+,
                        PP3_Strong,   PP3_Strong+,
                        PP3_VeryStrong, PP3_VeryStrong+}
```

Anything weaker (`PP3_Supporting`), anything benign-direction
(`BP4_*`), and anything blank ⇒ `FALSE`.

### "Does `single_gene PP3≥Moderate` mean AM-deleterious after Chen calibration?"

**Yes — that's exactly the operational definition.** Plain-language
statement of the primary carrier flag:

> A variant is called an **AM carrier (primary policy)** if Chen's
> per-gene, ACMG/AMP-anchored calibration assigns it `PP3_Moderate`
> or stronger evidence, using the *single_gene* calibration approach
> (see next section). That is the level at which AlphaMissense's
> score reaches the ACMG/AMP "Moderate" pathogenic threshold for that
> specific gene, given Chen's ClinVar-anchored fit.

A few corollaries worth being explicit about:

- The threshold is **gene-specific** in two ways: (i) Chen fit a
  separate logistic regression per gene where possible, so the AM
  score required to reach PP3_Moderate differs gene by gene;
  (ii) the threshold is not the global 0.864 default — that global
  number was Chen's pooled estimate and over-claimed in genes with
  noisier ClinVar truth sets.
- `is_AM_carrier_primary` is **a carrier definition, not a clinical
  classification**. ACMG/AMP 2015 does not permit a single in-silico
  tool to drive variant classification on its own; the column is
  meant to support *cohort-level association testing* and *candidate
  prioritization*, not lab reporting.

---

## `chen_calibration_approach` — how Chen got to a threshold

Two values:

- **`single_gene`** — Chen had enough labeled variants in this gene
  alone to fit a per-gene calibration. The PP3 thresholds for this
  gene rest entirely on that gene's own labeled data. This is the
  *strongest* calibration tier and is what the primary policy uses.
- **`domain_aggregate`** — Chen did **not** have enough labeled
  variants in this gene alone. Instead the calibration is borrowed
  from a *protein-domain class* (e.g. all kinase domains, all WD40
  repeats) pooled across genes. The threshold still applies, but the
  evidence is weaker — it's "variants in this domain class behave
  this way" rather than "variants in this gene behave this way."

How the pipeline uses this:

- `is_AM_carrier_primary` only considers `single_gene` rows. A
  `PP3_Moderate` call from a `domain_aggregate` row is **not**
  counted under the primary policy.
- `would_be_carrier_domain_aggregate` is the informational column
  that *would* flip TRUE under a looser policy that accepts
  `domain_aggregate` rows. Used only for sensitivity analyses.
- `threshold_source` records which path was taken, with three values:

| `threshold_source` | meaning |
|---|---|
| `single_gene` | variant in Chen, calibration_approach=single_gene; primary policy can act on it |
| `domain_aggregate` | variant in Chen, calibration_approach=domain_aggregate; primary policy skips it |
| `not_in_chen_table` | variant absent from Chen entirely; AM score retained for reference only |

---

## `chen_domain` — which protein domain class

A short tag like `WD40`, `kinase_domain`, `BRCT`, or `no_pfam` (if
the variant doesn't fall in a Pfam-annotated domain). Two uses:

- **For `domain_aggregate` rows**: this is the domain class whose
  pooled calibration drove the threshold. Knowing the class lets you
  audit whether the borrowing is biologically reasonable (e.g.
  borrowing across kinase domains is more defensible than borrowing
  across "no_pfam").
- **For `single_gene` rows**: still recorded for transparency, but
  the calibration didn't depend on it.

---

## `chen_vep_score` — the AM score Chen scored the variant with

The AlphaMissense pathogenicity score at the time Chen built the
calibration table. We carry it through for **auditability**:

- It should match `am_pathogenicity` for the same variant up to
  release drift. If they don't match, AlphaMissense has been
  re-released since Chen's calibration was fit — flag it.
- Lets a reviewer check: "for this PP3_Moderate call, what AM score
  is Chen's calibration claiming reaches the moderate threshold for
  this gene?"

---

## How these fields are produced

`scripts/00_prepare_refs.sh` downloads Chen's calibration TSV (from
Zenodo by default, or a local file if `calibration.local_file` is
set), filters to the 28 target genes, and writes
`refs/chen_calibration.target_genes.tsv.gz` (+ tabix index).

`scripts/02_annotate_am.lsf` calls `python/annotate_am.py`, which
joins every observed missense variant against both the AlphaMissense
subset (`refs/AlphaMissense_hg38.subset_targets.tsv.gz`) and the
Chen subset, on `(CHROM, POS, REF, ALT)`. The Chen fields above are
filled from the matching Chen row; if no Chen match exists,
`in_chen_table = FALSE` and the Chen fields are blank.

`scripts/03_call_carriers.lsf` applies the threshold policy
described above to derive `is_AM_carrier_primary` and
`would_be_carrier_domain_aggregate`.

The exact threshold-policy source is in `python/call_carriers.py`,
controlled by these `config/config.yaml` keys:

```yaml
calibration:
  min_evidence_strength: "Moderate"     # PP3 floor for primary policy
  primary_threshold: "gene_specific"    # not the legacy single 0.864
```

If a reviewer asks "what's the threshold for gene X?", the answer
is: it's the gene-specific Chen calibration encoded per-variant in
`chen_evidence`. There is no single number — the threshold is a
*decision rule* on Chen's labeled output, not an AM score cutoff.

---

## TL;DR for the rebuttal letter

- AlphaMissense gives a continuous pathogenicity score (0–1).
- A single global threshold (the 0.864 in our original analysis)
  treats every gene the same, which over-claims in noisy genes.
- Chen 2026 re-calibrated AM against ACMG/AMP evidence levels on a
  per-gene basis (or per-domain class where per-gene data was
  insufficient).
- For every observed variant in our 28-gene panel, the pipeline
  records Chen's calibrated evidence label. The **primary AM-carrier
  definition** is "Chen labeled this variant `PP3_Moderate` or
  stronger via the single_gene calibration approach."
- Variants outside Chen's table, or only callable via
  `domain_aggregate`, are tracked separately for sensitivity
  analyses — they do **not** contribute to the primary carrier set.
