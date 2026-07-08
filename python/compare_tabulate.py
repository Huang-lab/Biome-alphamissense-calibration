"""Emit the three deliverable tables (A, B, C) and the SUMMARY.

Inputs:
- Table A from step 03 (per-(sample,variant), with carrier columns)
- phenotype TSV (per-cohort id_column)
- PCs TSV (id column auto-discovered)
- ACMG P/LP TSV (filter to rows whose annotation contains "P/LP")

Outputs (under results/<cohort>/):
- A_variant_level_per_person.tsv      <- copied/forwarded from step 03 input
- B_ACMG_vs_AM_comparison.tsv         <- 3-way carrier list: union of ACMG-PLP,
                                         AM-primary, and ClinVar-P/LP carriers
                                         (all variant types incl. ClinVar indels)
- B_summary_by_gene.tsv               <- per-gene counts AM_only / ACMG_only / both
- C_regression_matrix.tsv             <- one row per sample, regression-ready

Plus results/SUMMARY.md (aggregated across cohorts) is appended to by the LSF
04 driver (this script writes a per-cohort summary; the driver concatenates).
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Set, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from util import LOG, die, load_config, read_tsv_dicts, report_join, resolve, write_tsv  # noqa: E402


def _truthy(v: str) -> bool:
    return (v or "").strip().upper() in ("TRUE", "T", "YES", "Y", "1")


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--a-table", required=True, help="A_variant_level_per_person.tsv from step 03")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    results_dir = resolve(cfg, cfg["paths"]["results_dir"])
    intermediate_dir = resolve(cfg, cfg["paths"]["intermediate_dir"])
    out_dir = os.path.join(results_dir, args.cohort)
    os.makedirs(out_dir, exist_ok=True)

    cohort_cfg = cfg["cohorts"][args.cohort]
    id_col = cohort_cfg["id_column"]
    pheno_path = cohort_cfg["phenotype_tsv"]
    acmg_path  = cohort_cfg["acmg_tsv"]
    pcs_path   = cfg["pcs"]["pcs_tsv"]
    pc_id_override = cfg["pcs"].get("id_column") or ""
    pc_cols = list(cfg["pcs"]["pc_columns"])
    acmg_substr = cfg["acmg"]["annotation_filter_substr"]
    acmg_anno_col = cfg["acmg"]["annotation_column"]
    groups: Dict[str, List[str]] = cfg["gene_groups"]
    target_genes = list(cfg["target_genes"])

    # ---- copy/forward Table A ----------------------------------------------
    a_dst = os.path.join(out_dir, "A_variant_level_per_person.tsv")
    if os.path.abspath(args.a_table) != os.path.abspath(a_dst):
        shutil.copyfile(args.a_table, a_dst)
    a_header, a_rows = read_tsv_dicts(a_dst)
    LOG.info("Table A: %d rows", len(a_rows))

    # ---- phenotype ---------------------------------------------------------
    if not os.path.isfile(pheno_path):
        die(f"phenotype TSV not found: {pheno_path}")
    pheno_header, pheno_rows = read_tsv_dicts(pheno_path)
    if id_col not in pheno_header:
        die(f"id_column {id_col!r} not in phenotype header: {pheno_header}")
    pheno_by_id: Dict[str, dict] = {(r.get(id_col) or "").strip(): r for r in pheno_rows if r.get(id_col)}
    LOG.info("phenotype: %d unique %s", len(pheno_by_id), id_col)

    # ---- join check: VCF samples (from A) vs phenotype ---------------------
    vcf_sample_ids = sorted({r["sample_id"] for r in a_rows})
    pheno_ids = list(pheno_by_id.keys())
    rep = report_join(vcf_sample_ids, pheno_ids, label=f"{args.cohort}:VCF↔phenotype on {id_col}")
    rep.log(f"{args.cohort}:VCF↔phenotype")

    # ---- ID bridge (cohort I's PCs key on MASKED_MRN; phenotype keys on
    # SINAI_ID; bridge file maps one to the other). Cohort II opts out. ------
    bridge_cfg = cfg.get("id_bridge") or {}
    bridge_used = False
    sinai_to_masked: Dict[str, str] = {}
    if bridge_cfg.get("use_for", {}).get(args.cohort, False):
        bridge_path = bridge_cfg.get("masked_mrn_map", "") or ""
        if not bridge_path or not os.path.isfile(bridge_path):
            LOG.warning("id_bridge.use_for.%s=true but bridge file missing at %r; "
                        "PC join will fall back to direct lookup", args.cohort, bridge_path)
        else:
            sinai_col = bridge_cfg.get("map_sinai_id_col", "RGNID")
            mrn_col = bridge_cfg.get("map_masked_mrn_col", "MASKED_MRN")
            _, br_rows = read_tsv_dicts(bridge_path)
            for br in br_rows:
                sid = (br.get(sinai_col) or "").strip()
                mrn = (br.get(mrn_col) or "").strip()
                if sid and mrn:
                    sinai_to_masked[sid] = mrn
            bridge_used = True
            LOG.info("id_bridge %s: loaded %d %s->%s pairs from %s",
                     args.cohort, len(sinai_to_masked), sinai_col, mrn_col, bridge_path)

    # ---- PCs ---------------------------------------------------------------
    pc_header: List[str] = []
    pc_rows: List[dict] = []
    pc_by_id: Dict[str, dict] = {}
    pc_id_col: Optional[str] = None
    rep_pc = None
    if not os.path.isfile(pcs_path):
        LOG.warning("PCs file not found at %s; C matrix will leave PC1..PC10 empty", pcs_path)
    else:
        pc_header, pc_rows = read_tsv_dicts(pcs_path)
        if pc_id_override and pc_id_override in pc_header:
            pc_id_col = pc_id_override
        else:
            # try to discover by maximum overlap with phenotype IDs
            best_col, best_overlap = None, 0
            phen_set = set(pheno_ids)
            for col in pc_header:
                vals = {(r.get(col) or "").strip() for r in pc_rows[:200]}
                overlap = len(vals & phen_set)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_col = col
            pc_id_col = best_col
            LOG.info("PC id column auto-discovered: %r (overlap=%d in first 200 rows)", pc_id_col, best_overlap)
        if pc_id_col:
            pc_by_id = {(r.get(pc_id_col) or "").strip(): r for r in pc_rows if r.get(pc_id_col)}
            # If we're routing through the bridge, the PC join keys are MASKED_MRN
            # values, so report match-rate against the bridge image of phenotype
            # IDs (not the SINAI_IDs themselves).
            if bridge_used:
                phen_in_pc_space = [sinai_to_masked.get(s, "") for s in pheno_ids]
                phen_in_pc_space = [s for s in phen_in_pc_space if s]
                rep_pc = report_join(
                    phen_in_pc_space, list(pc_by_id.keys()),
                    label=f"{args.cohort}:phenotype(via_bridge)↔PCs on {pc_id_col}",
                )
                rep_pc.log(f"{args.cohort}:phenotype(via_bridge)↔PCs")
            else:
                rep_pc = report_join(pheno_ids, list(pc_by_id.keys()),
                                     label=f"{args.cohort}:phenotype↔PCs on {pc_id_col}")
                rep_pc.log(f"{args.cohort}:phenotype↔PCs")

    # ---- ACMG --------------------------------------------------------------
    acmg_carrier_by_gene: Dict[str, Set[str]] = defaultdict(set)
    acmg_rows_kept = 0
    acmg_rows_excluded_ptv = 0
    if os.path.isfile(acmg_path):
        acmg_header, acmg_rows = read_tsv_dicts(acmg_path)
        if acmg_anno_col not in acmg_header:
            die(f"ACMG annotation column {acmg_anno_col!r} not in {acmg_path} header: {acmg_header}")
        # figure out gene + sample id columns by heuristics
        sample_col = next((c for c in acmg_header if c.upper() in (id_col.upper(), "SAMPLE", "SAMPLE_ID")), None)
        if sample_col is None:
            # fallback to id_col exactly
            sample_col = id_col if id_col in acmg_header else None
        gene_col = next((c for c in acmg_header if c.lower() in ("gene", "gene_name", "symbol")), None)
        if not sample_col or not gene_col:
            die(f"ACMG file {acmg_path} must contain sample-id column ({id_col!r} or SAMPLE) and gene column (got header {acmg_header})")
        for r in acmg_rows:
            anno = (r.get(acmg_anno_col) or "")
            if acmg_substr not in anno:
                if "PTV" in anno:
                    acmg_rows_excluded_ptv += 1
                continue
            acmg_carrier_by_gene[(r.get(gene_col) or "").upper()].add((r.get(sample_col) or "").strip())
            acmg_rows_kept += 1
        LOG.info("ACMG: kept %d P/LP rows; excluded %d PTV-only rows (PTV is not comparable to AM-missense)",
                 acmg_rows_kept, acmg_rows_excluded_ptv)
    else:
        LOG.warning("ACMG file not found at %s; B comparison will be empty", acmg_path)

    # ---- ClinVar P/LP >=2* carriers (standalone comparator; ALL variant types) -
    # Built by python/annotate_clinvar.py (LSF step 02c) from the all-variant QC
    # VCF, so it includes non-SNV (indel/MNV) P/LP carriers. Forwarded to results/
    # so the figure layer can read it. Two views:
    #   clinvar_carrier_by_gene : gene_u -> {sample}                (summary)
    #   clinvar_by_variant      : (sample,gene_u,chr,pos,ref,alt) -> {clnsig,stars}
    #                             (variant-level, for the 3-way Table B carrier list)
    clinvar_carrier_by_gene: Dict[str, Set[str]] = defaultdict(set)
    clinvar_by_variant: Dict[Tuple[str, str, str, str, str, str], dict] = {}
    clinvar_rows_kept = 0
    clinvar_src = os.path.join(intermediate_dir, args.cohort, "clinvar_carriers.tsv")
    clinvar_dst = os.path.join(out_dir, "clinvar_carriers.tsv")
    if os.path.isfile(clinvar_src):
        if os.path.abspath(clinvar_src) != os.path.abspath(clinvar_dst):
            shutil.copyfile(clinvar_src, clinvar_dst)
        _cv_header, cv_rows = read_tsv_dicts(clinvar_dst)
        for r in cv_rows:
            gene_u = (r.get("gene") or "").upper()
            sample = (r.get("sample_id") or "").strip()
            if not (gene_u and sample):
                continue
            clinvar_carrier_by_gene[gene_u].add(sample)
            vkey = (sample, gene_u, r.get("chr", ""), r.get("pos", ""),
                    r.get("ref", ""), r.get("alt", ""))
            clinvar_by_variant[vkey] = {
                "clnsig":       r.get("clinvar_clnsig", ""),
                "review_stars": r.get("clinvar_review_stars", ""),
            }
            clinvar_rows_kept += 1
        LOG.info("ClinVar: loaded %d P/LP≥2★ carrier rows from %s", clinvar_rows_kept, clinvar_src)
    else:
        LOG.warning("ClinVar carrier file not found at %s; ClinVar carriers will be empty "
                    "(run scripts/02c_annotate_clinvar.lsf + gather first)", clinvar_src)

    # ---- AM-primary + AM-global-0864 carriers by gene ----------------------
    am_primary_by_gene: Dict[str, Set[str]] = defaultdict(set)
    am_global_0864_by_gene: Dict[str, Set[str]] = defaultdict(set)
    for r in a_rows:
        gene_u = (r.get("gene") or "").upper()
        sample = r.get("sample_id", "")
        if _truthy(r.get("is_AM_carrier_primary", "")):
            am_primary_by_gene[gene_u].add(sample)
        if _truthy(r.get("is_AM_carrier_global_0864", "")):
            am_global_0864_by_gene[gene_u].add(sample)

    # ---- B_ACMG_vs_AM_comparison — 3-way carrier list ----------------------
    # One row per (sample, gene, variant) that is carried by ANY of the three
    # tracks: ACMG P/LP, AM_calibrated (primary), or ClinVar P/LP ≥2★. Includes
    # ALL variant types — ClinVar indels/MNVs come in via clinvar_by_variant even
    # though the AM Table A is SNV-only. Membership flags + ClinVar clnsig/stars
    # + a combined `category` label ("+"-joined of the tracks present).
    def _category(in_acmg: bool, in_am: bool, in_clinvar: bool) -> str:
        parts = []
        if in_acmg:    parts.append("ACMG_PLP")
        if in_am:      parts.append("AM")
        if in_clinvar: parts.append("ClinVar")
        return "+".join(parts) if parts else "none"

    b_header = list(a_header) + ["in_ACMG_PLP", "in_AM_primary", "in_ClinVar_PLP",
                                 "clinvar_clnsig", "clinvar_review_stars", "category"]
    b_rows: List[List[str]] = []
    seen_keys: Set[Tuple[str, str, str, str, str, str]] = set()  # (sample,gene,chr,pos,ref,alt)

    # 1) AM Table A rows (SNV): flag ACMG (gene-level) + ClinVar (variant-level).
    for r in a_rows:
        gene_u = (r.get("gene") or "").upper()
        sample = r.get("sample_id", "")
        vkey = (sample, gene_u, r.get("chr", ""), r.get("pos", ""),
                r.get("ref", ""), r.get("alt", ""))
        in_am = _truthy(r.get("is_AM_carrier_primary", ""))
        in_acmg = sample in acmg_carrier_by_gene.get(gene_u, set())
        cv = clinvar_by_variant.get(vkey)
        in_clinvar = cv is not None
        if not (in_am or in_acmg or in_clinvar):
            continue
        b_rows.append([r.get(h, "") for h in a_header] + [
            "yes" if in_acmg else "no", "yes" if in_am else "no",
            "yes" if in_clinvar else "no",
            cv["clnsig"] if cv else "", cv["review_stars"] if cv else "",
            _category(in_acmg, in_am, in_clinvar),
        ])
        seen_keys.add(vkey)

    # 2) ClinVar carriers not represented in Table A (e.g. non-SNV variants the
    #    SNV-only AM path never emitted). One row per ClinVar variant, with
    #    coords + clnsig/stars; AM fields blank; ACMG flagged at gene level.
    for vkey, cv in clinvar_by_variant.items():
        if vkey in seen_keys:
            continue
        sample, gene_u, chr_, pos_, ref_, alt_ = vkey
        in_acmg = sample in acmg_carrier_by_gene.get(gene_u, set())
        blanks = {h: "" for h in a_header}
        blanks.update({"sample_id": sample, "cohort": args.cohort, "gene": gene_u,
                       "chr": chr_, "pos": pos_, "ref": ref_, "alt": alt_,
                       "passes_QC": "TRUE", "is_AM_carrier_primary": "FALSE",
                       "would_be_carrier_domain_aggregate": "FALSE"})
        b_rows.append([blanks.get(h, "") for h in a_header] + [
            "yes" if in_acmg else "no", "no", "yes",
            cv["clnsig"], cv["review_stars"], _category(in_acmg, False, True),
        ])
        seen_keys.add(vkey)

    # 3) ACMG-only (sample,gene) pairs with NO AM/ClinVar variant on file.
    for gene_u, samples in acmg_carrier_by_gene.items():
        for sample in samples:
            if any(k[0] == sample and k[1] == gene_u for k in seen_keys):
                continue
            blanks = {h: "" for h in a_header}
            blanks.update({"sample_id": sample, "cohort": args.cohort, "gene": gene_u,
                           "passes_QC": "TRUE", "is_AM_carrier_primary": "FALSE",
                           "would_be_carrier_domain_aggregate": "FALSE"})
            b_rows.append([blanks.get(h, "") for h in a_header] + [
                "yes", "no", "no", "", "", "ACMG_PLP",
            ])

    b_path = os.path.join(out_dir, "B_ACMG_vs_AM_comparison.tsv")
    write_tsv(b_path, b_header, b_rows)

    # ---- B_summary_by_gene — per-gene carrier-row counts per track ---------
    gene_i = a_header.index("gene")
    n_h = len(a_header)  # in_ACMG_PLP, in_AM_primary, in_ClinVar_PLP are the next 3 cols
    counts: Dict[str, Counter] = defaultdict(Counter)
    for row in b_rows:
        gene_u = (row[gene_i] or "").upper()
        if row[n_h] == "yes":     counts[gene_u]["ACMG_PLP"] += 1
        if row[n_h + 1] == "yes": counts[gene_u]["AM"] += 1
        if row[n_h + 2] == "yes": counts[gene_u]["ClinVar"] += 1
        counts[gene_u]["any"] += 1
    bs_header = ["gene", "n_ACMG_PLP", "n_AM", "n_ClinVar_PLP", "n_any"]
    bs_rows = []
    for g in target_genes:
        c = counts.get(g.upper(), Counter())
        bs_rows.append([g, c.get("ACMG_PLP", 0), c.get("AM", 0),
                        c.get("ClinVar", 0), c.get("any", 0)])
    write_tsv(os.path.join(out_dir, "B_summary_by_gene.tsv"), bs_header, bs_rows)

    # ---- C_regression_matrix ----------------------------------------------
    # samples = phenotype rows (so cohort sample set is well-defined); join AM/ACMG carrier flags.
    # per-group columns: carrier_<grp>_{ACMG,AMprimary,AMonly,AM0864}
    group_names = list(groups.keys())
    c_header = ["sample_id", "cohort", "Group", "Age_at_diagnosis", "Age2", "GENDER", "genetically_determined"] \
               + pc_cols
    for grp in group_names:
        c_header += [
            f"carrier_{grp}_ACMG",
            f"carrier_{grp}_AMprimary",
            f"carrier_{grp}_AMonly",
            f"carrier_{grp}_AM0864",
        ]

    c_rows: List[List[str]] = []
    for sample, prow in pheno_by_id.items():
        age_raw = (prow.get("Age_at_diagnosis") or "").strip()
        try:
            age_f = float(age_raw)
            age2 = f"{age_f * age_f:g}"
        except ValueError:
            age2 = ""
        line = [sample, args.cohort,
                prow.get("Group", ""), age_raw, age2,
                prow.get("GENDER", ""), prow.get("genetically_determined", "")]
        # PCs — look up by the bridge image when configured, else direct.
        pc_lookup_id = sinai_to_masked.get(sample, "") if bridge_used else sample
        prec = pc_by_id.get(pc_lookup_id, {})
        for pc in pc_cols:
            line.append(prec.get(pc, ""))
        # group carriers
        for grp in group_names:
            grp_genes = [g.upper() for g in groups[grp]]
            acmg_hit = any(sample in acmg_carrier_by_gene.get(g, set()) for g in grp_genes)
            am_hit   = any(sample in am_primary_by_gene.get(g, set())   for g in grp_genes)
            am_only  = am_hit and not acmg_hit
            am0864   = any(sample in am_global_0864_by_gene.get(g, set()) for g in grp_genes)
            line += ["TRUE" if acmg_hit else "FALSE",
                     "TRUE" if am_hit   else "FALSE",
                     "TRUE" if am_only  else "FALSE",
                     "TRUE" if am0864   else "FALSE"]
        c_rows.append(line)

    write_tsv(os.path.join(out_dir, "C_regression_matrix.tsv"), c_header, c_rows)

    # ---- per-cohort summary fragment --------------------------------------
    # Aggregated into results/SUMMARY.md by the LSF 04 driver.
    frag = os.path.join(out_dir, "SUMMARY_fragment.md")
    with open(frag, "w") as fh:
        fh.write(f"## Cohort: {args.cohort}\n\n")
        fh.write(f"- VCF↔phenotype match rate ({id_col}): **{rep.match_rate:.3f}**  "
                 f"(matched={rep.matched}/{rep.left_n}; phenotype rows={rep.right_n})\n")
        if pc_id_col:
            fh.write(f"- PC id column: `{pc_id_col}`"
                     f"{'  (via id_bridge: SINAI_ID -> MASKED_MRN)' if bridge_used else ''}\n")
            if rep_pc is not None:
                fh.write(f"- phenotype↔PC match rate: **{rep_pc.match_rate:.3f}**  "
                         f"(matched={rep_pc.matched}/{rep_pc.left_n})\n")
        fh.write(f"- ACMG P/LP rows kept: **{acmg_rows_kept}**; PTV-only excluded: **{acmg_rows_excluded_ptv}**\n")
        fh.write(f"- ClinVar P/LP ≥2★ carrier rows: **{clinvar_rows_kept}** "
                 f"(all variant types incl. indels; normalized — see refs/REFERENCE_REPORT.md)\n")
        fh.write(f"- Table A rows (per sample × variant): **{len(a_rows)}**\n")
        fh.write(f"- Table B rows (union ACMG/AM/ClinVar carriers, all variant types): **{len(b_rows)}**\n")
        fh.write(f"- Table C rows (per sample, regression-ready): **{len(c_rows)}**\n\n")
        fh.write(f"### Carrier counts by gene group ({args.cohort})\n\n")
        fh.write("| group | ACMG | ClinVar | AM-primary | AM-only (excl. ACMG+ClinVar) |\n|---|---|---|---|---|\n")
        for grp in group_names:
            grp_genes = [g.upper() for g in groups[grp]]
            n_acmg = sum(1 for s in pheno_by_id if any(s in acmg_carrier_by_gene.get(g, set()) for g in grp_genes))
            n_cv   = sum(1 for s in pheno_by_id if any(s in clinvar_carrier_by_gene.get(g, set()) for g in grp_genes))
            n_am   = sum(1 for s in pheno_by_id if any(s in am_primary_by_gene.get(g, set())   for g in grp_genes))
            n_only = sum(1 for s in pheno_by_id if any(s in am_primary_by_gene.get(g, set())   for g in grp_genes)
                         and not any(s in acmg_carrier_by_gene.get(g, set())    for g in grp_genes)
                         and not any(s in clinvar_carrier_by_gene.get(g, set()) for g in grp_genes))
            fh.write(f"| {grp} | {n_acmg} | {n_cv} | {n_am} | {n_only} |\n")
        if rep.left_only_samples:
            fh.write(f"\n*VCF IDs not in phenotype (up to 5):* `{rep.left_only_samples[:5]}`\n")
        if rep.right_only_samples:
            fh.write(f"*phenotype IDs not in VCF (up to 5):* `{rep.right_only_samples[:5]}`\n")
    LOG.info("wrote per-cohort summary fragment -> %s", frag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
