"""Build the per-pipeline reference artifacts from the public sources.

Driven by `scripts/00_prepare_refs.sh`. Internet is required to fetch the
Zenodo calibration record and the gencode v32 GTF; on failure we emit a
PLACEHOLDER calibration with all-NA thresholds (and the AM-row min/max BED
fallback for coordinates), and warn loudly in the REFERENCE_REPORT.

Outputs (under refs/):
- chen_calibration.target_genes.tsv.gz   bgzip+tabix per-variant Chen labels (28-gene subset)
- chen_summary_by_gene.tsv               per-gene Chen coverage stats
- gene_transcript_map.tsv                gene_name | transcript_id | matched_by | am_transcript_used
- target_genes.exons.bed                 chrom\tstart\tend  (canonical-transcript exons, merged)
- AlphaMissense_hg38.subset_targets.tsv.gz   AM subset to target-gene transcripts; tabixed
- REFERENCE_REPORT.md                    per-gene Chen coverage table (hard gate)
"""
from __future__ import annotations

import argparse
import gzip
import os
import shutil
import subprocess
import sys
from typing import Dict, List, Optional, Tuple

import requests  # type: ignore[import-untyped]

# allow `python python/prepare_refs.py` from repo root
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from util import (  # noqa: E402
    LOG,
    chrom_norm,
    die,
    load_config,
    open_text,
    parse_gtf_attrs,
    resolve,
    stream_tsv_dicts,
    strip_version,
    write_tsv,
)


# ----------------------------------------------------------------------------
# small parsing helpers (used by subset_chen_calibration and subset_am)
# ----------------------------------------------------------------------------
def _open_table_text(path: str):
    """Open a TSV/CSV that may be plain or .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _find_col(lc: List[str], *keywords: str) -> Optional[int]:
    for i, h in enumerate(lc):
        if all(kw in h for kw in keywords):
            return i
    return None



# ----------------------------------------------------------------------------
# gencode parsing
# ----------------------------------------------------------------------------
def parse_gencode_map(path: str, target_genes: List[str]) -> List[Tuple[str, str]]:
    """Parse gencode.v32.transcriptID_genename.tsv -> deduped (gene, transcript_id) pairs
    restricted to target genes. Drops rows with empty transcript_id (the first
    row per gene in this file is the gene-level row).
    """
    tg = {g.upper() for g in target_genes}
    pairs: set[Tuple[str, str]] = set()
    n_raw = 0
    n_empty = 0
    for row in stream_tsv_dicts(path):
        n_raw += 1
        g = (row.get("gene_name") or "").strip().upper()
        t = (row.get("transcript_id") or "").strip()
        if g not in tg:
            continue
        if not t:
            n_empty += 1
            continue
        pairs.add((g, t))
    LOG.info("gencode map: %d raw rows; %d empty-transcript rows dropped; %d unique (gene,transcript) for the 28 panel", n_raw, n_empty, len(pairs))
    missing = sorted(tg - {g for (g, _) in pairs})
    if missing:
        LOG.warning("gencode map missing %d target gene(s): %s", len(missing), ", ".join(missing))
    return sorted(pairs)


def fetch_gencode_gtf(cfg: dict, refs_dir: str) -> Optional[str]:
    """Return path to a local copy of the gencode v32 GTF, or None on failure."""
    refs = cfg.get("references", {})
    local = (refs.get("gencode_gtf_local") or "").strip()
    if local and os.path.isfile(local):
        LOG.info("using local gencode GTF: %s", local)
        return local
    url = refs.get("gencode_gtf_url")
    if not url:
        return None
    dst = os.path.join(refs_dir, "gencode.v32.primary_assembly.annotation.gtf.gz")
    if os.path.isfile(dst):
        LOG.info("gencode GTF already present: %s", dst)
        return dst
    try:
        LOG.info("downloading gencode v32 GTF: %s", url)
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with open(dst, "wb") as fh:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    fh.write(chunk)
        LOG.info("downloaded -> %s (%.1f MB)", dst, os.path.getsize(dst) / 1e6)
        return dst
    except Exception as e:  # noqa: BLE001
        LOG.warning("gencode GTF download failed: %s; will fall back to AM-row min/max", e)
        return None


def exon_bed_from_gtf(gtf_path: str, transcript_to_gene: Dict[str, str]) -> Dict[str, List[Tuple[str, int, int]]]:
    """Return {gene -> [(chrom,start0,end), ...]} of exon intervals for the
    target transcripts. Coordinates are converted from GTF 1-based inclusive
    to BED 0-based half-open.

    transcript_to_gene maps versioned AND unversioned transcript ids.
    """
    out: Dict[str, List[Tuple[str, int, int]]] = {}
    with open_text(gtf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "exon":
                continue
            attrs = parse_gtf_attrs(parts[8])
            tid = attrs.get("transcript_id", "")
            gene = transcript_to_gene.get(tid) or transcript_to_gene.get(strip_version(tid))
            if not gene:
                continue
            try:
                start1 = int(parts[3])
                end1 = int(parts[4])
            except ValueError:
                continue
            out.setdefault(gene, []).append((chrom_norm(parts[0]), start1 - 1, end1))
    return out


def merge_intervals(ivs: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
    if not ivs:
        return []
    ivs = sorted(ivs)
    merged: List[Tuple[str, int, int]] = [ivs[0]]
    for c, s, e in ivs[1:]:
        pc, ps, pe = merged[-1]
        if c == pc and s <= pe:
            merged[-1] = (pc, ps, max(pe, e))
        else:
            merged.append((c, s, e))
    return merged


# ----------------------------------------------------------------------------
# Chen/Pejaver calibration table (variant-level) — subset + tabix
# ----------------------------------------------------------------------------
# Per Chen/Pejaver methodology, each variant in the published table is
# pre-assigned an ACMG-PP3/BP4 evidence label using the gene's (or its PFAM
# domain's) calibrated posterior curve. Downstream we apply by PER-VARIANT
# LOOKUP — we do NOT re-derive per-gene thresholds. See README + Tavtigian
# Bayesian framework (calib_step03.py in the Chen repo).

# Evidence labels that count as a "carrier" at a given minimum strength.
# `+` suffix is the Tavtigian intermediate strength (+3 points, between
# Moderate=+2 and Strong=+4). BP4_* are anti-pathogenic and never count.
_PP3_AT_LEAST_MODERATE = {
    "pp3_moderate", "pp3_moderate+",
    "pp3_strong",   "pp3_strong+",
    "pp3_verystrong", "pp3_very_strong",
}
_PP3_AT_LEAST_SUPPORTING = _PP3_AT_LEAST_MODERATE | {"pp3_supporting", "pp3_supporting+"}
_PP3_AT_LEAST_STRONG = {"pp3_strong", "pp3_strong+", "pp3_verystrong", "pp3_very_strong"}
_PP3_AT_LEAST_VERYSTRONG = {"pp3_verystrong", "pp3_very_strong"}


def evidence_set_for(min_strength: str) -> set:
    s = (min_strength or "").strip().lower().replace(" ", "")
    if s in ("supporting", "pp3supporting", "pp3_supporting"):
        return _PP3_AT_LEAST_SUPPORTING
    if s in ("strong", "pp3strong", "pp3_strong"):
        return _PP3_AT_LEAST_STRONG
    if s in ("verystrong", "very_strong", "pp3verystrong", "pp3_verystrong"):
        return _PP3_AT_LEAST_VERYSTRONG
    return _PP3_AT_LEAST_MODERATE  # default


def normalize_evidence(label: str) -> str:
    """Normalize an evidence label for set-membership comparison."""
    if not label:
        return ""
    return label.strip().lower().replace(" ", "").replace("-", "_")


def subset_chen_calibration(chen_csv: str, target_genes: List[str], out_path: str,
                            min_strength: str) -> Dict[str, dict]:
    """Stream the Chen variant-level calibration CSV, keep rows for target
    genes only, write a sorted bgzip+tabix indexed TSV.

    Output TSV columns (header "#"-prefixed for tabix):
      #chrom  pos  ref  alt  gene_symbol  evidence  points  calibration_approach  domain  vep_score

    Returns a per-gene coverage summary:
      { gene_upper: { n_in_chen, n_pp3_atleast_min, n_single_gene, n_domain_agg,
                      domains_seen (set), example_evidence (mode) } }
    """
    import csv as _csv
    refs_dir = os.path.dirname(out_path) or "."
    os.makedirs(refs_dir, exist_ok=True)
    tg_set = {g.upper() for g in target_genes}
    pos_set = evidence_set_for(min_strength)

    unsorted = out_path[:-3] + ".unsorted.tsv" if out_path.endswith(".gz") else out_path + ".unsorted.tsv"
    sorted_tmp = out_path[:-3] if out_path.endswith(".gz") else out_path + ".tmp"

    summary: Dict[str, dict] = {
        g.upper(): {"n_in_chen": 0, "n_pp3_atleast_min": 0,
                    "n_single_gene": 0, "n_domain_agg": 0,
                    "domains_seen": set(), "labels": Counter_dict()}
        for g in target_genes
    }

    n_total = 0
    n_kept = 0
    LOG.info("subsetting Chen calibration table → %d target genes from %s", len(tg_set), chen_csv)
    with _open_table_text(chen_csv) as fh, open(unsorted, "wt") as out_fh:
        rdr = _csv.reader(fh)
        header = next(rdr)
        # header has been seen (e.g. gene_symbol, ...). Locate columns case-insensitively.
        lc = [h.strip().lstrip("#").strip().lower() for h in header]
        gene_i  = _find_col(lc, "gene", "symbol") or _find_col(lc, "gene_name") or _find_col(lc, "gene")
        chrom_i = _find_col(lc, "chrom") if _find_col(lc, "chrom") is not None else _find_col(lc, "chr")
        pos_i   = _find_col(lc, "pos")
        ref_i   = _find_col(lc, "ref")
        alt_i   = _find_col(lc, "alt")
        ev_i    = _find_col(lc, "evidence")
        pts_i   = _find_col(lc, "points")
        ca_i    = _find_col(lc, "calibration", "approach") or _find_col(lc, "approach")
        dom_i   = _find_col(lc, "domain")
        score_i = _find_col(lc, "vep_score") or _find_col(lc, "vep", "score") or _find_col(lc, "score")
        for col, name in [(gene_i, "gene_symbol"), (chrom_i, "chrom"), (pos_i, "pos"),
                          (ref_i, "ref"), (alt_i, "alt"), (ev_i, "evidence"),
                          (ca_i, "calibration_approach")]:
            if col is None:
                die(f"Chen CSV missing required column: {name}; got header {header}")
        out_fh.write("#" + "\t".join(["chrom", "pos", "ref", "alt", "gene_symbol",
                                       "evidence", "points", "calibration_approach",
                                       "domain", "vep_score"]) + "\n")
        for r in rdr:
            n_total += 1
            try:
                g = r[gene_i].strip().upper()
            except IndexError:
                continue
            if g not in tg_set:
                continue
            chrom = chrom_norm(r[chrom_i].strip())
            pos = r[pos_i].strip()
            ref = r[ref_i].strip()
            alt = r[alt_i].strip()
            ev  = r[ev_i].strip()
            pts = r[pts_i].strip() if pts_i is not None and pts_i < len(r) else ""
            ca  = r[ca_i].strip()
            dom = r[dom_i].strip() if dom_i is not None and dom_i < len(r) else ""
            vep = r[score_i].strip() if score_i is not None and score_i < len(r) else ""
            out_fh.write("\t".join([chrom, pos, ref, alt, g, ev, pts, ca, dom, vep]) + "\n")
            n_kept += 1
            s = summary[g]
            s["n_in_chen"] += 1
            ca_l = ca.lower()
            if ca_l.startswith("single") or ca_l == "gene_specific":
                s["n_single_gene"] += 1
            elif ca_l.startswith("domain") or ca_l.startswith("agg"):
                s["n_domain_agg"] += 1
            if dom and dom.lower() not in ("no_pfam", "no-pfam", "none", ""):
                s["domains_seen"].add(dom)
            ev_n = normalize_evidence(ev)
            s["labels"][ev_n] += 1
            if ev_n in pos_set:
                s["n_pp3_atleast_min"] += 1

    LOG.info("Chen subset: %d rows scanned, %d kept (%.2f%%); now sorting by (chrom,pos)",
             n_total, n_kept, 100.0 * n_kept / max(n_total, 1))

    # Sort by chrom (lex) then pos (numeric) — keeps the "#"-prefixed header on top.
    # We use the shell `sort` to avoid loading the file in memory.
    subprocess.check_call(
        ["bash", "-c",
         f"(head -n 1 {shellquote(unsorted)} && tail -n +2 {shellquote(unsorted)} "
         f"| sort -t$'\\t' -k1,1 -k2,2n) > {shellquote(sorted_tmp)}"]
    )
    os.remove(unsorted)
    subprocess.check_call(["bgzip", "-f", sorted_tmp])
    if sorted_tmp + ".gz" != out_path:
        shutil.move(sorted_tmp + ".gz", out_path)
    subprocess.check_call(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", out_path])
    LOG.info("Chen calibration target subset -> %s (+ .tbi)", out_path)
    return summary


# ----------------------------------------------------------------------------
# ClinVar subset (target genes, P/LP >= 2 stars -> tabixed)
# ----------------------------------------------------------------------------
# NCBI CLNREVSTAT -> gold-star mapping.
_CLNREVSTAT_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
    "no_classification_for_the_single_variant": 0,
    "no_interpretation_for_the_single_variant": 0,
}
# Molecular-consequence (MC) SO terms treated as PTV/LoF (dropped when
# clinvar.exclude_ptv is true, to stay comparable to the ACMG comparator).
_CLINVAR_PTV_MC = {
    "nonsense", "frameshift_variant", "stop_gained", "stop_lost", "start_lost",
    "splice_acceptor_variant", "splice_donor_variant",
}


def _clnrevstat_stars(clnrevstat: str) -> int:
    key = (clnrevstat or "").strip().lower()
    if key in _CLNREVSTAT_STARS:
        return _CLNREVSTAT_STARS[key]
    # Defensive: any "conflicting" phrasing not enumerated above is 1 star.
    # Match "conflicting" (not "conflict") so "no_conflicts" is unaffected.
    if "conflicting" in key:
        return 1
    if key.startswith("criteria_provided"):
        return 1
    return 0


def _parse_info(info: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for field in info.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        elif field:
            out[field] = ""
    return out


def subset_clinvar(clinvar_vcf: str, target_genes: List[str], out_path: str,
                   clinvar_cfg: dict) -> Dict[str, dict]:
    """Stream the ClinVar VCF (GRCh38), keep P/LP variants in target genes that
    pass the review-star / conflicting / PTV filters, and write a sorted
    bgzip+tabix indexed TSV.

    Output TSV columns (header "#"-prefixed for tabix):
      #chrom  pos  ref  alt  gene_symbol  clnsig  review_stars  mc

    Returns per-gene coverage: { gene_upper: {"n_clinvar_plp_2star": int} }.
    """
    refs_dir = os.path.dirname(out_path) or "."
    os.makedirs(refs_dir, exist_ok=True)
    tg_set = {g.upper() for g in target_genes}
    min_stars = int(clinvar_cfg.get("min_review_stars", 2))
    sig_include = {s.strip().lower() for s in clinvar_cfg.get(
        "sig_include", ["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"])}
    exclude_conflicting = bool(clinvar_cfg.get("exclude_conflicting", True))
    exclude_ptv = bool(clinvar_cfg.get("exclude_ptv", True))

    summary: Dict[str, dict] = {g.upper(): {"n_clinvar_plp_2star": 0} for g in target_genes}

    unsorted = out_path[:-3] + ".unsorted.tsv" if out_path.endswith(".gz") else out_path + ".unsorted.tsv"
    sorted_tmp = out_path[:-3] if out_path.endswith(".gz") else out_path + ".tmp"

    n_total = 0
    n_kept = 0
    LOG.info("subsetting ClinVar VCF → %d target genes (min_stars=%d) from %s",
             len(tg_set), min_stars, clinvar_vcf)
    with _open_table_text(clinvar_vcf) as fh, open(unsorted, "wt") as out_fh:
        out_fh.write("#" + "\t".join(["chrom", "pos", "ref", "alt", "gene_symbol",
                                      "clnsig", "review_stars", "mc"]) + "\n")
        for line in fh:
            if not line or line.startswith("#"):
                continue
            n_total += 1
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            info = _parse_info(f[7])
            clnsig = info.get("CLNSIG", "")
            if clnsig.strip().lower() not in sig_include:
                continue
            clnrevstat = info.get("CLNREVSTAT", "")
            # NB: match "conflicting" (not "conflict") so the 2★ status
            # "criteria_provided,_multiple_submitters,_no_conflicts" is NOT
            # mistaken for a conflicting record.
            is_conflicting = ("conflicting" in clnrevstat.lower()
                              or "conflicting" in clnsig.lower())
            if exclude_conflicting and is_conflicting:
                continue
            stars = _clnrevstat_stars(clnrevstat)
            if stars < min_stars:
                continue
            mc = info.get("MC", "")
            if exclude_ptv:
                mc_terms = {p.split("|", 1)[1] for p in mc.split(",") if "|" in p}
                if mc_terms & _CLINVAR_PTV_MC:
                    continue
            # GENEINFO is "SYM:id" possibly "SYMA:1|SYMB:2"; keep the first
            # target-gene symbol we recognise.
            geneinfo = info.get("GENEINFO", "")
            gene = ""
            for tok in geneinfo.split("|"):
                sym = tok.split(":", 1)[0].strip().upper()
                if sym in tg_set:
                    gene = sym
                    break
            if not gene:
                continue
            chrom = chrom_norm(f[0].strip())
            pos = f[1].strip()
            ref = f[3].strip()
            alt = f[4].strip()
            if not ref or not alt or alt == ".":
                continue
            out_fh.write("\t".join([chrom, pos, ref, alt, gene, clnsig, str(stars), mc]) + "\n")
            n_kept += 1
            summary[gene]["n_clinvar_plp_2star"] += 1

    LOG.info("ClinVar subset: %d records scanned, %d kept; now sorting by (chrom,pos)",
             n_total, n_kept)
    subprocess.check_call(
        ["bash", "-c",
         f"(head -n 1 {shellquote(unsorted)} && tail -n +2 {shellquote(unsorted)} "
         f"| sort -t$'\\t' -k1,1 -k2,2n) > {shellquote(sorted_tmp)}"]
    )
    os.remove(unsorted)
    subprocess.check_call(["bgzip", "-f", sorted_tmp])
    if sorted_tmp + ".gz" != out_path:
        shutil.move(sorted_tmp + ".gz", out_path)
    subprocess.check_call(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", out_path])
    LOG.info("ClinVar P/LP ≥%d★ target subset -> %s (+ .tbi)", min_stars, out_path)
    return summary


def shellquote(p: str) -> str:
    import shlex
    return shlex.quote(p)


class Counter_dict(dict):
    """Tiny counter that returns 0 for missing keys and supports +=."""
    def __getitem__(self, k):
        return super().get(k, 0)
    def __setitem__(self, k, v):
        super().__setitem__(k, v)
    def mode(self) -> Optional[str]:
        if not self:
            return None
        return max(self.items(), key=lambda kv: kv[1])[0]


# ----------------------------------------------------------------------------
# AM subset (target transcripts -> tabixed)
# ----------------------------------------------------------------------------
def subset_am(am_path: str, target_transcripts: Dict[str, str], out_path: str) -> Tuple[int, int]:
    """Subset the AM TSV to rows whose transcript is in our target set.

    Returns (n_in, n_kept). Emits an UNgzipped intermediate then bgzips/tabixes.
    target_transcripts maps versioned and unversioned -> gene.
    """
    tmp = out_path[:-3] if out_path.endswith(".gz") else out_path + ".tmp"
    n_in = 0
    n_kept = 0
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    def _find(header_lc: Dict[str, int], *aliases: str) -> int:
        for a in aliases:
            if a.lower() in header_lc:
                return header_lc[a.lower()]
        raise KeyError(f"none of {aliases!r} in AM header")

    with open(tmp, "wt") as out_fh:
        # standard AM columns we care about. Header is "#"-prefixed so tabix
        # treats it as a comment and `tabix <file> <region>` skips it cleanly.
        out_fh.write("#" + "\t".join(["chrom", "pos", "ref", "alt", "uniprot_id", "transcript_id",
                                       "protein_variant", "am_pathogenicity", "am_class", "gene_name"]) + "\n")
        with open_text(am_path) as fh:
            header: Optional[List[str]] = None
            header_lc: Dict[str, int] = {}
            with_tid = with_chrom = with_pos = with_ref = with_alt = with_amp = -1
            with_pv = with_amc = with_upid = -1
            for line in fh:
                if not line.strip():
                    continue
                if header is None:
                    s = line.lstrip().rstrip("\n")
                    if s.startswith("#"):
                        cand = s.lstrip("#").strip()
                        if "\t" in cand:
                            header = cand.split("\t")
                    else:
                        header = s.split("\t")
                    if header is None:
                        continue
                    header_lc = {h.lower(): i for i, h in enumerate(header)}
                    with_tid   = _find(header_lc, "transcript_id")
                    with_chrom = _find(header_lc, "chrom", "#chrom")
                    with_pos   = _find(header_lc, "pos")
                    with_ref   = _find(header_lc, "ref")
                    with_alt   = _find(header_lc, "alt")
                    with_amp   = _find(header_lc, "am_pathogenicity")
                    with_pv    = header_lc.get("protein_variant", -1)
                    with_amc   = header_lc.get("am_class", -1)
                    with_upid  = header_lc.get("uniprot_id", -1)
                    continue
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                n_in += 1
                try:
                    tid = parts[with_tid]
                except IndexError:
                    continue
                gene = target_transcripts.get(tid) or target_transcripts.get(strip_version(tid))
                if not gene:
                    continue
                try:
                    chrom = chrom_norm(parts[with_chrom])
                    pos = parts[with_pos]
                    ref = parts[with_ref]
                    alt = parts[with_alt]
                    amp = parts[with_amp]
                except IndexError:
                    continue
                pv   = parts[with_pv]   if with_pv   >= 0 and with_pv   < len(parts) else ""
                amc  = parts[with_amc]  if with_amc  >= 0 and with_amc  < len(parts) else ""
                upid = parts[with_upid] if with_upid >= 0 and with_upid < len(parts) else ""
                out_fh.write("\t".join([chrom, pos, ref, alt, upid, tid, pv, amp, amc, gene]) + "\n")
                n_kept += 1

    # bgzip + tabix
    bgz = out_path
    subprocess.check_call(["bgzip", "-f", tmp])
    if tmp + ".gz" != bgz:
        shutil.move(tmp + ".gz", bgz)
    subprocess.check_call(["tabix", "-f", "-s", "1", "-b", "2", "-e", "2", bgz])
    LOG.info("AM subset: %d rows scanned, %d kept -> %s", n_in, n_kept, bgz)
    return n_in, n_kept


# ----------------------------------------------------------------------------
# ClinVar-only fast path (skip the slow AM/Chen ref rebuild)
# ----------------------------------------------------------------------------
def _resolve_clinvar_vcf(cfg: dict) -> str:
    """Env BIOAM_CLINVAR_VCF (normalized by 00) wins over config local path."""
    v = (os.environ.get("BIOAM_CLINVAR_VCF", "") or "").strip()
    if not v:
        v = (cfg.get("references", {}).get("clinvar_vcf_local") or "").strip()
    return v


def _build_clinvar_only(cfg: dict, refs_dir: str, target_genes: List[str]) -> int:
    """Build ONLY refs/clinvar_plp_2star.target_genes.tsv.gz + refs/CLINVAR_REPORT.md,
    leaving the AM/Chen/gencode/BED refs untouched. Fast: no AM-TSV streaming."""
    clinvar_cfg = cfg.get("clinvar", {}) or {}
    clinvar_subset_path = os.path.join(refs_dir, "clinvar_plp_2star.target_genes.tsv.gz")
    clinvar_vcf = _resolve_clinvar_vcf(cfg)
    if not (clinvar_vcf and os.path.isfile(clinvar_vcf)):
        die("--only-clinvar: no ClinVar VCF available (BIOAM_CLINVAR_VCF unset and "
            "references.clinvar_vcf_local empty/missing). Run 00_prepare_refs.sh so it "
            "downloads/normalizes ClinVar, or set references.clinvar_vcf_local.")
    summary = subset_clinvar(clinvar_vcf, target_genes, clinvar_subset_path, clinvar_cfg)
    n_total = sum(summary.get(g.upper(), {}).get("n_clinvar_plp_2star", 0) for g in target_genes)

    report = os.path.join(refs_dir, "CLINVAR_REPORT.md")
    with open(report, "w") as fh:
        fh.write("# BioMe — ClinVar P/LP ≥2★ subset (built with --only-clinvar)\n\n")
        fh.write(f"- ClinVar source: `{clinvar_vcf}`\n")
        fh.write(f"- Filter: CLNSIG ∈ {sorted(clinvar_cfg.get('sig_include', []))}; "
                 f"stars ≥ **{clinvar_cfg.get('min_review_stars', 2)}**; "
                 f"exclude_conflicting=**{clinvar_cfg.get('exclude_conflicting', True)}**; "
                 f"exclude_ptv=**{clinvar_cfg.get('exclude_ptv', False)}**\n")
        fh.write(f"- Total target-gene ClinVar P/LP ≥2★ variants: **{n_total}**\n\n")
        fh.write("| gene | n_ClinVar_PLP_≥2★ |\n|------|-------------------|\n")
        for g in target_genes:
            fh.write(f"| {g} | {summary.get(g.upper(), {}).get('n_clinvar_plp_2star', 0)} |\n")
    LOG.info("wrote %s", report)
    print(f"\n== ClinVar-only subset built: {n_total} target-gene P/LP ≥2★ variants -> "
          f"{clinvar_subset_path} ==")
    print(f"== review {report} ==")
    return 0


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Prepare AM/gencode/calibration references.")
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-network", action="store_true",
                    help="skip Zenodo + GTF network fetches; useful in tests")
    ap.add_argument("--only-clinvar", action="store_true",
                    help="ONLY (re)build the ClinVar P/LP subset; skip the Chen/AM/"
                         "gencode/BED work (which is slow and unchanged). Use after the "
                         "AM refs already exist and only ClinVar is new.")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    target_genes = list(cfg["target_genes"])
    target_genes_set = {g.upper() for g in target_genes}
    refs_dir = resolve(cfg, cfg["paths"]["refs_dir"])
    os.makedirs(refs_dir, exist_ok=True)

    if args.only_clinvar:
        return _build_clinvar_only(cfg, refs_dir, target_genes)

    # ---- Chen/Pejaver variant-level calibration table (per-variant lookup) --
    cal_cfg = cfg.get("calibration", {})
    min_strength = cal_cfg.get("min_evidence_strength", "Moderate")
    chen_subset_path = os.path.join(refs_dir, "chen_calibration.target_genes.tsv.gz")
    chen_csv = (os.environ.get("BIOAM_CALIBRATION_CSV", "") or "").strip()
    if not chen_csv:
        chen_csv = (cal_cfg.get("local_file") or "").strip()
    chen_summary: Dict[str, dict] = {}
    chen_source = "NOT_AVAILABLE"
    if chen_csv and os.path.isfile(chen_csv):
        chen_summary = subset_chen_calibration(chen_csv, target_genes, chen_subset_path, min_strength)
        chen_source = os.path.basename(chen_csv)
    else:
        LOG.warning("No Chen calibration CSV available (BIOAM_CALIBRATION_CSV unset and "
                    "calibration.local_file empty/missing). Step 02 will mark every "
                    "variant as not_in_Chen_table and no carriers will be called. "
                    "Run scripts/00_prepare_refs.sh on a node with internet OR drop the "
                    "file at calibration.local_file.")
        chen_summary = {g.upper(): {"n_in_chen": 0, "n_pp3_atleast_min": 0,
                                     "n_single_gene": 0, "n_domain_agg": 0,
                                     "domains_seen": set(), "labels": Counter_dict()}
                        for g in target_genes}

    # ---- ClinVar P/LP >= 2-star subset (standalone comparator) --------------
    clinvar_cfg = cfg.get("clinvar", {}) or {}
    clinvar_subset_path = os.path.join(refs_dir, "clinvar_plp_2star.target_genes.tsv.gz")
    # Env first: 00_prepare_refs.sh exports BIOAM_CLINVAR_VCF pointing at the
    # downloaded and/or left-normalized ClinVar VCF; that must win over the raw
    # config path so the subset inherits the normalization.
    clinvar_vcf = (os.environ.get("BIOAM_CLINVAR_VCF", "") or "").strip()
    if not clinvar_vcf:
        clinvar_vcf = (cfg.get("references", {}).get("clinvar_vcf_local") or "").strip()
    clinvar_summary: Dict[str, dict] = {g.upper(): {"n_clinvar_plp_2star": 0} for g in target_genes}
    clinvar_source = "NOT_AVAILABLE"
    if clinvar_vcf and os.path.isfile(clinvar_vcf):
        clinvar_summary = subset_clinvar(clinvar_vcf, target_genes, clinvar_subset_path, clinvar_cfg)
        clinvar_source = os.path.basename(clinvar_vcf)
    else:
        LOG.warning("No ClinVar VCF available (references.clinvar_vcf_local empty/missing and "
                    "BIOAM_CLINVAR_VCF unset). Step 02c will find no ClinVar carriers; the "
                    "ClinVar_PLP category will be empty. 00_prepare_refs.sh downloads it from "
                    "references.clinvar_vcf_url on a node with internet.")

    # ---- gencode map (gene -> transcripts) ----------------------------------
    gencode_map_path = resolve(cfg, cfg["references"]["gencode_transcript_map"])
    if not os.path.isfile(gencode_map_path):
        die(f"references.gencode_transcript_map not found: {gencode_map_path}")
    pairs = parse_gencode_map(gencode_map_path, target_genes)

    # ---- match each transcript to AM (exact, then version-stripped) ---------
    am_path = resolve(cfg, cfg["references"]["alphamissense_tsv"])
    if not os.path.isfile(am_path):
        die(f"references.alphamissense_tsv not found: {am_path}")

    # collect AM transcript set with a single streaming pass
    am_versioned: set[str] = set()
    am_unversioned: Dict[str, str] = {}  # stripped -> first versioned seen
    LOG.info("scanning AM file for transcript IDs (one pass): %s", am_path)
    with open_text(am_path) as fh:
        header: Optional[List[str]] = None
        tid_col: int = -1
        for line in fh:
            if not line.strip():
                continue
            if header is None:
                s = line.lstrip().rstrip("\n")
                if s.startswith("#"):
                    cand = s.lstrip("#").strip()
                    if "\t" in cand:
                        header = cand.split("\t")
                else:
                    header = s.split("\t")
                if header is None:
                    continue
                hlc = {h.lower(): i for i, h in enumerate(header)}
                if "transcript_id" not in hlc:
                    die(f"AM header missing transcript_id column: {header}")
                tid_col = hlc["transcript_id"]
                continue
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            try:
                tid = parts[tid_col]
            except IndexError:
                continue
            if tid not in am_versioned:
                am_versioned.add(tid)
                stripped = strip_version(tid)
                am_unversioned.setdefault(stripped, tid)
    LOG.info("AM has %d distinct transcript_ids", len(am_versioned))

    map_rows: List[List[str]] = []
    transcript_to_gene: Dict[str, str] = {}
    coverage: Dict[str, dict] = {}
    for gene, tid in pairs:
        if tid in am_versioned:
            map_rows.append([gene, tid, "exact", tid])
            transcript_to_gene[tid] = gene
            transcript_to_gene[strip_version(tid)] = gene
            coverage.setdefault(gene, {"am_found": True})
        elif strip_version(tid) in am_unversioned:
            used = am_unversioned[strip_version(tid)]
            map_rows.append([gene, tid, "version_stripped", used])
            transcript_to_gene[used] = gene
            transcript_to_gene[strip_version(used)] = gene
            coverage.setdefault(gene, {"am_found": True})
        else:
            map_rows.append([gene, tid, "not_found", ""])
    write_tsv(os.path.join(refs_dir, "gene_transcript_map.tsv"),
              ["gene_name", "transcript_id", "matched_by", "am_transcript_used"], map_rows)
    # any gene with no found transcript
    for g in target_genes:
        coverage.setdefault(g, {"am_found": False})

    # ---- AM subset to target transcripts ------------------------------------
    am_subset = os.path.join(refs_dir, "AlphaMissense_hg38.subset_targets.tsv.gz")
    subset_am(am_path, transcript_to_gene, am_subset)

    # ---- exon BED from gencode GTF (preferred); AM min/max fallback ---------
    bed_path = os.path.join(refs_dir, "target_genes.exons.bed")
    bed_source = "AM_minmax_fallback"
    coords_per_gene: Dict[str, Tuple[str, int, int]] = {}

    # --skip-network still allows a local GTF (set references.gencode_gtf_local)
    local_gtf = (cfg.get("references", {}).get("gencode_gtf_local") or "").strip()
    if args.skip_network and not (local_gtf and os.path.isfile(local_gtf)):
        gtf_path = None
    else:
        gtf_path = fetch_gencode_gtf(cfg, refs_dir)
    if gtf_path:
        try:
            bygene = exon_bed_from_gtf(gtf_path, transcript_to_gene)
            with open(bed_path, "w") as fh:
                for gene in sorted(bygene):
                    merged = merge_intervals(bygene[gene])
                    if merged:
                        c0 = merged[0][0]
                        lo = min(s for _, s, _ in merged)
                        hi = max(e for _, _, e in merged)
                        coords_per_gene[gene] = (c0, lo, hi)
                    for c, s, e in merged:
                        fh.write(f"{c}\t{s}\t{e}\t{gene}\n")
            bed_source = "gencode_exon"
            LOG.info("wrote exon BED -> %s (source=gencode_exon)", bed_path)
        except Exception as e:  # noqa: BLE001
            LOG.warning("exon BED from GTF failed (%s); falling back to AM min/max", e)
            gtf_path = None

    if not gtf_path or bed_source == "AM_minmax_fallback":
        # compute min/max POS per gene from the AM subset
        per_gene: Dict[str, List[Tuple[str, int]]] = {}
        import gzip as _gz
        with _gz.open(am_subset, "rt") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ci = hdr.index("chrom"); pi = hdr.index("pos"); gi = hdr.index("gene_name")
            for line in fh:
                p = line.rstrip("\n").split("\t")
                try:
                    per_gene.setdefault(p[gi], []).append((p[ci], int(p[pi])))
                except (ValueError, IndexError):
                    continue
        with open(bed_path, "w") as fh:
            for gene in sorted(per_gene):
                xs = per_gene[gene]
                if not xs:
                    continue
                c = xs[0][0]
                lo = min(p for _, p in xs) - 1
                hi = max(p for _, p in xs)
                coords_per_gene[gene] = (c, lo, hi)
                fh.write(f"{c}\t{lo}\t{hi}\t{gene}\n")
        if bed_source == "AM_minmax_fallback":
            LOG.warning("wrote BED via AM-row min/max FALLBACK -> %s. This includes introns; reviewer requested exon-based BED. Set references.gencode_gtf_local or re-run on a login node with internet.", bed_path)

    # ---- chen_summary_by_gene.tsv (per-gene Chen coverage stats) ------------
    summary_rows: List[List[str]] = []
    for g in target_genes:
        s = chen_summary.get(g.upper(), {"n_in_chen": 0, "n_pp3_atleast_min": 0,
                                          "n_single_gene": 0, "n_domain_agg": 0,
                                          "domains_seen": set(), "labels": Counter_dict()})
        ca_mode = "single_gene" if s["n_single_gene"] >= s["n_domain_agg"] else "domain_aggregate"
        if s["n_in_chen"] == 0:
            ca_mode = "none"
        summary_rows.append([
            g,
            str(s["n_in_chen"]),
            str(s["n_pp3_atleast_min"]),
            str(s["n_single_gene"]),
            str(s["n_domain_agg"]),
            ca_mode,
            ",".join(sorted(s["domains_seen"])) or "no_pfam",
        ])
    write_tsv(os.path.join(refs_dir, "chen_summary_by_gene.tsv"),
              ["gene", "n_variants_in_chen", "n_pp3_at_least_" + min_strength.lower().replace(" ", ""),
               "n_calibration_single_gene", "n_calibration_domain_aggregate",
               "calibration_approach_majority", "pfam_domains"], summary_rows)

    # ---- REFERENCE_REPORT.md (per-gene Chen coverage table + hard gate) -----
    n_with_pp3 = 0
    n_only_pp3_sup = 0
    n_no_chen = 0
    lines: List[str] = []
    lines.append("# BioMe AlphaMissense — Reference Report\n")
    lines.append(f"- Chen calibration source: **{chen_source}**")
    lines.append(f"- Min carrier evidence: **PP3_{min_strength}** (variant counts as carrier "
                 f"iff its Chen `evidence` label ≥ this)")
    lines.append(f"- AM file: `{am_path}`")
    lines.append(f"- Gencode map: `{gencode_map_path}`")
    lines.append(f"- Exon-BED source: **{bed_source}**" +
                 ("  ⚠️  *AM-row min/max fallback — drop a real GTF at `references.gencode_gtf_local` and rerun.*"
                  if bed_source != "gencode_exon" else ""))
    lines.append("")
    lines.append("Carriers are identified per-variant by **looking up each BioMe variant in Chen's "
                 "calibration table** (joined on chrom/pos/ref/alt) and reading the published "
                 "`evidence` label. We do **NOT** re-derive per-gene thresholds — Chen has already "
                 "applied its gene-specific or domain-aggregate calibration curve to every variant "
                 "in the table.")
    lines.append("")
    lines.append("**Note on the Chen subset.** `refs/chen_calibration.target_genes.tsv.gz` "
                 "contains *every* Chen row for the 28 target genes, including BP4_* "
                 "(anti-pathogenic) and PP3_Supporting rows — the `min_evidence_strength` "
                 "filter is applied downstream by `call_carriers.py`, not in the subset. "
                 "This means downstream analyses that need the full evidence distribution "
                 "(e.g. per-gene FPR at the global 0.864 threshold, false-negative supplementary "
                 "figures) can read the subset directly without re-fetching.")
    lines.append("")
    lines.append(f"## Target panel ({len(target_genes)} genes)")
    lines.append("")
    lines.append("| gene | AM_transcript_found | coordinates_resolved | n_variants_in_Chen | "
                 f"n_PP3_≥{min_strength} | calibration_approach | PFAM_domains |")
    lines.append("|------|---------------------|----------------------|--------------------|"
                 "----------------------|----------------------|--------------|")
    for g in target_genes:
        gu = g.upper()
        s = chen_summary.get(gu, {})
        am_found = coverage.get(g, {}).get("am_found", False)
        coords = "yes" if g in coords_per_gene else "no"
        n_in = s.get("n_in_chen", 0)
        n_pp3 = s.get("n_pp3_atleast_min", 0)
        n_sg = s.get("n_single_gene", 0)
        n_da = s.get("n_domain_agg", 0)
        if n_in == 0:
            ca_label = "—"
            n_no_chen += 1
        else:
            ca_label = "single_gene" if n_sg >= n_da else "domain_aggregate"
            if n_pp3 > 0:
                n_with_pp3 += 1
            else:
                n_only_pp3_sup += 1
        domains = sorted(s.get("domains_seen", set()))
        dom_str = ",".join(domains) if domains else "(none / no_pfam)"
        if len(dom_str) > 40:
            dom_str = dom_str[:37] + "..."
        lines.append(f"| {g} | {'yes' if am_found else 'no'} | {coords} | "
                     f"{n_in} | {n_pp3} | {ca_label} | {dom_str} |")
    lines.append("")
    # ---- ClinVar P/LP >= 2-star coverage ----------------------------------
    n_clinvar_total = sum(clinvar_summary.get(g.upper(), {}).get("n_clinvar_plp_2star", 0)
                          for g in target_genes)
    lines.append("")
    lines.append("## ClinVar P/LP ≥2★ subset (standalone comparator)")
    lines.append("")
    lines.append(f"- ClinVar source: **{clinvar_source}**")
    lines.append(f"- Filter: CLNSIG ∈ {sorted(clinvar_cfg.get('sig_include', []))}; "
                 f"CLNREVSTAT ≥ **{clinvar_cfg.get('min_review_stars', 2)}** gold stars; "
                 f"exclude_conflicting=**{clinvar_cfg.get('exclude_conflicting', True)}**; "
                 f"exclude_ptv=**{clinvar_cfg.get('exclude_ptv', True)}**")
    lines.append(f"- Total target-gene ClinVar P/LP ≥2★ variants: **{n_clinvar_total}**")
    lines.append("- Carriers are called from the ALL-VARIANT QC VCF (chr<N>.qc_allvar.vcf.gz), "
                 "so non-SNV P/LP variants (indels/MNVs) are included. Both sides are "
                 "left-normalized against references.reference_fasta for reliable indel "
                 "matching (best-effort exact match if no FASTA is set).")
    lines.append("")
    lines.append("| gene | n_ClinVar_PLP_≥2★ |")
    lines.append("|------|-------------------|")
    for g in target_genes:
        n_cv = clinvar_summary.get(g.upper(), {}).get("n_clinvar_plp_2star", 0)
        lines.append(f"| {g} | {n_cv} |")
    lines.append("")
    lines.append("## Coverage summary")
    lines.append(f"- genes with ≥1 carrier-eligible Chen variant (PP3_≥{min_strength}): **{n_with_pp3}**")
    lines.append(f"- genes with Chen rows but none reaching PP3_{min_strength}: **{n_only_pp3_sup}** "
                 f"(carriers must reach this threshold to be called primary)")
    lines.append(f"- genes with no Chen rows at all: **{n_no_chen}**")
    lines.append(f"- total: **{len(target_genes)}**")
    lines.append("")
    lines.append("## How carriers are called downstream")
    lines.append("")
    lines.append(f"For every BioMe variant that passes step-01 QC, step 02 (annotate_am.py) "
                 f"does a tabix lookup against `chen_calibration.target_genes.tsv.gz` on "
                 f"`(chrom, pos, ref, alt)`. Step 03 (call_carriers.py) sets:")
    lines.append("")
    lines.append(f"- `is_AM_carrier_primary = TRUE` iff the looked-up `evidence` label is in "
                 f"{{PP3_Moderate, PP3_Moderate+, PP3_Strong, PP3_Strong+, PP3_VeryStrong}} "
                 f"(or whatever `calibration.min_evidence_strength` resolves to).")
    lines.append(f"- `chen_calibration_approach` is copied from the Chen table "
                 f"(`single_gene` / `domain_aggregate`).")
    lines.append(f"- Variants **not** in the Chen table get `is_AM_carrier_primary = FALSE` "
                 f"and `threshold_source = not_in_Chen_table`. Their AM score is still recorded.")
    lines.append("")
    lines.append("## HARD GATE")
    lines.append("Review this report before submitting steps 01+. In particular:")
    lines.append(f"- if Chen calibration source is `NOT_AVAILABLE`, drop the file at "
                 f"`calibration.local_file` (or set `BIOAM_CALIBRATION_CSV`) and re-run "
                 f"`00_prepare_refs.sh`;")
    lines.append("- if exon BED source is `AM_minmax_fallback`, set `references.gencode_gtf_local` "
                 "to a local GTF and re-run.")
    report_path = os.path.join(refs_dir, "REFERENCE_REPORT.md")
    with open(report_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    LOG.info("wrote %s", report_path)

    print()
    print(f"== Chen coverage: genes_with_PP3_>={min_strength}={n_with_pp3}  "
          f"genes_with_only_subthreshold={n_only_pp3_sup}  genes_with_no_chen={n_no_chen} ==")
    print(f"== review {report_path} before submitting 01+ ==")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
