"""Apply the threshold policy to the gathered annotated table.

Primary: gene-specific calibrated threshold ONLY (locked decision A).
Domain-aggregate: RECORDED but NEVER used for primary call.
Uncovered: retained with am_pathogenicity, primary=FALSE.

Input : intermediate/<cohort>/am_annotated.tsv
        (header from python/annotate_am.py; one row per (sample, variant))
Output: intermediate/<cohort>/A_variant_level_per_person.tsv
        (super-set of the input + threshold columns; passes_QC is TRUE for all
         rows here since step 01 already masked failing GTs and step 02 only
         emitted non-missing/non-homref samples)
"""
from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from util import LOG, die, load_config, read_tsv_dicts, resolve, write_tsv  # noqa: E402


def load_thresholds(refs_dir: str) -> Dict[str, Tuple[Optional[float], Optional[float]]]:
    path = os.path.join(refs_dir, "calibration_thresholds.tsv")
    if not os.path.isfile(path):
        die(f"calibration table missing: {path} (run scripts/00_prepare_refs.sh)")
    _, rows = read_tsv_dicts(path)
    out: Dict[str, Tuple[Optional[float], Optional[float]]] = {}
    def _f(v: str) -> Optional[float]:
        v = (v or "").strip()
        if not v or v.upper() == "NA":
            return None
        try:
            return float(v)
        except ValueError:
            return None
    for r in rows:
        out[r["gene"].upper()] = (_f(r["gene_specific_threshold"]),
                                  _f(r["domain_aggregate_threshold"]))
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=None)
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--in", dest="inp", required=True, help="am_annotated.tsv")
    ap.add_argument("--out", required=True, help="A_variant_level_per_person.tsv")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    refs_dir = resolve(cfg, cfg["paths"]["refs_dir"])
    thresholds = load_thresholds(refs_dir)
    policy = cfg.get("calibration", {}).get("primary_threshold", "gene_specific")
    if policy != "gene_specific":
        die(f"calibration.primary_threshold must be 'gene_specific' (got {policy!r}); the spec forbids other modes")

    id_col = cfg["cohorts"][args.cohort]["id_column"]

    _, rows = read_tsv_dicts(args.inp)
    out_header = [
        "sample_id", "id_column_used", "cohort", "chr", "pos", "ref", "alt",
        "gene", "transcript", "protein_variant",
        "am_pathogenicity", "gene_specific_threshold", "domain_aggregate_threshold",
        "threshold_source", "genotype", "DP", "GQ", "passes_QC",
        "is_AM_carrier_primary", "would_be_carrier_domain_aggregate",
    ]
    out_rows: List[List[str]] = []
    n_primary = 0
    n_domain_only = 0
    n_uncov = 0
    for r in rows:
        gene = (r.get("gene") or "").upper()
        try:
            amp = float(r["am_pathogenicity"])
        except (KeyError, ValueError):
            amp = float("nan")
        gs, da = thresholds.get(gene, (None, None))
        if gs is not None:
            source = "gene_specific"; n_primary += 1
            is_primary = amp == amp and amp >= gs  # NaN-safe
            would_da = (da is not None) and (amp == amp and amp >= da)
        elif da is not None:
            source = "domain_aggregate_recorded_only"; n_domain_only += 1
            is_primary = False
            would_da = amp == amp and amp >= da
        else:
            source = "uncovered"; n_uncov += 1
            is_primary = False
            would_da = False
        out_rows.append([
            r.get("sample_id", ""), id_col, args.cohort,
            r.get("chr", ""), r.get("pos", ""), r.get("ref", ""), r.get("alt", ""),
            r.get("gene", ""), r.get("transcript", ""), r.get("protein_variant", ""),
            r.get("am_pathogenicity", ""),
            "NA" if gs is None else f"{gs:g}",
            "NA" if da is None else f"{da:g}",
            source,
            r.get("genotype", ""), r.get("DP", ""), r.get("GQ", ""),
            "TRUE",  # see module docstring
            "TRUE" if is_primary else "FALSE",
            "TRUE" if would_da else "FALSE",
        ])

    write_tsv(args.out, out_header, out_rows)
    LOG.info("%s call_carriers: rows=%d  bygene_class: primary=%d domain_only_recorded=%d uncovered=%d",
             args.cohort, len(out_rows), n_primary, n_domain_only, n_uncov)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
