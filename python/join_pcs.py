"""Standalone PC-join helper.

Joins PC1..PC10 into an existing C_regression_matrix.tsv. Runs in seconds
without LSF — call it after step 04 has produced the matrix.

The PC file is PLINK-style (whitespace-separated, columns ID1, ID2, PC1...).
The matrix and bridge are tab-separated. The bridge maps SINAI_ID -> MASKED_MRN
so cohort I can reach the PC file (which keys on MASKED_MRN); cohort II's
matrix already carries MASKED_MRN.

Usage:
    python3 python/join_pcs.py --cohort cohortI
    python3 python/join_pcs.py --cohort cohortII
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import pandas as pd
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent


def log(msg: str) -> None:
    print(f"[join_pcs] {msg}", file=sys.stderr, flush=True)


def load_cfg(path: Path) -> dict:
    with open(path) as f:
        return yaml.safe_load(f)


def pick_pc_id_col(pcs: pd.DataFrame, mrn_values: set[str]) -> str:
    """Find the PC-file column whose values overlap most with MRN values.

    PLINK convention is FID + IID; ID2 (IID) usually carries the sample id.
    We pick by overlap instead of assuming, so the helper is robust to either
    layout.
    """
    candidates = [c for c in pcs.columns if not c.upper().startswith("PC")]
    best_col, best_overlap = None, 0
    for c in candidates:
        vals = set(pcs[c].dropna().astype(str).str.strip().head(2000))
        overlap = len(vals & mrn_values)
        log(f"  candidate PC-id col {c!r}: overlap={overlap} (sample of 2000)")
        if overlap > best_overlap:
            best_overlap, best_col = overlap, c
    if not best_col or best_overlap == 0:
        raise SystemExit("No PC-file column overlaps the bridge MASKED_MRN values.")
    log(f"  picked {best_col!r} (overlap={best_overlap})")
    return best_col


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default=str(REPO_ROOT / "config" / "config.yaml"))
    ap.add_argument("--cohort", required=True, choices=["cohortI", "cohortII"])
    ap.add_argument("--matrix", default=None,
                    help="Override matrix path (default: results/<cohort>/C_regression_matrix.tsv)")
    ap.add_argument("--dry-run", action="store_true",
                    help="Report match rates without writing the matrix back.")
    args = ap.parse_args()

    cfg = load_cfg(Path(args.config))

    results_dir = Path(cfg["paths"]["results_dir"]) / args.cohort
    matrix_path = Path(args.matrix) if args.matrix else results_dir / "C_regression_matrix.tsv"
    pc_path = Path(cfg["pcs"]["pcs_tsv"])
    pc_cols_wanted = list(cfg["pcs"].get("pc_columns") or ["PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10"])

    bridge_cfg = cfg.get("id_bridge") or {}
    use_bridge = bool(bridge_cfg.get("use_for", {}).get(args.cohort, False))
    bridge_path = Path(bridge_cfg.get("masked_mrn_map", "")) if bridge_cfg else None
    sinai_col = bridge_cfg.get("map_sinai_id_col", "RGNID")
    mrn_col = bridge_cfg.get("map_masked_mrn_col", "MASKED_MRN")

    log(f"matrix: {matrix_path}")
    log(f"PCs:    {pc_path}")
    log(f"bridge: {bridge_path if use_bridge else '(not used for ' + args.cohort + ')'}")

    if not matrix_path.is_file():
        raise SystemExit(f"matrix not found: {matrix_path}")
    if not pc_path.is_file():
        raise SystemExit(f"PC file not found: {pc_path}")

    # ---- read inputs --------------------------------------------------------
    mat = pd.read_csv(matrix_path, sep="\t", dtype=str)
    log(f"matrix: {len(mat):,} rows, {len(mat.columns)} columns")

    # PC file is PLINK-style whitespace-separated
    pcs = pd.read_csv(pc_path, sep=r"\s+", engine="python", dtype=str)
    log(f"PCs:    {len(pcs):,} rows, columns = {list(pcs.columns)}")

    pc_cols_present = [c for c in pc_cols_wanted if c in pcs.columns]
    if not pc_cols_present:
        raise SystemExit(f"None of {pc_cols_wanted} found in PC file columns {list(pcs.columns)}")
    log(f"PC columns kept: {pc_cols_present}")

    # ---- ensure MASKED_MRN is on the matrix ---------------------------------
    # Matrix's per-sample identifier lives in `sample_id` (set by step 04).
    # For cohort I, values look like `SINAI_<digits>_<alphanum>` and need to be
    # bridged through Masked_mrn_map.txt to reach MASKED_MRN.
    # For cohort II, `sample_id` IS MASKED_MRN already (no bridge needed).
    MATRIX_ID = "sample_id"
    if MATRIX_ID not in mat.columns:
        raise SystemExit(
            f"matrix is missing the {MATRIX_ID!r} column. "
            f"Available columns: {list(mat.columns)[:20]}..."
        )

    if "MASKED_MRN" not in mat.columns:
        if not use_bridge:
            # Cohort II: sample_id is already MASKED_MRN. Promote it.
            mat["MASKED_MRN"] = mat[MATRIX_ID].astype(str).str.strip()
            log(f"no bridge configured for {args.cohort}; treating {MATRIX_ID} as MASKED_MRN directly")
        else:
            if not bridge_path or not bridge_path.is_file():
                raise SystemExit(f"bridge file required but not found: {bridge_path}")
            bridge = pd.read_csv(bridge_path, sep="\t", dtype=str)
            log(f"bridge: {len(bridge):,} rows, columns = {list(bridge.columns)}")
            if sinai_col not in bridge.columns or mrn_col not in bridge.columns:
                raise SystemExit(f"bridge missing {sinai_col!r} or {mrn_col!r}: {list(bridge.columns)}")
            b = bridge[[sinai_col, mrn_col]].rename(columns={sinai_col: MATRIX_ID, mrn_col: "MASKED_MRN"})
            b[MATRIX_ID] = b[MATRIX_ID].astype(str).str.strip()
            b = b.drop_duplicates(MATRIX_ID)
            mat[MATRIX_ID] = mat[MATRIX_ID].astype(str).str.strip()
            mat = mat.merge(b, on=MATRIX_ID, how="left")
            rate = mat["MASKED_MRN"].notna().mean()
            log(f"bridge join: {MATRIX_ID} -> MASKED_MRN match = {rate:.1%}")
            if rate < 0.5:
                log(f"WARNING: bridge match rate < 50% — check {MATRIX_ID} format vs {sinai_col} format")

    # ---- pick the PC-file column that joins on MASKED_MRN -------------------
    mrn_values = set(mat["MASKED_MRN"].dropna().astype(str).str.strip())
    pc_join_col = pick_pc_id_col(pcs, mrn_values)

    # ---- final merge --------------------------------------------------------
    pcs_slim = pcs[[pc_join_col] + pc_cols_present].copy()
    pcs_slim[pc_join_col] = pcs_slim[pc_join_col].astype(str).str.strip()
    pcs_slim = pcs_slim.drop_duplicates(pc_join_col)
    pcs_slim = pcs_slim.rename(columns={pc_join_col: "MASKED_MRN"})

    # drop any existing (empty) PC columns from the matrix before re-joining
    mat = mat.drop(columns=[c for c in mat.columns if c in pc_cols_present], errors="ignore")
    out = mat.merge(pcs_slim, on="MASKED_MRN", how="left")

    rate = out[pc_cols_present[0]].notna().mean()
    log(f"FINAL: {len(out):,} rows, PC match rate = {rate:.1%}")

    if args.dry_run:
        log("dry-run: not writing")
        return

    out.to_csv(matrix_path, sep="\t", index=False)
    log(f"wrote {matrix_path}")


if __name__ == "__main__":
    main()
