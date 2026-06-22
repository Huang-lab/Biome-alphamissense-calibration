"""Build the per-pipeline reference artifacts from the public sources.

Driven by `scripts/00_prepare_refs.sh`. Internet is required to fetch the
Zenodo calibration record and the gencode v32 GTF; on failure we emit a
PLACEHOLDER calibration with all-NA thresholds (and the AM-row min/max BED
fallback for coordinates), and warn loudly in the REFERENCE_REPORT.

Outputs (under refs/):
- calibration_thresholds.tsv          gene | gene_specific | domain_aggregate | notes
- gene_transcript_map.tsv             gene_name | transcript_id | matched_by | am_transcript_used
- target_genes.exons.bed              chrom\tstart\tend  (canonical-transcript exons, merged)
- AlphaMissense_hg38.subset_targets.tsv.gz   AM subset to target-gene transcripts; tabixed
- REFERENCE_REPORT.md                 28-gene table + coverage counts (hard gate)
"""
from __future__ import annotations

import argparse
import gzip
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
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
# calibration fetch / parse
# ----------------------------------------------------------------------------
def _try_parse_calibration_table(content: bytes, src_label: str, target_genes: List[str]) -> Optional[Dict[str, Dict[str, str]]]:
    """Attempt to parse a calibration TSV/CSV into {gene: {'gene_specific':..., 'domain_aggregate':..., 'notes':...}}.

    Strategy: read header, detect which columns look like gene/threshold/domain;
    don't assume column order. Return None if nothing recognizable.
    """
    text = content.decode("utf-8", errors="replace")
    # auto-detect delimiter
    first = text.splitlines()[0] if text else ""
    delim = "\t" if first.count("\t") >= first.count(",") else ","
    rows = [r.split(delim) for r in text.splitlines() if r.strip()]
    if len(rows) < 2:
        return None
    header = [h.strip() for h in rows[0]]
    lc = [h.lower() for h in header]

    def find_col(*keywords: str) -> Optional[int]:
        for i, h in enumerate(lc):
            if all(kw in h for kw in keywords):
                return i
        return None

    gene_i = find_col("gene")
    gs_i = find_col("gene", "specific") or find_col("gene_thresh") or find_col("calibrated", "gene")
    da_i = find_col("domain") if find_col("domain") is not None else find_col("aggregate")
    if gene_i is None or (gs_i is None and da_i is None):
        return None

    out: Dict[str, Dict[str, str]] = {}
    tg_set = {g.upper() for g in target_genes}
    for r in rows[1:]:
        if len(r) <= gene_i:
            continue
        g = r[gene_i].strip().upper()
        if g not in tg_set:
            continue
        gs = r[gs_i].strip() if (gs_i is not None and len(r) > gs_i) else ""
        da = r[da_i].strip() if (da_i is not None and len(r) > da_i) else ""
        out[g] = {
            "gene_specific": gs if gs not in ("", "NA", "nan", "None") else "NA",
            "domain_aggregate": da if da not in ("", "NA", "nan", "None") else "NA",
            "notes": f"parsed from {src_label}",
        }
    return out or None


def fetch_calibration(cfg: dict, target_genes: List[str]) -> Tuple[Dict[str, Dict[str, str]], str]:
    """Return (gene_to_thresholds, source_label_or_warning).

    Resolves in this order:
    1) calibration.local_file   (drop-in override)
    2) calibration.zenodo_api_url   (network)
    3) placeholder (all NA) + warning
    """
    cal_cfg = cfg.get("calibration", {})
    local = cal_cfg.get("local_file") or ""
    if local:
        if not os.path.isfile(local):
            LOG.warning("calibration.local_file set to %r but file does not exist; falling back to network", local)
        else:
            LOG.info("reading calibration from local override: %s", local)
            with open(local, "rb") as fh:
                parsed = _try_parse_calibration_table(fh.read(), f"local:{local}", target_genes)
            if parsed:
                return parsed, f"local file {local}"
            LOG.warning("could not auto-parse calibration columns in %s; falling back to placeholder", local)

    url = cal_cfg.get("zenodo_api_url")
    if url:
        try:
            LOG.info("fetching calibration record from %s", url)
            resp = requests.get(url, timeout=30)
            resp.raise_for_status()
            meta = resp.json()
            files = meta.get("files", []) or []
            LOG.info("Zenodo record has %d file(s)", len(files))
            for f in files:
                name = f.get("key") or f.get("filename") or ""
                link = (f.get("links") or {}).get("self") or f.get("download") or ""
                if not link:
                    continue
                LOG.info("inspecting %s", name)
                if not (name.endswith(".tsv") or name.endswith(".csv") or name.endswith(".txt")):
                    continue
                try:
                    sub = requests.get(link, timeout=60)
                    sub.raise_for_status()
                except Exception as e:  # noqa: BLE001
                    LOG.warning("download failed for %s: %s", name, e)
                    continue
                parsed = _try_parse_calibration_table(sub.content, f"zenodo:{name}", target_genes)
                if parsed:
                    return parsed, f"zenodo:{name}"
            LOG.warning("could not parse any calibration table from Zenodo record %s", cal_cfg.get("zenodo_record"))
        except Exception as e:  # noqa: BLE001
            LOG.warning("Zenodo fetch failed: %s", e)

    # placeholder
    LOG.warning("emitting PLACEHOLDER calibration (all NA) — drop the real file at calibration.local_file and re-run")
    return ({g.upper(): {"gene_specific": "NA", "domain_aggregate": "NA",
                        "notes": "PLACEHOLDER — calibration fetch failed; replace via config.calibration.local_file"}
             for g in target_genes},
            "PLACEHOLDER (fetch failed)")


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
# main
# ----------------------------------------------------------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="Prepare AM/gencode/calibration references.")
    ap.add_argument("--config", default=None)
    ap.add_argument("--skip-network", action="store_true",
                    help="skip Zenodo + GTF network fetches; useful in tests")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    target_genes = list(cfg["target_genes"])
    target_genes_set = {g.upper() for g in target_genes}
    refs_dir = resolve(cfg, cfg["paths"]["refs_dir"])
    os.makedirs(refs_dir, exist_ok=True)

    # ---- calibration --------------------------------------------------------
    if args.skip_network and not (cfg.get("calibration", {}).get("local_file")):
        LOG.warning("--skip-network and no local_file set; emitting placeholder calibration")
        cal_map, cal_src = ({g.upper(): {"gene_specific": "NA", "domain_aggregate": "NA",
                                          "notes": "PLACEHOLDER (skip-network)"} for g in target_genes},
                            "PLACEHOLDER (skip-network)")
    else:
        cal_map, cal_src = fetch_calibration(cfg, target_genes)

    cal_out = os.path.join(refs_dir, "calibration_thresholds.tsv")
    rows = []
    for g in target_genes:
        e = cal_map.get(g.upper(), {"gene_specific": "NA", "domain_aggregate": "NA", "notes": "missing from source"})
        rows.append([g, e.get("gene_specific", "NA"), e.get("domain_aggregate", "NA"), e.get("notes", "")])
    write_tsv(cal_out, ["gene", "gene_specific_threshold", "domain_aggregate_threshold", "notes"], rows)

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

    # ---- REFERENCE_REPORT.md (28-gene table + hard gate) --------------------
    n_gs = n_da_only = n_uncov = 0
    lines: List[str] = []
    lines.append("# BioMe AlphaMissense — Reference Report\n")
    lines.append(f"- Calibration source: **{cal_src}**")
    lines.append(f"- AM file: `{am_path}`")
    lines.append(f"- Gencode map: `{gencode_map_path}`")
    lines.append(f"- Exon-BED source: **{bed_source}**" +
                 ("  ⚠️  *AM-row min/max fallback — drop a real GTF at `references.gencode_gtf_local` and rerun.*"
                  if bed_source != "gencode_exon" else ""))
    lines.append("")
    lines.append(f"## Target panel ({len(target_genes)} genes)")
    lines.append("")
    lines.append("| gene | AM_transcript_found | coordinates_resolved | gene_specific_threshold | domain_aggregate_threshold | primary_status |")
    lines.append("|------|---------------------|----------------------|-------------------------|----------------------------|----------------|")
    for g in target_genes:
        gu = g.upper()
        cal = cal_map.get(gu, {"gene_specific": "NA", "domain_aggregate": "NA"})
        gs = cal.get("gene_specific", "NA")
        da = cal.get("domain_aggregate", "NA")
        am_found = coverage.get(g, {}).get("am_found", False)
        coords = "yes" if g in coords_per_gene else "no"
        if gs != "NA":
            status = "gene_specific"; n_gs += 1
        elif da != "NA":
            status = "domain_only_recorded"; n_da_only += 1
        else:
            status = "uncovered_dropped"; n_uncov += 1
        lines.append(f"| {g} | {'yes' if am_found else 'no'} | {coords} | {gs} | {da} | {status} |")
    lines.append("")
    lines.append("## Coverage summary")
    lines.append(f"- gene_specific calibrated: **{n_gs}**")
    lines.append(f"- domain-aggregate only (RECORDED, not used for primary call): **{n_da_only}**")
    lines.append(f"- uncovered (retained with score, primary=FALSE): **{n_uncov}**")
    lines.append(f"- total: **{len(target_genes)}**")
    lines.append("")
    lines.append("## HARD GATE")
    lines.append("Review this report before submitting steps 01+.  In particular:")
    lines.append("- if calibration source is PLACEHOLDER, drop the real file at `calibration.local_file` and re-run `00_prepare_refs.sh`;")
    lines.append("- if exon BED source is `AM_minmax_fallback`, set `references.gencode_gtf_local` to a local GTF and re-run.")
    report_path = os.path.join(refs_dir, "REFERENCE_REPORT.md")
    with open(report_path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    LOG.info("wrote %s", report_path)

    # final hard-gate print (also visible in stdout)
    print()
    print(f"== reference coverage: gene_specific={n_gs}  domain_only={n_da_only}  uncovered={n_uncov} ==")
    print(f"== review {report_path} before submitting 01+ ==")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
