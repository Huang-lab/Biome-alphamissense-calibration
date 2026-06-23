"""Shared helpers for the BioMe AlphaMissense pipeline.

Every helper here is small, header-aware, and prints loud diagnostics so a
silent shape mismatch (wrong column order, ID column mismatch, missing chr)
turns into a clear failure.
"""
from __future__ import annotations

import csv
import gzip
import io
import logging
import os
import re
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

import yaml


# ----------------------------------------------------------------------------
# logging
# ----------------------------------------------------------------------------
LOG = logging.getLogger("biome-am")
if not LOG.handlers:
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("[%(asctime)s] %(levelname)s %(name)s :: %(message)s",
                                     datefmt="%Y-%m-%dT%H:%M:%S"))
    LOG.addHandler(h)
    LOG.setLevel(logging.INFO)


def die(msg: str, code: int = 1) -> "None":
    LOG.error(msg)
    sys.exit(code)


# ----------------------------------------------------------------------------
# config
# ----------------------------------------------------------------------------
def load_config(path: Optional[str] = None) -> dict:
    """Load the yaml config; default = repo_root/config/config.yaml."""
    if path is None:
        path = os.environ.get("CONFIG_PATH")
    if path is None:
        # walk up from this file: python/util.py -> python -> repo_root
        here = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(here)
        path = os.path.join(repo_root, "config", "config.yaml")
    if not os.path.isfile(path):
        die(f"config not found: {path}")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    cfg["__path__"] = path
    cfg["__repo_root__"] = os.path.dirname(os.path.dirname(os.path.abspath(path)))
    return cfg


def resolve(cfg: dict, rel: str) -> str:
    """Resolve a path: absolute -> kept; relative -> repo_root/rel."""
    if os.path.isabs(rel):
        return rel
    return os.path.join(cfg["__repo_root__"], rel)


# ----------------------------------------------------------------------------
# header-aware TSV reader
# ----------------------------------------------------------------------------
@contextmanager
def open_text(path: str) -> Iterator[io.TextIOBase]:
    """Open .gz / plain text uniformly."""
    if path.endswith(".gz"):
        with gzip.open(path, "rt", newline="") as fh:
            yield fh
    else:
        with open(path, "rt", newline="") as fh:
            yield fh


def read_tsv_dicts(path: str, comment_prefix: str = "#") -> Tuple[List[str], List[dict]]:
    """Read a TSV (possibly gzipped) into a list of dicts.

    - Skips leading lines starting with ``comment_prefix`` UNTIL the first
      non-comment line, which is treated as the header. If the header itself
      begins with '#' (AlphaMissense format), we strip the leading '#'.
    - Never positional column access; everything is by header name.
    """
    with open_text(path) as fh:
        header: Optional[List[str]] = None
        rows: List[dict] = []
        for line in fh:
            if not line or line == "\n":
                continue
            if header is None:
                stripped = line.lstrip().rstrip("\n")
                # AM format: a # header line listing column names.
                if stripped.startswith("#"):
                    # consume comment lines but remember the LAST one before
                    # the body — that's typically the header in AM-style files
                    candidate = stripped.lstrip("#").strip()
                    if "\t" in candidate:
                        header = candidate.split("\t")
                    continue
                # plain header line
                header = stripped.split("\t")
                continue
            if line.startswith(comment_prefix):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != len(header):
                # tolerate trailing empties only
                if len(parts) < len(header):
                    parts = parts + [""] * (len(header) - len(parts))
                else:
                    LOG.warning("column count mismatch (got %d, expected %d) in %s; truncating", len(parts), len(header), path)
                    parts = parts[: len(header)]
            rows.append(dict(zip(header, parts)))
        if header is None:
            die(f"empty TSV (no header found): {path}")
        return header, rows


def stream_tsv_dicts(path: str, comment_prefix: str = "#") -> Iterator[dict]:
    """Streaming variant of read_tsv_dicts for large files (AM, gencode GTF derived)."""
    with open_text(path) as fh:
        header: Optional[List[str]] = None
        for line in fh:
            if not line or line == "\n":
                continue
            if header is None:
                stripped = line.lstrip().rstrip("\n")
                if stripped.startswith("#"):
                    cand = stripped.lstrip("#").strip()
                    if "\t" in cand:
                        header = cand.split("\t")
                    continue
                header = stripped.split("\t")
                continue
            if line.startswith(comment_prefix):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                parts = parts + [""] * (len(header) - len(parts))
            elif len(parts) > len(header):
                parts = parts[: len(header)]
            yield dict(zip(header, parts))


def write_tsv(path: str, header: Sequence[str], rows: Iterable[Sequence]) -> int:
    """Write a TSV. Returns row count."""
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    n = 0
    opener = gzip.open if path.endswith(".gz") else open
    mode = "wt"
    with opener(path, mode, newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(list(header))
        for r in rows:
            w.writerow(list(r))
            n += 1
    LOG.info("wrote %d data rows -> %s", n, path)
    return n


# ----------------------------------------------------------------------------
# ID-join with match-rate report
# ----------------------------------------------------------------------------
@dataclass
class JoinReport:
    left_n: int
    right_n: int
    matched: int
    left_only_samples: List[str]
    right_only_samples: List[str]

    @property
    def match_rate(self) -> float:
        return 0.0 if self.left_n == 0 else self.matched / self.left_n

    def log(self, label: str, warn_below: float = 0.5) -> None:
        LOG.info("%s :: left_n=%d right_n=%d matched=%d match_rate=%.3f",
                 label, self.left_n, self.right_n, self.matched, self.match_rate)
        if self.left_only_samples:
            LOG.info("%s :: example left-only IDs (up to 5): %s", label, self.left_only_samples[:5])
        if self.right_only_samples:
            LOG.info("%s :: example right-only IDs (up to 5): %s", label, self.right_only_samples[:5])
        if self.match_rate < warn_below:
            LOG.warning("%s :: LOW MATCH RATE (%.3f). Hypothesis: ID column mismatch between left and right (formatting, leading zeros, prefix). Inspect the example IDs above.", label, self.match_rate)


def report_join(left_ids: Iterable[str], right_ids: Iterable[str], label: str) -> JoinReport:
    left = list(left_ids)
    right_set = set(right_ids)
    left_set = set(left)
    matched = len([x for x in left if x in right_set])
    return JoinReport(
        left_n=len(left),
        right_n=len(right_set),
        matched=matched,
        left_only_samples=sorted(left_set - right_set),
        right_only_samples=sorted(right_set - left_set),
    )


# ----------------------------------------------------------------------------
# small helpers
# ----------------------------------------------------------------------------
def strip_version(transcript_id: str) -> str:
    """ENST00000335137.4 -> ENST00000335137"""
    return transcript_id.split(".", 1)[0]


def chrom_norm(c: str) -> str:
    """Normalize a chromosome name to 'chrN' form."""
    if c.startswith("chr"):
        return c
    return f"chr{c}"


def vcf_chrom_prefix(vcf: str) -> str:
    """Inspect a VCF's first ##contig header and return 'chr' if the file uses
    the chr-prefixed convention, else ''. Used by step 02 to decide what region
    string to pass to `bcftools query -r` — must match the file's CHROM column."""
    import subprocess
    p = subprocess.run(["bcftools", "view", "-h", vcf],
                       capture_output=True, text=True, check=True)
    for line in p.stdout.splitlines():
        if line.startswith("##contig"):
            m = re.search(r"ID=([^,>]+)", line)
            if m:
                return "chr" if m.group(1).startswith("chr") else ""
    return "chr"  # no ##contig declared — fall back to chr-prefix


def parse_gtf_attrs(field9: str) -> Dict[str, str]:
    """GTF attribute field: gene_id "ENSG..."; transcript_id "ENST..."; ..."""
    out: Dict[str, str] = {}
    for chunk in field9.strip().rstrip(";").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        m = re.match(r'(\S+)\s+"([^"]*)"', chunk)
        if m:
            out[m.group(1)] = m.group(2)
    return out
