"""Paper-ready figure generators for the BioMe AlphaMissense rebuttal.

Each ``figN_*.py`` exposes ``make(out_dir, **inputs) -> None`` and writes both
``figN.png`` (300 dpi) and ``figN.pdf`` (vector). All styling is centralized
in :mod:`common`; do not set rcParams inside individual figure modules.
"""
