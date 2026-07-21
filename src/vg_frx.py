"""
vg GFA Post-Processor for FRX.

Provides VGFRXNormalizer: a pre-filter that converts a vg-generated GFA into
an FRX-compatible GFA before annotation. vg emits several quirks that require
normalization:

- Extra S-line tags vg injects: LN:i: (length), RC:i: (read count), FC:i:
- P lines (GFA 1.0 path records) - passed through unchanged
- W lines (GFA 1.1 walk records, newer vg versions) - passed through unchanged
- Numeric node IDs (vg uses integers) - FRX handles these natively, no change needed

After normalization, the output GFA is valid input for GFAAnnotator or
StreamingGFAAnnotator without modification.

Usage (module):
    from vg_frx import VGFRXNormalizer
    norm = VGFRXNormalizer()
    stats = norm.normalize("vg_output.gfa", "normalized.gfa")

Usage (CLI):
    frx vg-import --gfa vg_output.gfa --out index.db --build-kmer-index
"""

import os
import sys
import tempfile
from typing import Optional

# VG-internal tags to strip from S-lines (they collide with or clutter FRX tags)
_VG_STRIP_PREFIXES = ('LN:i:', 'RC:i:', 'FC:i:', 'SN:Z:', 'SO:i:', 'SR:i:')


class VGNormalizationError(ValueError):
    """Raised when the input GFA cannot be normalized."""


class VGFRXNormalizer:
    """
    Normalizes a vg-generated GFA file for FRX annotation.

    The normalizer performs a single-pass transformation:
    1. Strips vg-internal tags from S lines (optionally, unless keep_vg_tags=True).
    2. Passes H, L, P, W, and all other lines through unchanged.
    3. Skips empty lines.
    4. Reports per-line-type statistics.

    vg GFA compatibility:
        GFA 1.0 (S/L/P) - fully supported
        GFA 1.1 (S/L/W) - W lines passed through, FRX ignores them
    """

    # These S-line tag prefixes are added by vg and mean nothing to FRX
    VG_TAG_PREFIXES = _VG_STRIP_PREFIXES

    def __init__(self, keep_vg_tags: bool = False):
        """
        Args:
            keep_vg_tags: If True, do not strip vg-internal tags.
                          Useful for debugging or when downstream tools need them.
        """
        self.keep_vg_tags = keep_vg_tags

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def normalize(self, input_path: str, output_path: str) -> dict:
        """
        Read input_path, strip vg internals, write to output_path.

        Args:
            input_path  : Path to vg-generated GFA file.
            output_path : Path to write the normalized GFA.

        Returns:
            dict with keys: s_lines, l_lines, p_lines, w_lines,
                            h_lines, tags_stripped, total_lines.

        Raises:
            FileNotFoundError : input_path does not exist.
            VGNormalizationError : file is not recognizable as GFA.
        """
        if not os.path.isfile(input_path):
            raise FileNotFoundError(f"vg GFA not found: {input_path}")

        stats = {
            'h_lines': 0,
            's_lines': 0,
            'l_lines': 0,
            'p_lines': 0,
            'w_lines': 0,
            'other_lines': 0,
            'tags_stripped': 0,
            'total_lines': 0,
        }

        seen_gfa = False
        with open(input_path, 'r', encoding='utf-8') as fin, \
             open(output_path, 'w', encoding='utf-8') as fout:

            for raw_line in fin:
                stats['total_lines'] += 1
                line = raw_line.rstrip('\n')

                if not line:
                    continue

                record = line.split('\t')
                rtype = record[0]

                if rtype == 'H':
                    seen_gfa = True
                    stats['h_lines'] += 1
                    fout.write(line + '\n')

                elif rtype == 'S':
                    seen_gfa = True
                    stats['s_lines'] += 1
                    out_line, stripped = self._normalize_s_line(record)
                    stats['tags_stripped'] += stripped
                    fout.write(out_line + '\n')

                elif rtype == 'L':
                    stats['l_lines'] += 1
                    fout.write(line + '\n')

                elif rtype == 'P':
                    stats['p_lines'] += 1
                    fout.write(line + '\n')

                elif rtype == 'W':
                    stats['w_lines'] += 1
                    fout.write(line + '\n')

                else:
                    stats['other_lines'] += 1
                    fout.write(line + '\n')

        if not seen_gfa:
            raise VGNormalizationError(
                f"'{input_path}' does not appear to be a GFA file "
                f"(no H or S records found)."
            )

        return stats

    def normalize_to_tempfile(self, input_path: str) -> tuple:
        """
        Normalize to a system temporary file.

        Returns:
            (tmp_path: str, stats: dict)
            Caller is responsible for deleting tmp_path when done.
        """
        fd, tmp_path = tempfile.mkstemp(suffix='_frx_normalized.gfa')
        os.close(fd)
        stats = self.normalize(input_path, tmp_path)
        return tmp_path, stats

    def detect_gfa_version(self, path: str) -> str:
        """
        Heuristically detect GFA version from file content.

        Returns:
            "1.1" if W-lines are found, "1.0" otherwise.
        """
        with open(path, 'r', encoding='utf-8') as f:
            for line in f:
                if line.startswith('W\t'):
                    return "1.1"
                if line.startswith('H\t'):
                    # Check VN:Z: header tag
                    for part in line.strip().split('\t')[1:]:
                        if part.startswith('VN:Z:'):
                            return part[5:].strip()
        return "1.0"

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _normalize_s_line(self, record: list) -> tuple:
        """
        Strip vg-internal tags from an S-line record.

        Args:
            record: List of tab-split fields from a GFA S-line.

        Returns:
            (normalized_line: str, tags_stripped: int)
        """
        if len(record) < 3:
            return '\t'.join(record), 0

        # Fields 0,1,2 are mandatory: record type, node_id, sequence
        mandatory = record[:3]
        optional = record[3:]

        if self.keep_vg_tags:
            return '\t'.join(record), 0

        kept = []
        stripped = 0
        for tag in optional:
            if any(tag.startswith(pfx) for pfx in self.VG_TAG_PREFIXES):
                stripped += 1
            else:
                kept.append(tag)

        return '\t'.join(mandatory + kept), stripped


# ======================================================================
# Standalone pipeline runner (used by frx vg-import)
# ======================================================================

def run_vg_import(
    gfa_path: str,
    out_db: str,
    seed: bytes = b"panindex_default_seed",
    derivation_root: str = "PangenomeRoot",
    streaming: bool = False,
    build_kmer_index: bool = False,
    keep_vg_tags: bool = False,
    verbose: bool = True,
) -> dict:
    """
    Full vg -> FRX pipeline:
      1. Normalize vg GFA
      2. Annotate (streaming or standard)
      3. Save .frx.db
      4. Optionally build k-mer index

    Args:
        gfa_path         : Path to vg-generated GFA file.
        out_db           : Destination .frx.db path.
        seed             : Pangenome seed bytes.
        derivation_root  : Root label for the ratchet hierarchy.
        streaming        : Use StreamingGFAAnnotator for low RAM usage.
        build_kmer_index : Build and save a k-mer index after annotation.
        keep_vg_tags     : Keep vg-internal tags in normalized GFA.
        verbose          : Print progress messages.

    Returns:
        dict with keys: normalized_stats, nodes_written, kmer_entries (or 0).
    """
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

    from streaming_annotator import make_annotator
    from persistence import db_stats

    norm = VGFRXNormalizer(keep_vg_tags=keep_vg_tags)

    if verbose:
        print(f"[frx vg-import] Normalizing vg GFA: {gfa_path}")

    tmp_path, norm_stats = norm.normalize_to_tempfile(gfa_path)

    if verbose:
        gfa_ver = norm.detect_gfa_version(gfa_path)
        print(f"  GFA version    : {gfa_ver}")
        print(f"  S lines        : {norm_stats['s_lines']}")
        print(f"  L lines        : {norm_stats['l_lines']}")
        print(f"  P lines        : {norm_stats['p_lines']}")
        print(f"  W lines        : {norm_stats['w_lines']}")
        print(f"  Tags stripped  : {norm_stats['tags_stripped']}")

    try:
        annotated_path = tmp_path + '_annotated.gfa'
        ann = make_annotator(tmp_path, seed=seed, force_streaming=streaming)

        if verbose:
            print(f"[frx vg-import] Annotating...")

        ann.annotate(tmp_path, annotated_path, derivation_root=derivation_root)
        nodes_written = ann.store.save(out_db)

        kmer_entries = 0
        if build_kmer_index:
            from kmer_index import KmerIndex
            if verbose:
                print(f"[frx vg-import] Building k-mer index (k=12)...")
            ki = KmerIndex.build(ann.store, k=12)
            kmer_entries = ki.save(out_db)
            if verbose:
                st = ki.stats()
                print(f"  Distinct k-mers: {st['distinct_kmers']:,}")

        st = db_stats(out_db)
        if verbose:
            print(f"[frx vg-import] Done.")
            print(f"  Nodes written  : {nodes_written}")
            print(f"  Unique tags    : {st['unique_tags']}")
            print(f"  DB size        : {st['db_size_bytes']:,} bytes")
            print(f"  Index saved to : {out_db}")

    finally:
        # Always clean up temp files
        for p in [tmp_path, tmp_path + '_annotated.gfa']:
            try:
                os.unlink(p)
            except FileNotFoundError:
                pass

    return {
        'normalized_stats': norm_stats,
        'nodes_written': nodes_written,
        'kmer_entries': kmer_entries,
    }
