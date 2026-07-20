"""
FRX FASTA+GFA Merger and Subsequence Query Tool.

Two problems solved here:

Problem 1 - Fragmented input formats:
  Many pangenome tools (Minigraph, PGGB, Bandage) produce GFA files where
  S-lines carry '*' as the sequence placeholder, with actual sequences stored
  in a separate FASTA. This module merges the two back into a single GFA where
  every S-line carries its real nucleotide sequence, which the full PanIndex
  pipeline can then address and index.

Problem 2 - No FASTA-style coordinate queries in GFA tools:
  Standard GFA tools index graph topology but not sequence sub-regions. You
  cannot ask "give me segment seg1 from position 100 to 200" the way
  samtools faidx does for FASTA. This module adds that:

    SubsequenceQuery.query("seg1:100-200")
    SubsequenceQuery.search_pattern("ATGCGT")

  Each extracted region also receives a stable HKDF-derived address, making
  the sub-region a first-class PanIndex object that can be stored and looked
  up without re-running the query.

CLI:
  python src/fasta_merge.py merge      genome.fasta graph.gfa merged.gfa
  python src/fasta_merge.py fasta2gfa  genome.fasta output.gfa
  python src/fasta_merge.py region     --gfa merged.gfa --region seg1:100-200
  python src/fasta_merge.py pattern    --gfa merged.gfa --pattern ATGCGT
"""

import sys
import os
import argparse
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore
from annotator import GFAAnnotator


# ======================================================================
# FASTA Parsing
# ======================================================================

class FastaParser:
    """
    Parse a multi-record FASTA file into a name -> sequence mapping.

    Handles:
    - Multi-line FASTA records (sequence split across multiple lines)
    - Description lines (takes only the first whitespace-separated token as name)
    - Sequences are uppercased on load
    """

    @staticmethod
    def parse(fasta_path: str) -> Dict[str, str]:
        """
        Returns dict of {record_name: sequence}.
        record_name is the first word on the '>' header line (no '>').
        """
        sequences: Dict[str, str] = {}
        current_name: Optional[str] = None
        current_seq: List[str] = []

        with open(fasta_path, 'r', encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if not line:
                    continue
                if line.startswith('>'):
                    if current_name is not None:
                        sequences[current_name] = ''.join(current_seq)
                    current_name = line[1:].split()[0]
                    current_seq = []
                elif current_name is not None:
                    current_seq.append(line.upper())

        if current_name is not None:
            sequences[current_name] = ''.join(current_seq)

        return sequences

    @staticmethod
    def stats(sequences: Dict[str, str]) -> Dict[str, int]:
        """Basic statistics for a parsed FASTA dict."""
        lengths = [len(s) for s in sequences.values()]
        return {
            'records': len(sequences),
            'total_bp': sum(lengths),
            'min_len': min(lengths) if lengths else 0,
            'max_len': max(lengths) if lengths else 0,
        }


# ======================================================================
# FASTA + GFA Merger
# ======================================================================

class FastaMerger:
    """
    Merges a FASTA sequence file with a GFA topology file.

    Segment name matching strategy (in order):
    a. Exact match on segment name from the GFA S-line
    b. Strip leading 's'/'S' prefix (e.g. 's1' -> '1') and retry

    Three S-line cases:
    1. Sequence is '*'            -> fill from FASTA by segment name
    2. Sequence already present   -> keep unchanged, count as kept
    3. Segment not found in FASTA -> warn, keep '*', add to missing list
    """

    def __init__(self):
        self.fasta_sequences: Dict[str, str] = {}
        self.missing: List[str] = []
        self.filled: int = 0
        self.kept: int = 0

    def _resolve_name(self, seg_id: str) -> Optional[str]:
        """Return the FASTA key that matches seg_id, or None."""
        if seg_id in self.fasta_sequences:
            return seg_id
        alt = seg_id.lstrip('sS')
        if alt and alt in self.fasta_sequences:
            return alt
        return None

    def merge(self, fasta_path: str, gfa_path: str, output_path: str) -> Dict[str, int]:
        """
        Merge FASTA sequences into a GFA file and write the result.

        Args:
            fasta_path  : Input FASTA file with sequences.
            gfa_path    : Input GFA file (may have '*' in S-lines).
            output_path : Output merged GFA file path.

        Returns:
            Stats dict: fasta_records, filled, kept, missing.
        """
        self.fasta_sequences = FastaParser.parse(fasta_path)
        st = FastaParser.stats(self.fasta_sequences)
        print(f"[FastaMerger] Loaded {st['records']} FASTA records "
              f"({st['total_bp']:,} bp total) from '{fasta_path}'")

        self.missing = []
        self.filled = 0
        self.kept = 0

        with open(gfa_path, 'r', encoding='utf-8') as fin, \
             open(output_path, 'w', encoding='utf-8') as fout:

            for line in fin:
                parts = line.rstrip('\n').split('\t')

                if parts[0] == 'S' and len(parts) >= 3:
                    seg_id = parts[1]
                    seq = parts[2]

                    if seq == '*':
                        resolved = self._resolve_name(seg_id)
                        if resolved:
                            parts[2] = self.fasta_sequences[resolved]
                            self.filled += 1
                        else:
                            self.missing.append(seg_id)
                    else:
                        self.kept += 1

                    fout.write('\t'.join(parts) + '\n')
                else:
                    fout.write(line)

        stats = {
            'fasta_records': len(self.fasta_sequences),
            'filled': self.filled,
            'kept': self.kept,
            'missing': len(self.missing),
        }

        print(f"[FastaMerger] Merge complete -> '{output_path}'")
        print(f"  Filled from FASTA : {self.filled}")
        print(f"  Already had seq   : {self.kept}")
        print(f"  Still missing     : {len(self.missing)}")
        if self.missing:
            shown = self.missing[:10]
            suffix = '...' if len(self.missing) > 10 else ''
            print(f"  Missing segments  : {shown}{suffix}")

        return stats

    def fasta_as_gfa(self, fasta_path: str, output_path: str) -> str:
        """
        Convert a standalone FASTA file directly to GFA 1.0.

        Each FASTA record becomes one S-line segment with no L-links.
        The resulting GFA is compatible with the full PanIndex pipeline,
        enabling address derivation and all three query modes on FASTA data.

        Args:
            fasta_path  : Input FASTA file.
            output_path : Output GFA file path.

        Returns:
            output_path
        """
        sequences = FastaParser.parse(fasta_path)
        st = FastaParser.stats(sequences)
        source_name = os.path.basename(fasta_path)

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(f"H\tVN:Z:1.0\tPN:Z:{source_name}\n")
            for name, seq in sequences.items():
                f.write(f"S\t{name}\t{seq}\tLN:i:{len(seq)}\n")

        print(f"[FastaMerger] FASTA -> GFA: {st['records']} records "
              f"({st['total_bp']:,} bp) -> '{output_path}'")
        return output_path


# ======================================================================
# Subsequence / Region Query
# ======================================================================

class SubsequenceQuery:
    """
    FASTA-style subsequence extraction and pattern search on a PanIndexStore.

    Adds two capabilities that standard GFA tools lack:

    1. Region query  : "segment_id:start-end" (0-based, end exclusive)
       Returns the subsequence plus a stable HKDF-derived region address.
       Each unique region gets a unique 32-byte address derived from its
       parent segment address and its coordinate string.

    2. Pattern search: exact nucleotide substring -> all occurrences
       across all segments. Uses KmerIndex for O(|P|/k) average search
       when available, falls back to O(N*L) linear scan otherwise.
    """

    def __init__(
        self,
        engine: PanIndexEngine,
        store: PanIndexStore,
        kmer_index=None,
    ):
        self.engine = engine
        self.store = store
        self.kmer_index = kmer_index  # KmerIndex instance or None

    # ------------------------------------------------------------------
    # Region extraction
    # ------------------------------------------------------------------

    def extract_region(self, segment_id: str, start: int, end: int) -> Optional[Dict]:
        """
        Extract a subsequence from a stored segment by 0-based coordinate.

        Args:
            segment_id : GFA segment name (S-line field 1).
            start      : Start position, 0-based inclusive.
            end        : End position, 0-based exclusive.

        Returns:
            Dict with keys: segment_id, start, end, length, subsequence,
                            region_address (bytes), parent_address (hex str),
                            region_context (str).
            None if segment not found or coordinates produce empty range.
        """
        node = self.store.get_node(segment_id)
        if node is None:
            return None

        seq = node['metadata'].get('seq', '')
        if not seq:
            return None

        start = max(0, start)
        end = min(len(seq), end)

        if start >= end:
            return None

        subseq = seq[start:end]

        # Derive a stable HKDF address for this exact sub-region.
        # region_context encodes segment + coordinates so two regions of the
        # same segment always produce distinct addresses.
        region_context = f"{segment_id}:{start}-{end}"
        parent_addr = bytes.fromhex(node['address'])
        region_address = self.engine.derive_ratchet_address(parent_addr, region_context)

        return {
            'segment_id': segment_id,
            'start': start,
            'end': end,
            'length': end - start,
            'subsequence': subseq,
            'region_address': region_address,
            'parent_address': node['address'],
            'region_context': region_context,
        }

    def query(self, region: str) -> Optional[Dict]:
        """
        Extract a region by FASTA-style string.

        Args:
            region: "segment_id:start-end" (e.g. "seg1:100-200")

        Returns:
            Same dict as extract_region, or None on miss.
        """
        seg_id, start, end = self._parse_region_string(region)
        return self.extract_region(seg_id, start, end)

    # ------------------------------------------------------------------
    # Pattern search
    # ------------------------------------------------------------------

    def search_pattern(self, pattern: str) -> List[Dict]:
        """
        Find all exact occurrences of a nucleotide pattern across every
        indexed segment in the store.

        Routing:
        - KmerIndex present (k=12 default): O(|pattern|/k) average via seeded search.
        - No KmerIndex: O(N * L) linear scan fallback.

        Args:
            pattern: Nucleotide sequence to search for (case-insensitive).

        Returns:
            List of hit dicts sorted by (segment_id, start).
        """
        if self.kmer_index is not None:
            return self.kmer_index.search(pattern.upper(), self.store, self.engine)

        # Linear scan fallback
        pattern = pattern.upper()
        hits: List[Dict] = []

        for node_id in self.store.all_nodes():
            node = self.store.get_node(node_id)
            seq = node['metadata'].get('seq', '')
            if not seq:
                continue

            pos = 0
            while True:
                idx = seq.find(pattern, pos)
                if idx == -1:
                    break

                end = idx + len(pattern)
                region_context = f"{node_id}:{idx}-{end}"
                parent_addr = bytes.fromhex(node['address'])
                region_address = self.engine.derive_ratchet_address(
                    parent_addr, region_context
                )

                hits.append({
                    'segment_id': node_id,
                    'start': idx,
                    'end': end,
                    'length': len(pattern),
                    'subsequence': pattern,
                    'region_address': region_address,
                    'parent_address': node['address'],
                    'region_context': region_context,
                    'tags': node['tags'],
                })

                pos = idx + 1

        hits.sort(key=lambda h: (h['segment_id'], h['start']))
        return hits

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_region_string(region: str) -> Tuple[str, int, int]:
        """Parse "segment_id:start-end" into (segment_id, start, end)."""
        if ':' not in region:
            raise ValueError(
                f"Region must be 'segment_id:start-end', got '{region}'"
            )
        seg_id, coords = region.rsplit(':', 1)
        if '-' not in coords:
            raise ValueError(
                f"Coordinates must be 'start-end', got '{coords}'"
            )
        start_s, end_s = coords.split('-', 1)
        try:
            return seg_id, int(start_s), int(end_s)
        except ValueError:
            raise ValueError(
                f"Coordinates must be integers, got '{coords}'"
            )

    def print_region_result(self, result: Dict):
        """Pretty-print a single region extraction result."""
        print(f"\n[Region] {result['region_context']}")
        print(f"  Length    : {result['length']} bp")
        print(f"  Sequence  : {result['subsequence']}")
        print(f"  Address   : {result['region_address'].hex()[:24]}...")
        print(f"  Parent    : {result['parent_address'][:24]}...")

    def print_pattern_results(self, hits: List[Dict], pattern: str):
        """Pretty-print pattern search results."""
        print(f"\n[Pattern Search] '{pattern}' -> {len(hits)} hit(s)")
        if not hits:
            print("  No matches found.")
            return
        for hit in hits:
            tags = ', '.join(hit['tags']) if hit['tags'] else 'none'
            print(f"  {hit['region_context']}")
            print(f"    Seq  : {hit['subsequence']}")
            print(f"    Addr : {hit['region_address'].hex()[:24]}...")
            print(f"    Tags : {tags}")


# ======================================================================
# Shared GFA loader
# ======================================================================

def _load_gfa_into_store(
    gfa_path: str, seed: bytes
) -> Tuple[PanIndexEngine, PanIndexStore, str]:
    """
    Run the annotator pipeline on a GFA file to populate engine + store.
    Writes an annotated copy alongside the input.
    Returns (engine, store, annotated_output_path).
    """
    base, ext = os.path.splitext(gfa_path)
    annotated_path = base + '_frx_annotated' + ext
    print(f"[FRX] Loading GFA and building index: {gfa_path}")
    ann = GFAAnnotator(seed=seed)
    ann.annotate(gfa_path, annotated_path)
    print(f"[FRX] Index built. {len(ann.store)} segments indexed.")
    return ann.engine, ann.store, annotated_path


# ======================================================================
# CLI
# ======================================================================

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='fasta_merge',
        description='FRX FASTA+GFA Merger and Subsequence Query Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
commands:
  merge      Merge FASTA sequences into a GFA file (fills * placeholders)
  fasta2gfa  Convert a standalone FASTA directly to GFA
  region     Extract a subsequence by FASTA-style coordinate string
  pattern    Search for a nucleotide pattern across all graph segments

examples:
  python src/fasta_merge.py merge genome.fasta graph.gfa merged.gfa
  python src/fasta_merge.py merge genome.fasta graph.gfa merged.gfa --annotate
  python src/fasta_merge.py fasta2gfa genome.fasta output.gfa --annotate
  python src/fasta_merge.py region  --gfa merged.gfa --region seg1:100-200
  python src/fasta_merge.py pattern --gfa merged.gfa --pattern ATGCGT
        """
    )

    sub = p.add_subparsers(dest='command', metavar='command')

    # merge
    pm = sub.add_parser('merge', help='Merge FASTA into GFA')
    pm.add_argument('fasta',  help='Input FASTA file')
    pm.add_argument('gfa',    help='Input GFA file (may have * sequences)')
    pm.add_argument('output', help='Output merged GFA file')
    pm.add_argument('--annotate', action='store_true',
                    help='Run PanIndex annotation after merge')
    pm.add_argument('--seed', default='panindex_default_seed',
                    help='Pangenome seed string for annotation')

    # fasta2gfa
    pf = sub.add_parser('fasta2gfa', help='Convert FASTA to GFA')
    pf.add_argument('fasta',  help='Input FASTA file')
    pf.add_argument('output', help='Output GFA file')
    pf.add_argument('--annotate', action='store_true',
                    help='Run PanIndex annotation after conversion')
    pf.add_argument('--seed', default='panindex_default_seed')

    # region
    pr = sub.add_parser('region', help='Extract subsequence by coordinate')
    pr.add_argument('--gfa',    required=True, help='GFA file with sequences embedded')
    pr.add_argument('--region', required=True,
                    help='Region string: segment_id:start-end (0-based, end exclusive)')
    pr.add_argument('--seed', default='panindex_default_seed')

    # pattern
    pp = sub.add_parser('pattern', help='Search pattern across all segments')
    pp.add_argument('--gfa',     required=True, help='GFA file with sequences embedded')
    pp.add_argument('--pattern', required=True, help='Nucleotide sequence pattern')
    pp.add_argument('--seed', default='panindex_default_seed')

    return p


def main():
    parser = _build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == 'merge':
        merger = FastaMerger()
        merger.merge(args.fasta, args.gfa, args.output)
        if args.annotate:
            base, ext = os.path.splitext(args.output)
            annotated = base + '_annotated' + ext
            ann = GFAAnnotator(seed=args.seed.encode())
            ann.annotate(args.output, annotated)
            ann.print_index_summary()
            print(f"\nAnnotated GFA written to: {annotated}")

    elif args.command == 'fasta2gfa':
        merger = FastaMerger()
        merger.fasta_as_gfa(args.fasta, args.output)
        if args.annotate:
            base, ext = os.path.splitext(args.output)
            annotated = base + '_annotated' + ext
            ann = GFAAnnotator(seed=args.seed.encode())
            ann.annotate(args.output, annotated)
            ann.print_index_summary()
            print(f"\nAnnotated GFA written to: {annotated}")

    elif args.command == 'region':
        engine, store, _ = _load_gfa_into_store(args.gfa, args.seed.encode())
        sq = SubsequenceQuery(engine, store)
        result = sq.query(args.region)
        if result:
            sq.print_region_result(result)
        else:
            print(f"No result for region '{args.region}'. "
                  f"Check segment name and coordinate bounds.")
            sys.exit(1)

    elif args.command == 'pattern':
        engine, store, _ = _load_gfa_into_store(args.gfa, args.seed.encode())
        sq = SubsequenceQuery(engine, store)
        hits = sq.search_pattern(args.pattern)
        sq.print_pattern_results(hits, args.pattern)


if __name__ == '__main__':
    main()
