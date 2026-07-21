"""
frx - Unified CLI for FractalIndex (FRX).

All FRX operations available through a single command after pip install:

  frx annotate   <input.gfa> <output.gfa> [--seed STR]
  frx save-index --gfa <annotated.gfa> --out <index.db> [--seed STR]
  frx merge      <genome.fasta> <graph.gfa> <merged.gfa> [--annotate] [--seed STR]
  frx fasta2gfa  <genome.fasta> <output.gfa> [--annotate] [--seed STR]
  frx query      --index <index.db> (--path STR | --tag STR | --sequence STR)
  frx region     --index <index.db> --region <seg:start-end>
  frx pattern    --index <index.db> --pattern <ATGCGT>
  frx stats      --index <index.db>
  frx hgt-sim

Existing 'python src/X.py' scripts continue to work unchanged.
This CLI is an additional, unified entry point.
"""

import sys
import os
import argparse

# Add src/ to path so this module can be run both as a script and as
# the installed console entry point frx = "src.cli:main"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))


# ======================================================================
# Subcommand implementations
# ======================================================================

def cmd_annotate(args):
    """Parse a GFA, compute PanIndex addresses, write annotated GFA."""
    from annotator import GFAAnnotator

    seed = args.seed.encode()
    ann = GFAAnnotator(seed=seed)
    print(f"[frx annotate] {args.input} -> {args.output}")
    ann.annotate(args.input, args.output, derivation_root=args.root)
    ann.print_index_summary()
    print(f"\nDone. {len(ann.store)} nodes annotated.")


def cmd_save_index(args):
    """
    Annotate a GFA and save the resulting PanIndexStore to a SQLite db.

    Supports two annotation paths:
    - Standard (default): fast in-memory GFAAnnotator
    - Streaming (--streaming): O(1) RAM StreamingGFAAnnotator for large files

    With --build-kmer-index, also builds and saves a k-mer inverted index
    (k=12) enabling O(|pattern|/k) average-case pattern search.
    """
    from streaming_annotator import make_annotator
    from kmer_index import KmerIndex
    from persistence import db_stats

    seed = args.seed.encode()
    ann = make_annotator(
        args.gfa,
        seed=seed,
        force_streaming=args.streaming,
    )

    print(f"[frx save-index] Annotating: {args.gfa}")
    ann.annotate(args.gfa, args.gfa + '.frx_annotated.gfa',
                 derivation_root=args.root)

    print(f"[frx save-index] Saving index to: {args.out}")
    n = ann.store.save(args.out)

    if args.build_kmer_index:
        print(f"[frx save-index] Building k-mer index (k=12)...")
        ki = KmerIndex.build(ann.store, k=12)
        ki_rows = ki.save(args.out)
        st_ki = ki.stats()
        print(f"  K-mer entries  : {ki_rows:,}")
        print(f"  Distinct k-mers: {st_ki['distinct_kmers']:,}")

    st = db_stats(args.out)
    print(f"  Nodes written : {n}")
    print(f"  Unique tags   : {st['unique_tags']}")
    print(f"  DB size       : {st['db_size_bytes']:,} bytes")
    print(f"Done. Index saved to '{args.out}'")


def cmd_merge(args):
    """Merge a FASTA file into a GFA file (fills * sequence placeholders)."""
    from fasta_merge import FastaMerger, GFAAnnotator

    merger = FastaMerger()
    merger.merge(args.fasta, args.gfa, args.output)

    if args.annotate:
        base, ext = os.path.splitext(args.output)
        annotated = base + '_annotated' + ext
        ann = GFAAnnotator(seed=args.seed.encode())
        ann.annotate(args.output, annotated)
        ann.print_index_summary()
        print(f"\nAnnotated GFA: {annotated}")

        if args.save_index:
            db_path = base + '.frx.db'
            ann.store.save(db_path)
            print(f"Index saved  : {db_path}")


def cmd_fasta2gfa(args):
    """Convert a standalone FASTA file to GFA 1.0 format."""
    from fasta_merge import FastaMerger, GFAAnnotator

    merger = FastaMerger()
    merger.fasta_as_gfa(args.fasta, args.output)

    if args.annotate:
        base, ext = os.path.splitext(args.output)
        annotated = base + '_annotated' + ext
        ann = GFAAnnotator(seed=args.seed.encode())
        ann.annotate(args.output, annotated)
        ann.print_index_summary()
        print(f"\nAnnotated GFA: {annotated}")

        if args.save_index:
            db_path = base + '.frx.db'
            ann.store.save(db_path)
            print(f"Index saved  : {db_path}")


def cmd_query(args):
    """Query a saved FRX index by ratchet path, tag, or sequence similarity."""
    from index import PanIndexStore
    from engine import PanIndexEngine
    from query import PanIndexQuery

    if not any([args.path, args.tag, args.sequence]):
        print("Error: provide --path, --tag, or --sequence", file=sys.stderr)
        sys.exit(1)

    print(f"[frx query] Loading index: {args.index}")
    store = PanIndexStore.load(args.index)
    engine = PanIndexEngine(pangenome_seed=args.seed.encode())
    q = PanIndexQuery(engine, store)

    organism = args.organism if args.organism else None

    if args.path:
        q.query_by_path(args.path).print()
    if args.tag:
        nodes = store.lookup_by_tag(args.tag, organism=organism)
        print(f"\n[Tag Query] '{args.tag}' -> {len(nodes)} hit(s)")
        if organism:
            print(f"  Organism scope: {organism}")
        for nid in nodes:
            print(f"  {nid}")
    if args.sequence:
        q.query_by_similarity(args.sequence, threshold=args.threshold).print()


def cmd_region(args):
    """Extract a subsequence from a saved index by FASTA-style coordinates."""
    from index import PanIndexStore
    from engine import PanIndexEngine
    from fasta_merge import SubsequenceQuery

    print(f"[frx region] Loading index: {args.index}")
    store = PanIndexStore.load(args.index)
    engine = PanIndexEngine(pangenome_seed=args.seed.encode())
    sq = SubsequenceQuery(engine, store)

    result = sq.query(args.region)
    if result:
        sq.print_region_result(result)
    else:
        print(f"No result for region '{args.region}'. "
              f"Check segment name and coordinate bounds.")
        sys.exit(1)


def cmd_pattern(args):
    """Search for a nucleotide pattern across all segments in a saved index."""
    from index import PanIndexStore
    from engine import PanIndexEngine
    from fasta_merge import SubsequenceQuery
    from kmer_index import KmerIndex

    print(f"[frx pattern] Loading index: {args.index}")
    store = PanIndexStore.load(args.index)
    engine = PanIndexEngine(pangenome_seed=args.seed.encode())

    kmer_idx = None
    if KmerIndex.is_present(args.index):
        print(f"[frx pattern] K-mer index found. Using seeded search.")
        kmer_idx = KmerIndex.load(args.index)
    else:
        print(f"[frx pattern] No k-mer index. Using linear scan.")

    sq = SubsequenceQuery(engine, store, kmer_index=kmer_idx)
    hits = sq.search_pattern(args.pattern)
    sq.print_pattern_results(hits, args.pattern)


def cmd_stats(args):
    """Print statistics about a saved FRX index without fully loading it."""
    from persistence import db_stats
    import os

    st = db_stats(args.index)
    size_kb = st['db_size_bytes'] / 1024
    print(f"\n[frx stats] {args.index}")
    print(f"  Nodes       : {st['nodes']:,}")
    print(f"  Tag entries : {st['tags']:,}")
    print(f"  Unique tags : {st['unique_tags']:,}")
    print(f"  DB size     : {size_kb:.1f} KB")


def cmd_hgt_sim(args):
    """Run the built-in HGT simulation demo (K. pneumoniae -> E. coli blaTEM transfer)."""
    from hgt_handler import HGTSimulation

    sim = HGTSimulation()
    sim.run()


def cmd_variant(args):
    """Compare two indexed nodes with hash-based variant detection (Mode 4)."""
    from index import PanIndexStore
    from engine import PanIndexEngine
    from query import PanIndexQuery

    print(f"[frx variant] Loading index: {args.index}")
    store = PanIndexStore.load(args.index)
    engine = PanIndexEngine(pangenome_seed=args.seed.encode())
    q = PanIndexQuery(engine, store)

    result = q.query_variant_diff(
        args.node1,
        args.node2,
        similarity_threshold=args.threshold,
    )
    result.print()


def cmd_vg_import(args):
    """Normalize a vg-generated GFA and build an FRX index from it."""
    from vg_frx import run_vg_import

    run_vg_import(
        gfa_path=args.gfa,
        out_db=args.out,
        seed=args.seed.encode(),
        derivation_root=args.root,
        streaming=args.streaming,
        build_kmer_index=args.build_kmer_index,
        keep_vg_tags=args.keep_vg_tags,
        verbose=True,
    )


# ======================================================================
# Parser construction
# ======================================================================

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='frx',
        description='FractalIndex (FRX) - Content-Addressed Pangenome Indexing',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
examples:
  frx annotate   input.gfa  output.gfa
  frx save-index --gfa output.gfa --out pangenome.frx.db
  frx query      --index pangenome.frx.db --tag AMR:blaTEM
  frx query      --index pangenome.frx.db --path PangenomeRoot/1
  frx region     --index pangenome.frx.db --region 1:0-4
  frx pattern    --index pangenome.frx.db --pattern CAAATAAG
  frx stats      --index pangenome.frx.db
  frx merge      genome.fasta graph.gfa merged.gfa --annotate --save-index
  frx fasta2gfa  genome.fasta output.gfa --annotate
  frx vg-import  --gfa vg_output.gfa --out index.db --build-kmer-index
  frx hgt-sim
        """
    )

    sub = p.add_subparsers(dest='command', metavar='command')
    sub.required = True

    # -- annotate --
    pa = sub.add_parser('annotate',
                        help='Annotate a GFA with PanIndex ratchet addresses')
    pa.add_argument('input',  help='Input GFA file')
    pa.add_argument('output', help='Output annotated GFA file')
    pa.add_argument('--seed', default='panindex_default_seed',
                    help='Pangenome seed string (default: panindex_default_seed)')
    pa.add_argument('--root', default='PangenomeRoot',
                    help='Derivation root label (default: PangenomeRoot)')
    pa.set_defaults(func=cmd_annotate)

    # -- save-index --
    ps = sub.add_parser('save-index',
                        help='Annotate a GFA and save the index to a SQLite db')
    ps.add_argument('--gfa',  required=True, help='Input GFA file')
    ps.add_argument('--out',  required=True, help='Output .db file path')
    ps.add_argument('--seed', default='panindex_default_seed')
    ps.add_argument('--root', default='PangenomeRoot')
    ps.add_argument('--streaming', action='store_true',
                    help='Use streaming annotator (O(1) RAM, for large GFA files)')
    ps.add_argument('--build-kmer-index', action='store_true', dest='build_kmer_index',
                    help='Build and save a k-mer index for fast pattern search')
    ps.set_defaults(func=cmd_save_index)

    # -- merge --
    pm = sub.add_parser('merge',
                        help='Merge FASTA sequences into a GFA file')
    pm.add_argument('fasta',  help='Input FASTA file')
    pm.add_argument('gfa',    help='Input GFA file (may have * sequences)')
    pm.add_argument('output', help='Output merged GFA file')
    pm.add_argument('--annotate',   action='store_true',
                    help='Run PanIndex annotation after merge')
    pm.add_argument('--save-index', action='store_true',
                    help='Save the index to a .frx.db file (requires --annotate)')
    pm.add_argument('--seed', default='panindex_default_seed')
    pm.set_defaults(func=cmd_merge)

    # -- fasta2gfa --
    pf = sub.add_parser('fasta2gfa',
                        help='Convert a standalone FASTA to GFA')
    pf.add_argument('fasta',  help='Input FASTA file')
    pf.add_argument('output', help='Output GFA file')
    pf.add_argument('--annotate',   action='store_true')
    pf.add_argument('--save-index', action='store_true',
                    help='Save the index after annotation')
    pf.add_argument('--seed', default='panindex_default_seed')
    pf.set_defaults(func=cmd_fasta2gfa)

    # -- query --
    pq = sub.add_parser('query',
                        help='Query a saved FRX index')
    pq.add_argument('--index',    required=True, help='Path to .frx.db file')
    pq.add_argument('--path',     help='Ratchet path (e.g. PangenomeRoot/1)')
    pq.add_argument('--tag',      help='Anubandha tag (e.g. AMR:blaTEM)')
    pq.add_argument('--sequence', help='Nucleotide sequence for similarity search')
    pq.add_argument('--threshold', type=float, default=0.5,
                    help='MinHash similarity threshold (default: 0.5)')
    pq.add_argument('--organism', default='',
                    help='Scope tag queries to a specific organism/strain')
    pq.add_argument('--seed', default='panindex_default_seed')
    pq.set_defaults(func=cmd_query)

    # -- region --
    pr = sub.add_parser('region',
                        help='Extract a subsequence by coordinate')
    pr.add_argument('--index',  required=True, help='Path to .frx.db file')
    pr.add_argument('--region', required=True,
                    help='Region string: segment_id:start-end (0-based, end exclusive)')
    pr.add_argument('--seed', default='panindex_default_seed')
    pr.set_defaults(func=cmd_region)

    # -- pattern --
    pp = sub.add_parser('pattern',
                        help='Search a nucleotide pattern across all segments')
    pp.add_argument('--index',   required=True, help='Path to .frx.db file')
    pp.add_argument('--pattern', required=True, help='Nucleotide sequence pattern')
    pp.add_argument('--seed', default='panindex_default_seed')
    pp.set_defaults(func=cmd_pattern)

    # -- stats --
    pst = sub.add_parser('stats',
                         help='Print statistics about a saved FRX index')
    pst.add_argument('--index', required=True, help='Path to .frx.db file')
    pst.set_defaults(func=cmd_stats)

    # -- hgt-sim --
    ph = sub.add_parser('hgt-sim',
                        help='Run the built-in HGT simulation demo')
    ph.set_defaults(func=cmd_hgt_sim)

    # -- variant --
    pv = sub.add_parser('variant',
                        help='Hash-based variant comparison between two nodes (Mode 4)')
    pv.add_argument('--index',  required=True, help='Path to .frx.db file')
    pv.add_argument('--node1',  required=True, help='First node ID')
    pv.add_argument('--node2',  required=True, help='Second node ID')
    pv.add_argument('--threshold', type=float, default=0.5,
                    help='MinHash similarity threshold for VARIANT verdict (default: 0.5)')
    pv.add_argument('--seed', default='panindex_default_seed')
    pv.set_defaults(func=cmd_variant)

    # -- vg-import --
    pvg = sub.add_parser('vg-import',
                         help='Import a vg-generated GFA, normalize, annotate, and build index')
    pvg.add_argument('--gfa',  required=True, help='vg-generated GFA file')
    pvg.add_argument('--out',  required=True, help='Output .frx.db file path')
    pvg.add_argument('--seed', default='panindex_default_seed')
    pvg.add_argument('--root', default='PangenomeRoot')
    pvg.add_argument('--streaming', action='store_true',
                     help='Use streaming annotator for large GFA files')
    pvg.add_argument('--build-kmer-index', action='store_true', dest='build_kmer_index',
                     help='Build k-mer index after annotation')
    pvg.add_argument('--keep-vg-tags', action='store_true', dest='keep_vg_tags',
                     help='Keep vg-internal tags (LN:i:, RC:i:, etc.) in output')
    pvg.set_defaults(func=cmd_vg_import)

    return p


# ======================================================================
# Entry point
# ======================================================================

def main():
    parser = _build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
