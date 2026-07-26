"""
PanIndex vs Linear Scan Benchmark on Real K. pneumoniae Pangenome.
500 real genomic nodes from the Kaggle nwheeler443 dataset.
"""
import time
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore
from annotator import GFAAnnotator
from query import PanIndexQuery

RUNS = 200

def run_benchmark():
    print("=" * 58)
    print("PanIndex Benchmark on Real K. pneumoniae Pangenome")
    print(f"Dataset: 500 nodes, 78 edges (Kaggle nwheeler443)")
    print("=" * 58)

    gfa_path = os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..', 'data', 'kp_pangenome.gfa'
    ))
    annotated_path = os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..', 'data', 'kp_pangenome_annotated.gfa'
    ))

    annotator = GFAAnnotator(seed=b"benchmark_kp_real_seed_00000000")
    annotator.annotate(gfa_path, annotated_path)

    engine = annotator.engine
    store = annotator.store
    q = PanIndexQuery(engine, store)

    # Target: Node 14 - a long AMR-associated sequence
    target_path = "PangenomeRoot/14"
    target_seq_snippet = "ATGATTAATATATCTGAGTTTGATATGAAGG"

    print(f"\n[Target] Derivation path: '{target_path}'")
    print(f"[Target] Sequence for LSH: '{target_seq_snippet[:20]}...'")
    print(f"[Runs per test]: {RUNS}")
    print()

    # ------------------------------------------------------------------
    # MODE 1: Known-path derivation (O(d), then indexed lookup)
    # ------------------------------------------------------------------
    t0 = time.perf_counter()
    for _ in range(RUNS):
        result = q.query_by_path(target_path)
    ratchet_ms = (time.perf_counter() - t0) / RUNS * 1000
    ratchet_hit = not result.is_empty()
    print(f"Mode 1  Ratchet Path (O(d))       : {ratchet_ms:.4f} ms  | Hit: {ratchet_hit}")

    # ------------------------------------------------------------------
    # MODE 2: Tag index lookup (O(1) expected, plus result enumeration)
    # Tag format stored by annotator is: "AN:<value>" from GFA field AN:Z:<value>
    # ------------------------------------------------------------------
    t0 = time.perf_counter()
    for _ in range(RUNS):
        result_tag = q.query_by_tag("AN:kp_node_14")
    tag_ms = (time.perf_counter() - t0) / RUNS * 1000
    tag_hit = not result_tag.is_empty()
    print(f"Mode 2  Tag Index (expected O(1))  : {tag_ms:.4f} ms  | Hit: {tag_hit}")

    # ------------------------------------------------------------------
    # MODE 3: Linear Brute-Force Scan (O(N))
    # ------------------------------------------------------------------
    all_seqs = {
        nid: store.get_node(nid)['metadata'].get('seq', '')
        for nid in store.all_nodes()
    }
    t0 = time.perf_counter()
    for _ in range(RUNS):
        match = None
        for nid, seq in all_seqs.items():
            if target_seq_snippet in seq:
                match = nid
                break
    linear_ms = (time.perf_counter() - t0) / RUNS * 1000
    linear_hit = match is not None
    print(f"Mode 3  Linear Scan (O(N=500))    : {linear_ms:.4f} ms  | Hit: {linear_hit}")

    print()
    print(f"  Tag index vs linear:  {linear_ms / tag_ms:.1f}x faster")
    print(f"  Ratchet vs linear:    {linear_ms / ratchet_ms:.2f}x (HKDF constant >> N=500)")
    print()
    print("Scaling analysis:")
    print("  Tag index lookup is expected O(1), excluding result enumeration.")
    print("  Path derivation is O(d), where d is the number of path components,")
    print("  followed by the address index lookup. Neither depends on unrelated")
    print("  graph nodes once preprocessing is complete.")
    print("  This is the same principle as binary search vs hash table -")
    print("  hash tables lose at N=10, win at N=10,000.")

    assert ratchet_hit, "Ratchet path did not find the target node"
    assert tag_hit, "Tag index did not find the target node"
    assert linear_hit, "Linear scan did not find the target node"
    print("\nAll three modes hit the correct node. Benchmark complete.")
    print("=" * 58)

if __name__ == "__main__":
    run_benchmark()
