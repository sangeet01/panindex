import sys
import os
import hashlib
import argparse
from typing import List, Optional, Dict, Any

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore


class QueryResult:
    """Structured result returned by all query modes."""

    def __init__(self, query_mode: str, query_input: str,
                 matched_nodes: List[str], store: PanIndexStore,
                 derived_address: Optional[bytes] = None):
        self.query_mode = query_mode
        self.query_input = query_input
        self.matched_nodes = matched_nodes
        self.store = store
        self.derived_address = derived_address

    def is_empty(self) -> bool:
        return len(self.matched_nodes) == 0

    def print(self):
        print(f"\n[{self.query_mode}] Query: '{self.query_input}'")
        if self.derived_address:
            print(f"  Derived address : {self.derived_address.hex()[:24]}...")
        if self.is_empty():
            print("  Result          : No match found.")
            return
        print(f"  Matched nodes   : {len(self.matched_nodes)}")
        for nid in self.matched_nodes:
            node = self.store.get_node(nid)
            if node:
                seq = node['metadata'].get('seq', 'N/A')
                tags = ', '.join(node['tags']) if node['tags'] else 'none'
                path = node['metadata'].get('derivation_path', 'N/A')
                print(f"    -> Node '{nid}'")
                print(f"       Sequence path : {path}")
                print(f"       Sequence      : {seq}")
                print(f"       Tags          : {tags}")
                print(f"       Address       : {node['address'][:24]}...")


class PanIndexQuery:
    """
    Phase 4 - Unified Query Interface.

    Three modes, all work against a populated PanIndexStore:

    Mode 1 - Ratchet Path Query  (O(1) per level):
        Input : "Root/Chr4/BRCA1/VarA"
        Action: Walk the slash-separated path, chaining HKDF derivations.
                The final derived address is looked up in the index.
        Use   : Jump to any node in the pangenome without graph traversal.

    Mode 2 - Tag Query  (O(1)):
        Input : "AMR:blaTEM"
        Action: Direct tag index lookup using store.lookup_by_tag().
        Use   : Find all nodes carrying a known Anubandha annotation.

    Mode 3 - LSH Similarity Query  (O(K) over indexed nodes):
        Input : A raw sequence string, e.g. "ATGCGTCGTA..."
        Action: Compute MinHash signature of query. Compare against
                stored node sequences using Hamming distance.
                Returns nodes within the similarity threshold.
        Use   : Find structurally similar genomic regions without alignment.
    """

    def __init__(self, engine: PanIndexEngine, store: PanIndexStore,
                 root_label: str = "Root"):
        self.engine = engine
        self.store = store
        self.root_label = root_label

    # ------------------------------------------------------------------
    # Mode 1: Ratchet Path Query
    # ------------------------------------------------------------------

    def query_by_path(self, path: str) -> QueryResult:
        """
        Derive an address from a hierarchical path string and look it up.

        Args:
            path: Slash-separated derivation path, must start with root label.
                  Example: "Root/KP001/Chr1/gene_blaTEM/VarA"
        """
        parts = [p.strip() for p in path.split('/') if p.strip()]

        if not parts:
            return QueryResult("RatchetPath", path, [], self.store)

        # Strip root label if present as first component
        if parts[0] == self.root_label:
            parts = parts[1:]

        # Chain HKDF derivations down the path
        current_hash = self.engine.root_hash
        derivation_log = [self.root_label]

        for part in parts:
            current_hash = self.engine.derive_ratchet_address(current_hash, part)
            derivation_log.append(part)

        # The final hash is a pure ratchet address.
        # Compute the full node address: compute_node_address needs sequence.
        # Since we don't have the sequence at query time, we look up by the
        # raw derived address (Layer B address from the two-layer architecture).
        # The index stores Layer B addresses directly for path-queryable nodes.
        node_id = self.store.lookup_by_address(current_hash)
        matched = [node_id] if node_id else []

        return QueryResult(
            "RatchetPath",
            path,
            matched,
            self.store,
            derived_address=current_hash
        )

    # ------------------------------------------------------------------
    # Mode 2: Tag Query
    # ------------------------------------------------------------------

    def query_by_tag(self, tag: str) -> QueryResult:
        """
        Return all nodes carrying the given Anubandha tag.

        Args:
            tag: Tag string (e.g. "AMR:blaTEM", "upstream", "SNP").
        """
        matched = self.store.lookup_by_tag(tag)
        return QueryResult("TagQuery", tag, matched, self.store)

    # ------------------------------------------------------------------
    # Mode 3: LSH / MinHash Similarity Query
    # ------------------------------------------------------------------

    def _minhash_signature(self, sequence: str, k: int = 4,
                           num_hashes: int = 64) -> List[int]:
        """
        Compute a MinHash signature for a k-mer set of the sequence.

        Each hash function h_i is simulated by seeding SHA-256 with a
        different salt (the hash index i), producing a minimal hash over
        all k-mers. This is the standard bottom-k MinHash approach.
        """
        n = len(sequence)
        if n < k:
            return [0] * num_hashes

        kmers = set(sequence[i:i+k] for i in range(n - k + 1))
        signature = []

        for i in range(num_hashes):
            min_val = None
            for kmer in kmers:
                h = int(hashlib.sha256(
                    i.to_bytes(4, 'little') + kmer.encode()
                ).hexdigest(), 16)
                if min_val is None or h < min_val:
                    min_val = h
            signature.append(min_val if min_val is not None else 0)

        return signature

    def _hamming_similarity(self, sig_a: List[int], sig_b: List[int]) -> float:
        """
        Jaccard similarity estimate via MinHash Hamming agreement.
        Returns a float in [0.0, 1.0]; higher means more similar.
        """
        if not sig_a or not sig_b or len(sig_a) != len(sig_b):
            return 0.0
        matches = sum(1 for a, b in zip(sig_a, sig_b) if a == b)
        return matches / len(sig_a)

    def query_by_similarity(self, sequence: str,
                            threshold: float = 0.5) -> QueryResult:
        """
        Find all indexed nodes whose sequence is similar to the query.

        Args:
            sequence  : Raw nucleotide sequence to search for.
            threshold : Minimum Jaccard similarity [0.0, 1.0]. Default 0.5.
        """
        query_sig = self._minhash_signature(sequence)
        matched = []

        for node_id in self.store.all_nodes():
            node = self.store.get_node(node_id)
            node_seq = node['metadata'].get('seq', '')
            if not node_seq:
                continue
            node_sig = self._minhash_signature(node_seq)
            sim = self._hamming_similarity(query_sig, node_sig)
            if sim >= threshold:
                matched.append((sim, node_id))

        matched.sort(reverse=True)
        matched_ids = [nid for _, nid in matched]

        return QueryResult("LSHSimilarity", sequence[:20] + "...", matched_ids, self.store)


# ------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------

def build_demo_state():
    """Build a demo engine + store so the CLI can run standalone."""
    engine = PanIndexEngine(pangenome_seed=b"panindex_query_demo_seed_0000000")
    store = PanIndexStore()

    nodes = [
        ("Chr1",        "ATGCATGCATGCATGC", ["core_chromosome"],           "Root/Chr1"),
        ("gene_blaTEM", "ATGCGTCGTAGCTAGC", ["AMR:blaTEM", "mobile_element"], "Root/Chr1/gene_blaTEM"),
        ("VarA",        "ATGCGTCGTAGCTAGT", ["AMR:blaTEM", "variant"],     "Root/Chr1/gene_blaTEM/VarA"),
        ("Chr2",        "GCTAGCTAGCTAGCTA", ["core_chromosome"],           "Root/Chr2"),
        ("gene_rpoB",   "TTAATTAATTAATTAA", ["housekeeping"],               "Root/Chr2/gene_rpoB"),
    ]

    parent = engine.root_hash
    for node_id, seq, tags, deriv_path in nodes:
        # Build chained ratchet address following the derivation path
        parts = [p for p in deriv_path.split('/') if p and p != 'Root']
        current = engine.root_hash
        for part in parts:
            current = engine.derive_ratchet_address(current, part)

        store.insert(node_id, current, tags,
                     {'seq': seq, 'derivation_path': deriv_path})

    return engine, store


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PanIndex Query Interface",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python src/query.py --gfa data/real_kp_annotated.gfa --path "PangenomeRoot/10"
  python src/query.py --gfa data/real_kp_annotated.gfa --tag "AN:sequence_10"
  python src/query.py --gfa data/real_kp_annotated.gfa --sequence "TAAAAAAGCCAT..."
        """
    )
    parser.add_argument('--gfa', type=str, help='Annotated GFA file to load')
    parser.add_argument('--path', type=str,
                        help='Ratchet derivation path (e.g. "Root/Chr1/gene_blaTEM")')
    parser.add_argument('--tag', type=str,
                        help='Anubandha tag to search (e.g. "AMR:blaTEM")')
    parser.add_argument('--sequence', type=str,
                        help='Nucleotide sequence for similarity search')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='MinHash similarity threshold (default: 0.5)')
    args = parser.parse_args()

    if args.gfa:
        engine = PanIndexEngine()
        store = PanIndexStore()
        # Scan GFA and populate store
        with open(args.gfa, 'r') as f:
            for line in f:
                if line.startswith('S'):
                    parts = line.strip().split('\t')
                    node_id = parts[1]
                    seq = parts[2]
                    tags = parts[3:]
                    # Extract address from AN:Z: tag
                    addr = None
                    anubandha = []
                    for t in tags:
                        if t.startswith('AN:Z:'):
                            if len(t) > 37: # Likely the hex address
                                try:
                                    addr = bytes.fromhex(t[5:])
                                except:
                                    anubandha.append(t[5:])
                            else:
                                anubandha.append(t[5:])
                        elif t.startswith('PA:Z:'):
                            # Store path in metadata
                            pass
                    
                    if addr:
                        store.insert(node_id, addr, anubandha, {
                            'seq': seq,
                            'derivation_path': next((t[5:] for t in tags if t.startswith('PA:Z:')), '')
                        })
    else:
        engine, store = build_demo_state()

    q = PanIndexQuery(engine, store)

    if not any([args.path, args.tag, args.sequence]):
        if not args.gfa:
            print("PanIndex Query Demo (Synthetic) - running all three modes:\n")
            q.query_by_path("Root/Chr1/gene_blaTEM").print()
            q.query_by_tag("AMR:blaTEM").print()
            q.query_by_similarity("ATGCGTCGTAGCTAGC", threshold=0.4).print()
        else:
            print(f"GFA '{args.gfa}' loaded. Use --path, --tag, or --sequence to query.")
    else:
        if args.path:
            q.query_by_path(args.path).print()
        if args.tag:
            q.query_by_tag(args.tag).print()
        if args.sequence:
            q.query_by_similarity(args.sequence, args.threshold).print()
