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


class VariantDiffResult:
    """
    Result of a Mode 4 hash-based variant comparison.

    Verdicts
    --------
    IDENTICAL  : Addresses match exactly. Same sequence + same context.
    VARIANT    : Addresses differ. Same organism. MinHash similarity >= 0.5.
    NOVEL      : Addresses differ. Different organism or similarity < 0.5.
    NOT_FOUND  : One or both node IDs absent from the index.
    """

    VERDICTS = ('IDENTICAL', 'VARIANT', 'NOVEL', 'NOT_FOUND')

    def __init__(self, node_id_a: str, node_id_b: str, verdict: str,
                 similarity: float, addr_a: Optional[str],
                 addr_b: Optional[str], same_organism: bool):
        self.node_id_a = node_id_a
        self.node_id_b = node_id_b
        self.verdict = verdict
        self.similarity = similarity
        self.addr_a = addr_a
        self.addr_b = addr_b
        self.same_organism = same_organism

    def print(self):
        print(f"\n[Variant Diff] '{self.node_id_a}' vs '{self.node_id_b}'")
        print(f"  Verdict       : {self.verdict}")
        if self.verdict == 'NOT_FOUND':
            return
        print(f"  Similarity    : {self.similarity:.4f}")
        print(f"  Same organism : {self.same_organism}")
        if self.addr_a:
            print(f"  Addr A        : {self.addr_a[:24]}...")
        if self.addr_b:
            print(f"  Addr B        : {self.addr_b[:24]}...")


class PanIndexQuery:
    """
    Phase 4 - Unified Query Interface.

    Four modes, all work against a populated PanIndexStore:

    Mode 1 - Ratchet Path Query  (O(d) for path depth d):
        Input : "Root/Chr4/BRCA1/VarA"
        Action: Walk the slash-separated path, chaining HKDF derivations.
                The final derived address is looked up in the index.
        Use   : Jump to any node in the pangenome without graph traversal.

    Mode 2 - Tag Query  (expected O(1) posting-list lookup + results):
        Input : "AMR:blaTEM"
        Action: Direct tag index lookup using store.lookup_by_tag().
        Use   : Find all nodes carrying a known Anubandha annotation.

    Mode 3 - LSH Similarity Query  (O(K) over indexed nodes):
        Input : A raw sequence string, e.g. "ATGCGTCGTA..."
        Action: Compute MinHash signature of query. Compare against
                stored node sequences using Hamming distance.
                Returns nodes within the similarity threshold.
        Use   : Find structurally similar genomic regions without alignment.

    Mode 4 - Hash-Based Variant Diff  (O(1) address compare + O(k) MinHash):
        Input : Two node IDs, e.g. "gene_blaTEM", "VarA"
        Action: Compare ratchet addresses. If different, compute MinHash
                similarity and check organism membership.
        Use   : Determine IDENTICAL / VARIANT / NOVEL relationship without
                sequence alignment. Unique to FRX.
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

        Uses the standard single-hash + linear-congruential permutation trick:
        for each k-mer, compute one SHA-256 hash, then derive num_hashes
        permuted values via h_i(x) = (a_i * h(x) + b_i) mod P.
        This is O(|kmers| + num_hashes) rather than O(|kmers| * num_hashes).

        Permutation coefficients are seeded deterministically so signatures
        are reproducible across processes.
        """
        n = len(sequence)
        if n < k:
            return [0] * num_hashes

        # Mersenne prime for the hash universe
        _P = (1 << 61) - 1

        # Deterministic permutation coefficients (seeded, never 0)
        import struct
        _a = []
        _b = []
        for i in range(num_hashes):
            seed_a = hashlib.sha256(struct.pack('<QQ', i, 0xDEADBEEF)).digest()
            seed_b = hashlib.sha256(struct.pack('<QQ', i, 0xCAFEBABE)).digest()
            _a.append((int.from_bytes(seed_a[:8], 'little') % (_P - 1)) + 1)
            _b.append(int.from_bytes(seed_b[:8], 'little') % _P)

        # One SHA-256 per k-mer
        kmer_hashes = [
            int(hashlib.sha256(sequence[i:i + k].encode()).hexdigest(), 16) % _P
            for i in range(n - k + 1)
        ]

        if not kmer_hashes:
            return [0] * num_hashes

        # For each permutation, take the minimum over all k-mer hashes
        signature = [
            min((_a[i] * h + _b[i]) % _P for h in kmer_hashes)
            for i in range(num_hashes)
        ]
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
            if not node or not node.get('metadata'):
                continue
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
    # Mode 4: Hash-Based Variant Diff
    # ------------------------------------------------------------------

    def query_variant_diff(
        self,
        node_id_a: str,
        node_id_b: str,
        similarity_threshold: float = 0.5,
    ) -> VariantDiffResult:
        """
        Compare two indexed nodes and classify their relationship.

        Algorithm:
        1. Look up both nodes. If either is missing -> NOT_FOUND.
        2. Compare ratchet addresses directly. Equal -> IDENTICAL.
        3. Compute MinHash similarity of their sequences.
        4. Check organism membership.
        5. Similarity >= threshold AND same organism -> VARIANT.
           Otherwise -> NOVEL.

        Address comparisons are O(1); MinHash comparison is O(k), where k is
        the signature length. Query and result enumeration costs still apply.
        No sequence alignment is performed.

        Args:
            node_id_a          : ID of the first node.
            node_id_b          : ID of the second node.
            similarity_threshold: MinHash Jaccard threshold for VARIANT. Default 0.5.

        Returns:
            VariantDiffResult with verdict, similarity score, and addresses.
        """
        node_a = self.store.get_node(node_id_a)
        node_b = self.store.get_node(node_id_b)

        if node_a is None or node_b is None:
            missing = node_id_a if node_a is None else node_id_b
            print(f"  [Variant Diff] Node not found: '{missing}'")
            return VariantDiffResult(
                node_id_a, node_id_b,
                verdict='NOT_FOUND',
                similarity=0.0,
                addr_a=node_a['address'] if node_a else None,
                addr_b=node_b['address'] if node_b else None,
                same_organism=False,
            )

        addr_a = node_a['address']
        addr_b = node_b['address']
        org_a = node_a.get('organism', '')
        org_b = node_b.get('organism', '')
        same_organism = (org_a == org_b)

        # Step 2: Exact address match -> IDENTICAL
        if addr_a == addr_b:
            return VariantDiffResult(
                node_id_a, node_id_b,
                verdict='IDENTICAL',
                similarity=1.0,
                addr_a=addr_a,
                addr_b=addr_b,
                same_organism=same_organism,
            )

        # Step 3+4: MinHash similarity + organism check
        seq_a = node_a['metadata'].get('seq', '')
        seq_b = node_b['metadata'].get('seq', '')
        sig_a = self._minhash_signature(seq_a)
        sig_b = self._minhash_signature(seq_b)
        similarity = self._hamming_similarity(sig_a, sig_b)

        if similarity >= similarity_threshold and same_organism:
            verdict = 'VARIANT'
        else:
            verdict = 'NOVEL'

        return VariantDiffResult(
            node_id_a, node_id_b,
            verdict=verdict,
            similarity=similarity,
            addr_a=addr_a,
            addr_b=addr_b,
            same_organism=same_organism,
        )


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
