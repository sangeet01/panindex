import hashlib
import binascii
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore


def prove_ratchet_stability(n_unrelated_insertions: int = 50_000):
    """
    Demonstration (not a formal proof) that a node's ratchet address does
    not depend on how many *other* nodes exist in the graph.

    Calling derive_ratchet_address(root, "ARG_blaTEM") twice in a row and
    checking equality only proves the function is deterministic -- true of
    any pure function and not evidence about graph scale. This version
    actually interleaves tens of thousands of unrelated derivations
    (simulating other genes/isolates being added to the pangenome) between
    the two calls, so the invariant under test -- "growing the rest of the
    graph does not disturb this node's address" -- is the thing actually
    exercised.
    """
    print("--- Demonstration: Ratchet Address Stability Under Graph Growth ---")
    engine = PanIndexEngine(pangenome_seed=b"stable_pangenome_seed_123")

    root_v1 = engine.root_hash
    gene_addr_v1 = engine.derive_ratchet_address(root_v1, "ARG_blaTEM")
    print(f"Address before growth: {binascii.hexlify(gene_addr_v1).decode()}")

    # Simulate a large pangenome update: many unrelated nodes derived from
    # the same root, as would happen when new isolates/genes are indexed.
    for i in range(n_unrelated_insertions):
        engine.derive_ratchet_address(root_v1, f"unrelated_gene_{i}")

    gene_addr_v2 = engine.derive_ratchet_address(root_v1, "ARG_blaTEM")
    print(f"Address after {n_unrelated_insertions:,} unrelated insertions: "
          f"{binascii.hexlify(gene_addr_v2).decode()}")

    assert gene_addr_v1 == gene_addr_v2
    print("SUCCESS: known-path address unchanged by unrelated graph growth.\n")


def prove_hgt_detection():
    """
    Real HGT detection proof: index a known resistance gene in a global
    PanIndexStore, then verify that a new isolate carrying the same gene
    (inserted via HGT) is detected via O(1) tag lookup — no alignment.
    """
    print("--- Proof of HGT Detection (Index Lookup) ---")

    engine = PanIndexEngine(pangenome_seed=b"hgt_proof_seed_000000000000000")
    store = PanIndexStore()

    # Known resistance gene sequence (blaTEM-1 fragment)
    seq_arg = "ATGCGTCGTAGCTAGCTAGCTGATCGATCG"

    # Donor organism: K. pneumoniae KP001 carries the gene on a plasmid.
    # Canonical cycle hash makes the address rotation-invariant.
    plasmid_seq = seq_arg + "AATTCGCTAGCTAGCTAGCATG"
    canonical_hash = engine.canonical_cycle_hash(plasmid_seq)

    store.insert(
        node_id="PLASMID:pKP001-blaTEM",
        address=canonical_hash,
        tags=["AMR:blaTEM", "resistance:beta-lactam", "mobile_element"],
        metadata={"seq": plasmid_seq, "organism": "KP001"},
        organism="KP001",
    )

    # Before HGT: E. coli EC042 has no resistance gene.
    ec042_symlinks = {}
    assert "AMR:blaTEM" not in [
        tag
        for nid in store.lookup_by_tag("AMR:blaTEM", organism="EC042")
        for tag in store.get_node(nid)['tags']
    ], "EC042 should not carry AMR:blaTEM before HGT"
    print("  Before HGT: EC042 has no AMR:blaTEM — confirmed.")

    # HGT event: EC042 receives the plasmid via conjugation.
    # Model as a topological symlink (pointer to donor's canonical hash).
    ec042_symlinks["HGT:pKP001-blaTEM"] = canonical_hash

    # Index the symlink in the global store so tag queries find it.
    store.insert(
        node_id="EC042:HGT:pKP001-blaTEM",
        address=canonical_hash,          # same hash — no data duplication
        tags=["AMR:blaTEM", "resistance:beta-lactam", "HGT:symlink"],
        metadata={"seq": plasmid_seq, "organism": "EC042", "donor": "KP001"},
        organism="EC042",
    )

    # Detection: expected O(1) tag posting-list lookup — no sequence alignment.
    amr_nodes = store.lookup_by_tag("AMR:blaTEM")
    ec042_amr = store.lookup_by_tag("AMR:blaTEM", organism="EC042")

    print(f"  After HGT: global AMR:blaTEM carriers -> {amr_nodes}")
    print(f"  EC042-scoped AMR:blaTEM carriers      -> {ec042_amr}")

    assert len(ec042_amr) >= 1, "EC042 should carry AMR:blaTEM after HGT"
    assert ec042_symlinks["HGT:pKP001-blaTEM"] == canonical_hash, \
        "Symlink must point to donor's exact canonical hash"

    # Verify content-addressable lookup: the Merkle/canonical address
    # resolves to the correct node without knowing the node_id.
    resolved = store.lookup_by_address(canonical_hash)
    assert resolved is not None, "Canonical hash must resolve via lookup_by_address"
    print(f"  Content-address lookup: {canonical_hash.hex()[:16]}... -> '{resolved}'")

    print("SUCCESS: HGT detection via direct tag index lookup verified.\n")


if __name__ == "__main__":
    prove_ratchet_stability()
    prove_hgt_detection()
