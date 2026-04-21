import bisect
from typing import List, Optional, Dict, Any, Tuple


class PanIndexStore:
    """
    PanIndex B-Tree-style Index.

    Two internal structures:
    1. Address Index   -> sorted list of (address_hex, node_id) pairs.
       Provides O(log K) lookup by ratchet address.
    2. Tag Index       -> dict of tag_string -> [node_id, ...].
       Provides O(1) lookup by Anubandha tag.
    """

    def __init__(self):
        # Sorted list of address strings (hex), parallel list of node_ids
        self._addr_keys: List[str] = []
        self._addr_vals: List[str] = []

        # Tag -> list of node_ids
        self._tag_index: Dict[str, List[str]] = {}

        # node_id -> full metadata dict
        self._node_store: Dict[str, Dict[str, Any]] = {}

    # ------------------------------------------------------------------
    # Write operations
    # ------------------------------------------------------------------

    def insert(self, node_id: str, address: bytes, tags: List[str], metadata: Dict[str, Any] = None):
        """
        Insert a node and its PanIndex address into the index.

        Args:
            node_id   : GFA segment ID (e.g. '1', '42').
            address   : 32-byte PanIndex address (bytes).
            tags      : List of Anubandha tag strings (e.g. ['AMR:blaTEM', 'Chr4']).
            metadata  : Any additional data to store (sequence, derivation path, etc.).
        """
        addr_hex = address.hex()

        # Sorted insert into address index (B-Tree leaf approximation)
        pos = bisect.bisect_left(self._addr_keys, addr_hex)
        if pos < len(self._addr_keys) and self._addr_keys[pos] == addr_hex:
            # Address already exists - update
            self._addr_vals[pos] = node_id
        else:
            self._addr_keys.insert(pos, addr_hex)
            self._addr_vals.insert(pos, node_id)

        # Tag index
        for tag in tags:
            self._tag_index.setdefault(tag, [])
            if node_id not in self._tag_index[tag]:
                self._tag_index[tag].append(node_id)

        # Node store
        self._node_store[node_id] = {
            'address': addr_hex,
            'tags': tags,
            'metadata': metadata or {},
        }

    # ------------------------------------------------------------------
    # Read operations
    # ------------------------------------------------------------------

    def lookup_by_address(self, address: bytes) -> Optional[str]:
        """
        O(log K) lookup of a node_id by its exact PanIndex address.
        Returns node_id or None if not found.
        """
        addr_hex = address.hex()
        pos = bisect.bisect_left(self._addr_keys, addr_hex)
        if pos < len(self._addr_keys) and self._addr_keys[pos] == addr_hex:
            return self._addr_vals[pos]
        return None

    def lookup_by_tag(self, tag: str) -> List[str]:
        """
        O(1) lookup of all node_ids carrying a given Anubandha tag.
        """
        return self._tag_index.get(tag, [])

    def range_lookup(self, address_start: bytes, address_end: bytes) -> List[Tuple[str, str]]:
        """
        O(log K + R) range scan between two addresses (lexicographic).
        Useful for zoomed-out region queries.
        Returns list of (address_hex, node_id) pairs.
        """
        lo = bisect.bisect_left(self._addr_keys, address_start.hex())
        hi = bisect.bisect_right(self._addr_keys, address_end.hex())
        return list(zip(self._addr_keys[lo:hi], self._addr_vals[lo:hi]))

    def get_node(self, node_id: str) -> Optional[Dict[str, Any]]:
        """Retrieve full node metadata."""
        return self._node_store.get(node_id)

    def all_nodes(self) -> List[str]:
        return list(self._node_store.keys())

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    def stats(self) -> Dict[str, int]:
        return {
            'total_nodes': len(self._node_store),
            'total_tags': sum(len(v) for v in self._tag_index.values()),
            'unique_tags': len(self._tag_index),
        }

    def __len__(self):
        return len(self._node_store)

    def __repr__(self):
        s = self.stats()
        return (
            f"PanIndexStore("
            f"nodes={s['total_nodes']}, "
            f"unique_tags={s['unique_tags']}, "
            f"tag_entries={s['total_tags']})"
        )


if __name__ == "__main__":
    import hashlib
    from engine import PanIndexEngine

    engine = PanIndexEngine(pangenome_seed=b"demo_seed_panindex_0123456789ab")
    store = PanIndexStore()

    # Simulate indexing three nodes
    nodes = [
        ("1", "ACTG",  ["upstream"]),
        ("2", "A",     ["SNP", "AMR:blaTEM"]),
        ("3", "T",     ["SNP"]),
        ("4", "GGGC",  ["downstream"]),
    ]

    parent = engine.root_hash
    for node_id, seq, tags in nodes:
        addr = engine.derive_ratchet_address(parent, node_id)
        final_addr = engine.compute_node_address(seq, addr)
        store.insert(node_id, final_addr, tags, {'seq': seq})
        parent = final_addr  # chain for demo purposes

    print(store)

    # Lookup by address
    addr_node2 = store.get_node("2")['address']
    found = store.lookup_by_address(bytes.fromhex(addr_node2))
    print(f"Lookup by address -> Node: {found}")

    # Lookup by tag
    amr_nodes = store.lookup_by_tag("AMR:blaTEM")
    print(f"Lookup by tag 'AMR:blaTEM' -> Nodes: {amr_nodes}")
    assert "2" in amr_nodes
    print("SUCCESS: Index lookups verified.")
