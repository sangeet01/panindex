import bisect
from typing import List, Optional, Dict, Any, Tuple, Union


class PanIndexStore:
    """
    PanIndex B-Tree-style Index.

    Two internal structures:
    1. Address Index   -> sorted list of (address_hex, node_id) pairs.
       Provides O(log K) lookup by ratchet address.
     2. Tag Index       -> dict of tag_string -> [node_id, ...].
         Provides expected O(1) lookup of a tag posting list.
    """

    def __init__(self):
        # Sorted list of address strings (hex), parallel list of node_ids
        self._addr_keys: List[str] = []
        self._addr_vals: List[str] = []

        # Tag -> list of node_ids
        self._tag_index: Dict[str, List[str]] = {}

        # Content ID -> all placements containing that sequence and strand
        self._content_index: Dict[str, List[str]] = {}

        # node_id -> full metadata dict
        self._node_store: Dict[str, Dict[str, Any]] = {}

    # ------------------------------------------------------------------
    # Write operations
    # ------------------------------------------------------------------

    def insert(
        self,
        node_id: str,
        address: bytes,
        tags: List[str],
        metadata: Dict[str, Any] = None,
        organism: str = "",
        merkle_addr: str = "",
        content_id: str = "",
        index_alias: bool = False,
    ):
        """
        Insert a node and its PanIndex address into the index.

        Args:
            node_id    : GFA segment ID (e.g. '1', '42').
            address    : 32-byte ratchet path address (bytes).
            tags       : List of Anubandha tag strings.
            metadata   : Any additional data (sequence, derivation path, etc.).
            organism   : Organism/strain identifier for namespace scoping.
            merkle_addr: Hex string of the content-derived Merkle address.
                         When provided, the node is also indexed under this
                         address so lookup_by_address() works for both the
                         ratchet path address and the Merkle address.
            content_id : Hex string identifying sequence and strand content.
                         Multiple node placements may share this identity.
        """
        addr_hex = address.hex()

        # Sorted insert into address index (B-Tree leaf approximation)
        pos = bisect.bisect_left(self._addr_keys, addr_hex)
        if pos < len(self._addr_keys) and self._addr_keys[pos] == addr_hex:
            self._addr_vals[pos] = node_id
        else:
            self._addr_keys.insert(pos, addr_hex)
            self._addr_vals.insert(pos, node_id)

        # Also index by Merkle address when it differs from the ratchet address
        if merkle_addr and merkle_addr != addr_hex:
            mpos = bisect.bisect_left(self._addr_keys, merkle_addr)
            if not (mpos < len(self._addr_keys) and self._addr_keys[mpos] == merkle_addr):
                self._addr_keys.insert(mpos, merkle_addr)
                self._addr_vals.insert(mpos, node_id)

        is_alias = metadata.get('is_alias', False) if metadata else False

        # Tag index
        if not is_alias:
            for tag in tags:
                self._tag_index.setdefault(tag, [])
                if node_id not in self._tag_index[tag]:
                    self._tag_index[tag].append(node_id)

        if content_id and not is_alias:
            self._content_index.setdefault(content_id, [])
            if node_id not in self._content_index[content_id]:
                self._content_index[content_id].append(node_id)

        # Node store
        self._node_store[node_id] = {
            'address': addr_hex,
            'merkle_addr': merkle_addr,
            'content_id': content_id,
            'tags': tags,
            'metadata': metadata or {},
            'organism': organism,
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

    def lookup_by_tag(self, tag: str, organism: Optional[str] = None) -> List[str]:
        """
        Expected O(1) lookup of the posting list for an Anubandha tag;
        returning all matches costs O(R) for R results.

        Args:
            tag      : Anubandha tag string (e.g. 'AMR:blaTEM').
            organism : If specified, return only nodes from this organism.
                       None (default) returns nodes from all organisms.

        Returns:
            List of matching node_ids.
        """
        nodes = self._tag_index.get(tag, [])
        if organism is None:
            return nodes
        return [
            n for n in nodes
            if self._node_store.get(n, {}).get('organism', '') == organism
        ]

    def lookup_by_content_id(self, content_id: Union[bytes, str]) -> List[str]:
        """Return every placement carrying a sequence content identity."""
        key = content_id.hex() if isinstance(content_id, bytes) else content_id
        return list(self._content_index.get(key, []))

    def range_lookup(self, address_start: bytes, address_end: bytes) -> List[Tuple[str, str]]:
        """
        O(log K + R) range scan between two addresses (lexicographic).
        Useful for zoomed-out region queries.
        Returns list of (address_hex, node_id) pairs.
        """
        lo = bisect.bisect_left(self._addr_keys, address_start.hex())
        hi = bisect.bisect_right(self._addr_keys, address_end.hex())
        return list(zip(self._addr_keys[lo:hi], self._addr_vals[lo:hi]))

    def get_node(
        self,
        node_id: str,
        organism: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        """
        Retrieve full node metadata.

        Args:
            node_id  : Segment ID to look up.
            organism : If specified, return None if the node's organism
                       does not match. None (default) returns any match.

        Returns:
            Node record dict or None if not found / organism mismatch.
        """
        record = self._node_store.get(node_id)
        if record is None:
            return None
        if organism is not None and record.get('organism', '') != organism:
            return None
        return record

    def all_nodes(self, organism: Optional[str] = None) -> List[str]:
        """
        Return all node_ids, optionally scoped to a specific organism.
        """
        if organism is None:
            return list(self._node_store.keys())
        return [
            nid for nid, rec in self._node_store.items()
            if rec.get('organism', '') == organism
        ]

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: str) -> int:
        """
        Persist this store to a SQLite database file.

        Args:
            path: Output file path (e.g. 'pangenome.frx.db').
                  Created or overwritten.

        Returns:
            Number of nodes written.
        """
        from persistence import save_store
        n = save_store(self, path)
        return n

    @classmethod
    def load(cls, path: str) -> 'PanIndexStore':
        """
        Load a PanIndexStore from a SQLite database file saved by save().

        Args:
            path: Path to a .db file written by save().

        Returns:
            Populated PanIndexStore ready for all query operations.

        Raises:
            FileNotFoundError if path does not exist.
        """
        from persistence import load_store
        return load_store(path)

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
