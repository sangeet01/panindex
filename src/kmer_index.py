"""
FRX K-mer Inverted Index.

Replaces the O(N * L) linear pattern scan in SubsequenceQuery with an
O(|pattern| / k) average-case search via k-mer seeding.

Algorithm
---------
Build: one pass over all sequences in a PanIndexStore.
       Extract every k-mer -> record (node_id, position) pairs.

Search (len >= k):
  1. Seed: look up the first k-mer of the query pattern.
  2. Verify: check if the full pattern matches at each candidate position.
  3. Return verified hits in the same dict format as search_pattern().

Search (len < k): fall back to linear scan (rare in practice since
  real genomic queries are almost always >= 12 bp).

Persistence: stored in the same .frx.db SQLite file as the main index,
  in a 'kmer_index' table. Old .db files without this table load fine
  (is_present() returns False, caller falls back to linear scan).

Default k=12: 4^12 = 16.7M possible k-mers, low collision rate at scale.
"""

import sqlite3
import os
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING

if TYPE_CHECKING:
    from index import PanIndexStore
    from engine import PanIndexEngine


_DDL = """
CREATE TABLE IF NOT EXISTS kmer_index (
    kmer    TEXT NOT NULL,
    node_id TEXT NOT NULL,
    pos     INTEGER NOT NULL,
    PRIMARY KEY (kmer, node_id, pos)
);
CREATE INDEX IF NOT EXISTS idx_kmer ON kmer_index(kmer);
"""

_KMER_META_TABLE = """
CREATE TABLE IF NOT EXISTS kmer_meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""


class KmerIndex:
    """
    K-mer inverted index for fast nucleotide pattern search.

    Build from a PanIndexStore, save to .frx.db, load back for queries.
    Transparent fallback to linear scan for short patterns.

    Usage (build + save):
        ki = KmerIndex.build(store, k=12)
        ki.save("pangenome.frx.db")

    Usage (load + search):
        ki = KmerIndex.load("pangenome.frx.db")
        hits = ki.search("ATGCGTCGTA", store, engine)
    """

    def __init__(self, k: int = 12):
        self.k = k
        # kmer_string -> [(node_id, position), ...]
        self._index: Dict[str, List[Tuple[str, int]]] = {}

    # ------------------------------------------------------------------
    # Build
    # ------------------------------------------------------------------

    @classmethod
    def build(cls, store: "PanIndexStore", k: int = 12) -> "KmerIndex":
        """
        Build a KmerIndex from a populated PanIndexStore.

        One pass over all sequences: extracts every k-mer and records
        (node_id, position) for each occurrence.

        Args:
            store : Populated PanIndexStore.
            k     : K-mer length. Default 12.

        Returns:
            Populated KmerIndex ready for save() or search().
        """
        ki = cls(k=k)
        idx = ki._index

        for node_id in store.all_nodes():
            node = store.get_node(node_id)
            seq = node['metadata'].get('seq', '')
            if len(seq) < k:
                continue
            for pos in range(len(seq) - k + 1):
                kmer = seq[pos:pos + k]
                if kmer not in idx:
                    idx[kmer] = []
                idx[kmer].append((node_id, pos))

        return ki

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, db_path: str) -> int:
        """
        Write the kmer index into an existing FRX SQLite .db file.

        Creates the kmer_index table if absent.
        Clears and rewrites all kmer rows (idempotent).

        Args:
            db_path : Path to an existing .frx.db file.

        Returns:
            Number of (kmer, node_id, pos) rows written.
        """
        conn = sqlite3.connect(db_path)
        try:
            conn.executescript(_DDL)
            conn.executescript(_KMER_META_TABLE)
            conn.execute("DELETE FROM kmer_index")
            conn.execute("DELETE FROM kmer_meta")

            rows = []
            for kmer, entries in self._index.items():
                for node_id, pos in entries:
                    rows.append((kmer, node_id, pos))

            conn.executemany(
                "INSERT OR IGNORE INTO kmer_index (kmer, node_id, pos) "
                "VALUES (?, ?, ?)",
                rows
            )
            conn.execute(
                "INSERT INTO kmer_meta (key, value) VALUES ('k', ?)",
                (str(self.k),)
            )
            conn.commit()
            return len(rows)
        finally:
            conn.close()

    @classmethod
    def load(cls, db_path: str, k: int = 12) -> "KmerIndex":
        """
        Load a KmerIndex from a .frx.db file.

        Args:
            db_path : Path to a .db file with a kmer_index table.
            k       : K-mer length used at build time (read from db if stored).

        Returns:
            Populated KmerIndex.

        Raises:
            FileNotFoundError if db_path does not exist.
        """
        if not os.path.isfile(db_path):
            raise FileNotFoundError(f"FRX index not found: {db_path}")

        conn = sqlite3.connect(db_path)
        ki = None
        try:
            # Read k from metadata if stored
            try:
                row = conn.execute(
                    "SELECT value FROM kmer_meta WHERE key='k'"
                ).fetchone()
                if row:
                    k = int(row[0])
            except Exception:
                pass

            ki = cls(k=k)
            cursor = conn.execute(
                "SELECT kmer, node_id, pos FROM kmer_index"
            )
            for kmer, node_id, pos in cursor:
                if kmer not in ki._index:
                    ki._index[kmer] = []
                ki._index[kmer].append((node_id, pos))
        finally:
            conn.close()

        return ki

    @staticmethod
    def is_present(db_path: str) -> bool:
        """
        Return True if the .db file contains a populated kmer_index table.
        Does not raise; returns False on any error.
        """
        if not os.path.isfile(db_path):
            return False
        conn = sqlite3.connect(db_path)
        try:
            row = conn.execute(
                "SELECT name FROM sqlite_master "
                "WHERE type='table' AND name='kmer_index'"
            ).fetchone()
            if not row:
                return False
            (count,) = conn.execute("SELECT COUNT(*) FROM kmer_index").fetchone()
            return count > 0
        except Exception:
            return False
        finally:
            conn.close()

    # ------------------------------------------------------------------
    # Search
    # ------------------------------------------------------------------

    def search(
        self,
        pattern: str,
        store: "PanIndexStore",
        engine: "PanIndexEngine",
    ) -> List[Dict]:
        """
        Find all exact occurrences of pattern across all indexed segments.

        Routing:
        - len(pattern) >= k : k-mer seed + exact verify (fast path)
        - len(pattern) <  k : linear scan fallback

        Args:
            pattern : Nucleotide sequence (case-insensitive).
            store   : PanIndexStore providing sequences + addresses.
            engine  : PanIndexEngine for deriving region addresses.

        Returns:
            List of hit dicts (same schema as SubsequenceQuery.search_pattern()).
            Sorted by (segment_id, start).
        """
        pattern = pattern.upper()
        plen = len(pattern)

        if plen < self.k:
            return self._linear_search(pattern, store, engine)

        seed_kmer = pattern[:self.k]
        candidates = self._index.get(seed_kmer, [])

        hits = []
        seen = set()

        for node_id, seed_pos in candidates:
            node = store.get_node(node_id)
            if node is None:
                continue
            seq = node['metadata'].get('seq', '')
            if not seq:
                continue

            # Exact verification of full pattern
            if seq[seed_pos:seed_pos + plen] == pattern:
                key = (node_id, seed_pos)
                if key in seen:
                    continue
                seen.add(key)

                end = seed_pos + plen
                region_context = f"{node_id}:{seed_pos}-{end}"
                parent_addr = bytes.fromhex(node['address'])
                region_address = engine.derive_ratchet_address(
                    parent_addr, region_context
                )
                hits.append({
                    'segment_id': node_id,
                    'start': seed_pos,
                    'end': end,
                    'length': plen,
                    'subsequence': pattern,
                    'region_address': region_address,
                    'parent_address': node['address'],
                    'region_context': region_context,
                    'tags': node['tags'],
                })

        hits.sort(key=lambda h: (h['segment_id'], h['start']))
        return hits

    def _linear_search(
        self,
        pattern: str,
        store: "PanIndexStore",
        engine: "PanIndexEngine",
    ) -> List[Dict]:
        """Fallback O(N * L) scan for patterns shorter than k."""
        hits = []
        for node_id in store.all_nodes():
            node = store.get_node(node_id)
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
                region_address = engine.derive_ratchet_address(
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
    # Introspection
    # ------------------------------------------------------------------

    def stats(self) -> Dict[str, int]:
        """Return build statistics."""
        total_entries = sum(len(v) for v in self._index.values())
        return {
            'k': self.k,
            'distinct_kmers': len(self._index),
            'total_entries': total_entries,
        }

    def __len__(self) -> int:
        return sum(len(v) for v in self._index.values())
