import sys
import os
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore
from kmer_index import KmerIndex


# ======================================================================
# Shared fixtures
# ======================================================================

def _build_store(k: int = 4) -> tuple:
    """Build a small engine + store with predictable sequences."""
    engine = PanIndexEngine(pangenome_seed=b"test_kmer_index_seed_000000000")
    store = PanIndexStore()
    segments = [
        ("seg1", "ATGCGTCGTAGCTAGC"),
        ("seg2", "GCTAGCTAGCATGCGT"),
        ("seg3", "TTTTAAAACCCCGGGG"),
        ("seg4", "ATGCATGCATGCATGC"),
    ]
    parent = engine.root_hash
    for node_id, seq in segments:
        addr = engine.derive_ratchet_address(parent, node_id)
        store.insert(node_id, addr, [], {'seq': seq})
        parent = addr
    return engine, store, segments


# ======================================================================
# KmerIndex.build
# ======================================================================

class TestKmerIndexBuild(unittest.TestCase):
    """Build must index all k-mers from all sequences."""

    def setUp(self):
        _, self.store, self.segments = _build_store()
        self.ki = KmerIndex.build(self.store, k=4)

    def test_stats_returns_k(self):
        st = self.ki.stats()
        self.assertEqual(st['k'], 4)

    def test_distinct_kmers_positive(self):
        st = self.ki.stats()
        self.assertGreater(st['distinct_kmers'], 0)

    def test_total_entries_positive(self):
        st = self.ki.stats()
        self.assertGreater(st['total_entries'], 0)

    def test_len_matches_total_entries(self):
        st = self.ki.stats()
        self.assertEqual(len(self.ki), st['total_entries'])

    def test_known_kmer_indexed(self):
        # "ATGC" appears in seg1, seg2, seg4
        self.assertIn('ATGC', self.ki._index)

    def test_segment_shorter_than_k_not_indexed(self):
        engine = PanIndexEngine(pangenome_seed=b"short_seg_test_seed_0000000000")
        store = PanIndexStore()
        addr = engine.derive_ratchet_address(engine.root_hash, "short")
        store.insert("short", addr, [], {'seq': 'ATG'})  # len=3 < k=4
        ki = KmerIndex.build(store, k=4)
        self.assertEqual(len(ki), 0)

    def test_empty_sequence_not_indexed(self):
        engine = PanIndexEngine(pangenome_seed=b"empty_seq_test_seed_000000000")
        store = PanIndexStore()
        addr = engine.derive_ratchet_address(engine.root_hash, "empty")
        store.insert("empty", addr, [], {'seq': ''})
        ki = KmerIndex.build(store, k=4)
        self.assertEqual(len(ki), 0)


# ======================================================================
# KmerIndex.search - fast path (len >= k)
# ======================================================================

class TestKmerIndexSearch(unittest.TestCase):
    """Search must find all occurrences and return correct coordinates."""

    def setUp(self):
        self.engine, self.store, _ = _build_store()
        self.ki = KmerIndex.build(self.store, k=4)

    def test_exact_sequence_found(self):
        hits = self.ki.search("ATGCGTCGTAGCTAGC", self.store, self.engine)
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]['segment_id'], 'seg1')
        self.assertEqual(hits[0]['start'], 0)

    def test_partial_pattern_found_in_multiple_segments(self):
        # "ATGC" appears in seg1, seg2, seg4
        hits = self.ki.search("ATGC", self.store, self.engine)
        found_segs = {h['segment_id'] for h in hits}
        self.assertIn('seg1', found_segs)
        self.assertIn('seg4', found_segs)

    def test_absent_pattern_returns_empty(self):
        # "ACGTACGT" does not appear in any of the four test sequences:
        # seg1: ATGCGTCGTAGCTAGC
        # seg2: GCTAGCTAGCATGCGT
        # seg3: TTTTAAAACCCCGGGG
        # seg4: ATGCATGCATGCATGC
        hits = self.ki.search("ACGTACGT", self.store, self.engine)
        self.assertEqual(hits, [])

    def test_case_insensitive(self):
        hits_upper = self.ki.search("ATGC", self.store, self.engine)
        hits_lower = self.ki.search("atgc", self.store, self.engine)
        self.assertEqual(len(hits_upper), len(hits_lower))

    def test_results_sorted_by_segment_then_position(self):
        hits = self.ki.search("ATGC", self.store, self.engine)
        keys = [(h['segment_id'], h['start']) for h in hits]
        self.assertEqual(keys, sorted(keys))

    def test_hit_address_is_32_bytes(self):
        hits = self.ki.search("ATGC", self.store, self.engine)
        self.assertGreater(len(hits), 0)
        self.assertEqual(len(hits[0]['region_address']), 32)

    def test_hit_coordinates_correct(self):
        hits = self.ki.search("ATGCGTCGTAGCTAGC", self.store, self.engine)
        self.assertEqual(len(hits), 1)
        h = hits[0]
        self.assertEqual(h['start'], 0)
        self.assertEqual(h['end'], 16)
        self.assertEqual(h['length'], 16)

    def test_no_duplicate_hits(self):
        hits = self.ki.search("ATGC", self.store, self.engine)
        seen = set()
        for h in hits:
            key = (h['segment_id'], h['start'])
            self.assertNotIn(key, seen, f"Duplicate hit: {key}")
            seen.add(key)


# ======================================================================
# KmerIndex.search - short pattern fallback (len < k)
# ======================================================================

class TestKmerIndexFallback(unittest.TestCase):
    """Patterns shorter than k must fall back to linear scan."""

    def setUp(self):
        self.engine, self.store, _ = _build_store()
        # Build with k=8 so "ATG" (len=3) triggers fallback
        self.ki = KmerIndex.build(self.store, k=8)

    def test_short_pattern_finds_results(self):
        # "ATG" is 3 bp, k=8 -> linear scan
        hits = self.ki.search("ATG", self.store, self.engine)
        self.assertGreater(len(hits), 0)

    def test_short_pattern_absent_returns_empty(self):
        hits = self.ki.search("ZZZ", self.store, self.engine)
        self.assertEqual(hits, [])


# ======================================================================
# KmerIndex save / load round-trip
# ======================================================================

class TestKmerIndexPersistence(unittest.TestCase):
    """Save -> load must produce identical search results."""

    def setUp(self):
        self.engine, self.store, _ = _build_store()
        self.ki_original = KmerIndex.build(self.store, k=4)
        fd, self.db_path = tempfile.mkstemp(suffix='.frx.db')
        os.close(fd)
        # Write main store first (required table structure)
        self.store.save(self.db_path)
        self.ki_original.save(self.db_path)

    def tearDown(self):
        try:
            os.unlink(self.db_path)
        except FileNotFoundError:
            pass

    def test_is_present_after_save(self):
        self.assertTrue(KmerIndex.is_present(self.db_path))

    def test_is_present_false_for_missing_file(self):
        self.assertFalse(KmerIndex.is_present("/tmp/does_not_exist.db"))

    def test_load_preserves_k(self):
        ki_loaded = KmerIndex.load(self.db_path, k=4)
        self.assertEqual(ki_loaded.k, 4)

    def test_load_same_search_results(self):
        ki_loaded = KmerIndex.load(self.db_path)
        hits_orig = self.ki_original.search("ATGC", self.store, self.engine)
        hits_load = ki_loaded.search("ATGC", self.store, self.engine)
        # Same number of hits and same coordinates
        self.assertEqual(len(hits_orig), len(hits_load))
        for ho, hl in zip(hits_orig, hits_load):
            self.assertEqual(ho['segment_id'], hl['segment_id'])
            self.assertEqual(ho['start'], hl['start'])

    def test_distinct_kmer_count_preserved(self):
        ki_loaded = KmerIndex.load(self.db_path)
        self.assertEqual(
            self.ki_original.stats()['distinct_kmers'],
            ki_loaded.stats()['distinct_kmers'],
        )

    def test_is_present_false_before_kmer_table(self):
        import sqlite3
        fd, path = tempfile.mkstemp(suffix='.db')
        os.close(fd)
        # Write a db with only the nodes table (no kmer_index)
        conn = sqlite3.connect(path)
        conn.execute("CREATE TABLE nodes (id TEXT)")
        conn.close()
        try:
            self.assertFalse(KmerIndex.is_present(path))
        finally:
            os.unlink(path)


if __name__ == '__main__':
    unittest.main(verbosity=2)
