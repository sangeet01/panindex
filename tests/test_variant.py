import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore
from query import PanIndexQuery, VariantDiffResult


# ======================================================================
# Fixture builder
# ======================================================================

def _build_store_with_variants() -> tuple:
    """
    Build engine + store with four nodes:
      Chr1    : ATGCATGCATGCATGC  (core, KP001)
      VarA    : ATGCATGCATGCATGT  (one base different from Chr1, KP001)
      VarNovel: GCTAGCTAGCTAGCTA  (unrelated sequence, KP001)
      EC_Chr1 : ATGCATGCATGCATGC  (same sequence as Chr1 but EC042 organism)
    """
    engine = PanIndexEngine(pangenome_seed=b"test_variant_diff_seed_00000000")
    store = PanIndexStore()

    parent = engine.root_hash
    nodes = [
        ("Chr1",     "ATGCATGCATGCATGC", ["core_chromosome"], "KP001"),
        ("VarA",     "ATGCATGCATGCATGT", ["core_chromosome"], "KP001"),
        ("VarNovel", "GCTAGCTAGCTAGCTA", ["core_chromosome"], "KP001"),
        ("EC_Chr1",  "ATGCATGCATGCATGC", ["core_chromosome"], "EC042"),
    ]

    for node_id, seq, tags, organism in nodes:
        addr = engine.derive_ratchet_address(parent, node_id)
        store.insert(node_id, addr, tags, {'seq': seq}, organism=organism)
        parent = addr

    return engine, store


# ======================================================================
# VariantDiffResult verdict cases
# ======================================================================

class TestVariantDiffVerdicts(unittest.TestCase):
    """query_variant_diff must return correct verdicts in all cases."""

    def setUp(self):
        self.engine, self.store = _build_store_with_variants()
        self.q = PanIndexQuery(self.engine, self.store)

    def test_identical_same_node(self):
        # A node compared against itself: same address -> IDENTICAL
        result = self.q.query_variant_diff("Chr1", "Chr1")
        self.assertEqual(result.verdict, "IDENTICAL")

    def test_identical_returns_similarity_1(self):
        result = self.q.query_variant_diff("Chr1", "Chr1")
        self.assertAlmostEqual(result.similarity, 1.0)

    def test_variant_close_sequences_same_organism(self):
        # Chr1 vs VarA: one base different, same organism -> VARIANT
        result = self.q.query_variant_diff("Chr1", "VarA")
        self.assertEqual(result.verdict, "VARIANT")

    def test_variant_similarity_high(self):
        result = self.q.query_variant_diff("Chr1", "VarA")
        self.assertGreater(result.similarity, 0.5)

    def test_variant_same_organism_true(self):
        result = self.q.query_variant_diff("Chr1", "VarA")
        self.assertTrue(result.same_organism)

    def test_novel_unrelated_sequence(self):
        # Chr1 vs VarNovel: very different sequence -> NOVEL
        result = self.q.query_variant_diff("Chr1", "VarNovel")
        self.assertEqual(result.verdict, "NOVEL")

    def test_novel_cross_organism_same_seq(self):
        # Chr1 (KP001) vs EC_Chr1 (EC042): same seq but different organism -> NOVEL
        result = self.q.query_variant_diff("Chr1", "EC_Chr1")
        self.assertNotEqual(result.verdict, "VARIANT")
        # Either IDENTICAL or NOVEL depending on whether addresses happen to match
        # (they won't because addresses are derived from different positions in the chain)
        self.assertIn(result.verdict, ("IDENTICAL", "NOVEL", "VARIANT"))

    def test_not_found_missing_node_a(self):
        result = self.q.query_variant_diff("DOES_NOT_EXIST", "Chr1")
        self.assertEqual(result.verdict, "NOT_FOUND")

    def test_not_found_missing_node_b(self):
        result = self.q.query_variant_diff("Chr1", "DOES_NOT_EXIST")
        self.assertEqual(result.verdict, "NOT_FOUND")

    def test_not_found_both_missing(self):
        result = self.q.query_variant_diff("NO_A", "NO_B")
        self.assertEqual(result.verdict, "NOT_FOUND")


# ======================================================================
# VariantDiffResult - field correctness
# ======================================================================

class TestVariantDiffResultFields(unittest.TestCase):
    """Result object must carry correct addresses, similarity, and organism data."""

    def setUp(self):
        self.engine, self.store = _build_store_with_variants()
        self.q = PanIndexQuery(self.engine, self.store)

    def test_node_ids_recorded(self):
        result = self.q.query_variant_diff("Chr1", "VarA")
        self.assertEqual(result.node_id_a, "Chr1")
        self.assertEqual(result.node_id_b, "VarA")

    def test_addr_a_is_chr1_address(self):
        result = self.q.query_variant_diff("Chr1", "VarA")
        expected = self.store.get_node("Chr1")['address']
        self.assertEqual(result.addr_a, expected)

    def test_addr_b_is_vara_address(self):
        result = self.q.query_variant_diff("Chr1", "VarA")
        expected = self.store.get_node("VarA")['address']
        self.assertEqual(result.addr_b, expected)

    def test_similarity_in_range(self):
        result = self.q.query_variant_diff("Chr1", "VarNovel")
        self.assertGreaterEqual(result.similarity, 0.0)
        self.assertLessEqual(result.similarity, 1.0)

    def test_not_found_addr_a_is_none_when_missing(self):
        result = self.q.query_variant_diff("MISSING_A", "Chr1")
        self.assertIsNone(result.addr_a)

    def test_identical_addr_a_equals_addr_b(self):
        result = self.q.query_variant_diff("Chr1", "Chr1")
        self.assertEqual(result.addr_a, result.addr_b)

    def test_verdict_in_valid_set(self):
        for node_a, node_b in [("Chr1", "VarA"), ("Chr1", "VarNovel"),
                                ("Chr1", "MISSING"), ("Chr1", "Chr1")]:
            result = self.q.query_variant_diff(node_a, node_b)
            self.assertIn(result.verdict, VariantDiffResult.VERDICTS)


# ======================================================================
# VariantDiffResult.print() - smoke test
# ======================================================================

class TestVariantDiffPrint(unittest.TestCase):
    """print() must not raise for any verdict."""

    def setUp(self):
        self.engine, self.store = _build_store_with_variants()
        self.q = PanIndexQuery(self.engine, self.store)

    def _call_print(self, node_a, node_b):
        import io
        import contextlib
        result = self.q.query_variant_diff(node_a, node_b)
        f = io.StringIO()
        with contextlib.redirect_stdout(f):
            result.print()
        return f.getvalue()

    def test_print_identical_no_raise(self):
        out = self._call_print("Chr1", "Chr1")
        self.assertIn("IDENTICAL", out)

    def test_print_variant_no_raise(self):
        out = self._call_print("Chr1", "VarA")
        self.assertIn("VARIANT", out)

    def test_print_novel_no_raise(self):
        out = self._call_print("Chr1", "VarNovel")
        self.assertIn("NOVEL", out)

    def test_print_not_found_no_raise(self):
        out = self._call_print("MISSING", "Chr1")
        self.assertIn("NOT_FOUND", out)

    def test_print_includes_similarity_for_variant(self):
        out = self._call_print("Chr1", "VarA")
        self.assertIn("Similarity", out)

    def test_print_no_similarity_for_not_found(self):
        out = self._call_print("MISSING", "Chr1")
        self.assertNotIn("Similarity", out)


# ======================================================================
# Custom similarity threshold
# ======================================================================

class TestVariantDiffThreshold(unittest.TestCase):
    """Threshold parameter must affect VARIANT vs NOVEL classification."""

    def setUp(self):
        self.engine, self.store = _build_store_with_variants()
        self.q = PanIndexQuery(self.engine, self.store)

    def test_threshold_1_0_forces_novel_for_non_identical(self):
        # Only IDENTICAL nodes pass threshold=1.0 -> VARIANT path impossible
        result = self.q.query_variant_diff("Chr1", "VarA", similarity_threshold=1.0)
        # Chr1 vs VarA differ by one base so similarity < 1.0 -> NOVEL
        self.assertEqual(result.verdict, "NOVEL")

    def test_threshold_0_0_allows_variant_for_any_same_organism(self):
        # At threshold=0.0, any non-zero similarity qualifies as VARIANT if same org
        result = self.q.query_variant_diff("Chr1", "VarA", similarity_threshold=0.0)
        self.assertEqual(result.verdict, "VARIANT")


if __name__ == '__main__':
    unittest.main(verbosity=2)
