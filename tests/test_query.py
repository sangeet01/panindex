import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore
from query import PanIndexQuery, build_demo_state


class TestRatchetPathQuery(unittest.TestCase):
    """
    Mode 1: Ratchet path must derive a stable address and match the correct node.
    """

    def setUp(self):
        self.engine, self.store = build_demo_state()
        self.q = PanIndexQuery(self.engine, self.store)

    def test_exact_path_match(self):
        result = self.q.query_by_path("Root/Chr1/gene_blaTEM")
        self.assertFalse(result.is_empty())
        self.assertIn("gene_blaTEM", result.matched_nodes)

    def test_nested_path_match(self):
        result = self.q.query_by_path("Root/Chr1/gene_blaTEM/VarA")
        self.assertFalse(result.is_empty())
        self.assertIn("VarA", result.matched_nodes)

    def test_wrong_path_no_match(self):
        # Chr7 was never indexed
        result = self.q.query_by_path("Root/Chr7/gene_blaTEM")
        self.assertTrue(result.is_empty())

    def test_same_gene_different_chromosome_different_address(self):
        # Chr1/gene_blaTEM and Chr2/gene_blaTEM must produce different derived addresses
        r1 = self.q.query_by_path("Root/Chr1/gene_blaTEM")
        r2 = self.q.query_by_path("Root/Chr2/gene_blaTEM")
        self.assertNotEqual(r1.derived_address, r2.derived_address)

    def test_derived_address_is_32_bytes(self):
        result = self.q.query_by_path("Root/Chr1")
        self.assertEqual(len(result.derived_address), 32)

    def test_path_stability_across_calls(self):
        # Same path must always derive the same address
        r1 = self.q.query_by_path("Root/Chr1/gene_blaTEM")
        r2 = self.q.query_by_path("Root/Chr1/gene_blaTEM")
        self.assertEqual(r1.derived_address, r2.derived_address)

    def test_empty_path_returns_empty(self):
        result = self.q.query_by_path("")
        self.assertTrue(result.is_empty())


class TestTagQuery(unittest.TestCase):
    """
    Mode 2: Tag query must return all and only nodes with that annotation.
    """

    def setUp(self):
        self.engine, self.store = build_demo_state()
        self.q = PanIndexQuery(self.engine, self.store)

    def test_amr_tag_returns_multiple(self):
        result = self.q.query_by_tag("AMR:blaTEM")
        self.assertFalse(result.is_empty())
        self.assertGreaterEqual(len(result.matched_nodes), 2)

    def test_housekeeping_tag_returns_correct_node(self):
        result = self.q.query_by_tag("housekeeping")
        self.assertIn("gene_rpoB", result.matched_nodes)

    def test_nonexistent_tag_returns_empty(self):
        result = self.q.query_by_tag("NONEXISTENT:XYZ")
        self.assertTrue(result.is_empty())

    def test_variant_tag(self):
        result = self.q.query_by_tag("variant")
        self.assertIn("VarA", result.matched_nodes)

    def test_mobile_element_tag(self):
        result = self.q.query_by_tag("mobile_element")
        self.assertIn("gene_blaTEM", result.matched_nodes)


class TestLSHSimilarityQuery(unittest.TestCase):
    """
    Mode 3: MinHash similarity must find close sequences and reject distant ones.
    """

    def setUp(self):
        self.engine, self.store = build_demo_state()
        self.q = PanIndexQuery(self.engine, self.store)

    def test_exact_sequence_match(self):
        # Exact sequence of gene_blaTEM should be top match
        result = self.q.query_by_similarity("ATGCGTCGTAGCTAGC", threshold=0.5)
        self.assertIn("gene_blaTEM", result.matched_nodes)

    def test_one_snp_variant_found(self):
        # VarA differs from gene_blaTEM by exactly 1 SNP at the last base
        # Should be found at a moderate threshold
        result = self.q.query_by_similarity("ATGCGTCGTAGCTAGT", threshold=0.3)
        self.assertIn("VarA", result.matched_nodes)

    def test_unrelated_sequence_not_matched(self):
        # A completely unrelated sequence should not match at high threshold
        result = self.q.query_by_similarity("CCCCGGGGCCCCGGGG", threshold=0.9)
        self.assertNotIn("gene_blaTEM", result.matched_nodes)

    def test_threshold_filters_correctly(self):
        # At threshold 1.0 (perfect), only exact match should pass
        result_strict = self.q.query_by_similarity("ATGCGTCGTAGCTAGC", threshold=0.99)
        result_loose = self.q.query_by_similarity("ATGCGTCGTAGCTAGC", threshold=0.1)
        # Loose threshold must return at least as many as strict
        self.assertGreaterEqual(len(result_loose.matched_nodes),
                                len(result_strict.matched_nodes))

    def test_minhash_signature_length(self):
        sig = self.q._minhash_signature("ATCGATCG", k=4, num_hashes=64)
        self.assertEqual(len(sig), 64)

    def test_hamming_similarity_identical(self):
        sig = self.q._minhash_signature("ATCGATCGATCG")
        sim = self.q._hamming_similarity(sig, sig)
        self.assertAlmostEqual(sim, 1.0)

    def test_hamming_similarity_empty(self):
        sim = self.q._hamming_similarity([], [])
        self.assertEqual(sim, 0.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
