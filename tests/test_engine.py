import sys
import os
import hashlib
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore
from parser import GFAParser
from meta_layer import PaninianRuleEngine, Rule


class TestCommutativity(unittest.TestCase):
    """
    Invariant: Commutative XOR neighborhood hash must be order-invariant.
    Splitting a node must not change addresses derived from the stable seed.
    """

    def setUp(self):
        self.engine = PanIndexEngine(pangenome_seed=b"test_seed_commutativity_0000000")

    def test_xor_order_invariant(self):
        ha = hashlib.sha256(b"neighbor_a").digest()
        hb = hashlib.sha256(b"neighbor_b").digest()
        hc = hashlib.sha256(b"neighbor_c").digest()

        result_abc = self.engine.compute_commutative_hash([ha, hb, hc])
        result_cba = self.engine.compute_commutative_hash([hc, hb, ha])
        result_bac = self.engine.compute_commutative_hash([hb, ha, hc])

        self.assertEqual(result_abc, result_cba)
        self.assertEqual(result_abc, result_bac)

    def test_single_neighbor(self):
        ha = hashlib.sha256(b"solo").digest()
        result = self.engine.compute_commutative_hash([ha])
        self.assertEqual(result, ha)

    def test_empty_neighborhood(self):
        result = self.engine.compute_commutative_hash([])
        self.assertEqual(result, b'\x00' * 32)


class TestFractalRatchet(unittest.TestCase):
    """
    Invariant: Ratchet derivation is deterministic and stable.
    The same path always produces the same address regardless of external state.
    """

    def setUp(self):
        self.engine = PanIndexEngine(pangenome_seed=b"test_seed_ratchet_stability_000")

    def test_same_path_same_address(self):
        addr_1 = self.engine.derive_ratchet_address(self.engine.root_hash, "Chr4")
        addr_2 = self.engine.derive_ratchet_address(self.engine.root_hash, "Chr4")
        self.assertEqual(addr_1, addr_2)

    def test_different_context_different_address(self):
        addr_chr4 = self.engine.derive_ratchet_address(self.engine.root_hash, "Chr4")
        addr_chr7 = self.engine.derive_ratchet_address(self.engine.root_hash, "Chr7")
        self.assertNotEqual(addr_chr4, addr_chr7)

    def test_hierarchical_derivation_stability(self):
        # Full derivation path: Root -> Chr4 -> BRCA1 -> VarA
        chr4  = self.engine.derive_ratchet_address(self.engine.root_hash, "Chr4")
        brca1 = self.engine.derive_ratchet_address(chr4, "BRCA1")
        vara  = self.engine.derive_ratchet_address(brca1, "VarA")

        # Derive again - must be identical
        chr4_v2  = self.engine.derive_ratchet_address(self.engine.root_hash, "Chr4")
        brca1_v2 = self.engine.derive_ratchet_address(chr4_v2, "BRCA1")
        vara_v2  = self.engine.derive_ratchet_address(brca1_v2, "VarA")

        self.assertEqual(vara, vara_v2)

    def test_address_length_is_32_bytes(self):
        addr = self.engine.derive_ratchet_address(self.engine.root_hash, "Chr1")
        self.assertEqual(len(addr), 32)


class TestCanonicalCycleHash(unittest.TestCase):
    """
    Invariant: A circular DNA sequence produces the same hash regardless
    of which base-pair is treated as the start (rotation invariance).
    """

    def setUp(self):
        self.engine = PanIndexEngine()

    def test_rotation_invariant(self):
        plasmid = "ATCGATCGATCG"
        # Generate all rotations
        rotations = [plasmid[i:] + plasmid[:i] for i in range(len(plasmid))]
        hashes = [self.engine.canonical_cycle_hash(r) for r in rotations]
        # All rotations must produce the same canonical hash
        self.assertEqual(len(set(hashes)), 1)

    def test_different_plasmids_different_hash(self):
        plasmid_a = "ATCGATCGATCG"
        plasmid_b = "GCTAGCTAGCTA"
        ha = self.engine.canonical_cycle_hash(plasmid_a)
        hb = self.engine.canonical_cycle_hash(plasmid_b)
        self.assertNotEqual(ha, hb)

    def test_canonical_hash_length(self):
        h = self.engine.canonical_cycle_hash("ATCG")
        self.assertEqual(len(h), 32)


class TestPanIndexStore(unittest.TestCase):
    """
    Invariant: Index must correctly insert, retrieve by address, and retrieve by tag.
    """

    def setUp(self):
        self.engine = PanIndexEngine(pangenome_seed=b"test_seed_index_00000000000000")
        self.store = PanIndexStore()

        # Insert three nodes
        nodes = [
            ("n1", "ACTG", ["upstream"]),
            ("n2", "A",    ["SNP", "AMR:blaTEM"]),
            ("n3", "T",    ["SNP"]),
        ]
        parent = self.engine.root_hash
        for nid, seq, tags in nodes:
            derived = self.engine.derive_ratchet_address(parent, nid)
            addr = self.engine.compute_node_address(seq, derived)
            self.store.insert(nid, addr, tags, {'seq': seq})
            self._last_addrs = getattr(self, '_last_addrs', {})
            self._last_addrs[nid] = addr
            parent = addr

    def test_lookup_by_address(self):
        found = self.store.lookup_by_address(self._last_addrs["n2"])
        self.assertEqual(found, "n2")

    def test_lookup_by_tag_single(self):
        nodes = self.store.lookup_by_tag("AMR:blaTEM")
        self.assertIn("n2", nodes)

    def test_lookup_by_tag_multi(self):
        nodes = self.store.lookup_by_tag("SNP")
        self.assertIn("n2", nodes)
        self.assertIn("n3", nodes)

    def test_lookup_nonexistent_address(self):
        fake = hashlib.sha256(b"not_in_index").digest()
        result = self.store.lookup_by_address(fake)
        self.assertIsNone(result)

    def test_lookup_nonexistent_tag(self):
        result = self.store.lookup_by_tag("NONEXISTENT")
        self.assertEqual(result, [])

    def test_stats(self):
        s = self.store.stats()
        self.assertEqual(s['total_nodes'], 3)
        self.assertGreater(s['unique_tags'], 0)


class TestPaninianRuleEngine(unittest.TestCase):
    """
    Invariant: Higher-precedence (Apavada) rules override general (Utsarga) rules.
    """

    def setUp(self):
        self.engine = PaninianRuleEngine()
        self.engine.add_rule(Rule(
            name="General",
            precedence=1,
            condition=lambda n, c: True,
            action=lambda n, c: "general"
        ))
        self.engine.add_rule(Rule(
            name="AMR_Override",
            precedence=10,
            condition=lambda n, c: "AMR" in n.get('tags', []),
            action=lambda n, c: "amr_pinpoint"
        ))

    def test_general_rule_fires(self):
        node = {'tags': []}
        self.assertEqual(self.engine.resolve(node, {}), "general")

    def test_apavada_overrides_utsarga(self):
        node = {'tags': ['AMR']}
        self.assertEqual(self.engine.resolve(node, {}), "amr_pinpoint")

    def test_no_rules_returns_default(self):
        empty_engine = PaninianRuleEngine()
        node = {'tags': ['AMR']}
        self.assertEqual(empty_engine.resolve(node, {}), "default_resolution")


class TestGFAParser(unittest.TestCase):
    """
    Invariant: Parser must assign unique non-None addresses to every node in a GFA.
    """

    def setUp(self):
        self.engine = PanIndexEngine(pangenome_seed=b"test_seed_parser_000000000000a")
        self.gfa_path = os.path.join(
            os.path.dirname(__file__), '..', 'test.gfa'
        )

    def test_all_nodes_get_addresses(self):
        parser = GFAParser(self.engine)
        parser.parse_file(self.gfa_path)
        parser.compute_all_addresses()

        for nid, data in parser.nodes.items():
            self.assertIsNotNone(data['panindex_addr'],
                                 f"Node {nid} has no address")
            self.assertEqual(len(data['panindex_addr']), 32)

    def test_addresses_are_unique(self):
        parser = GFAParser(self.engine)
        parser.parse_file(self.gfa_path)
        parser.compute_all_addresses()

        addrs = [data['panindex_addr'] for data in parser.nodes.values()]
        self.assertEqual(len(addrs), len(set(addrs)),
                         "Duplicate addresses detected")


if __name__ == "__main__":
    unittest.main(verbosity=2)
