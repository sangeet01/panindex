import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore


# ======================================================================
# Shared fixture builder
# ======================================================================

def _build_multi_organism_store() -> tuple:
    """
    Build a store with nodes from two organisms.

    KP001 nodes: 'chr1_kp', 'chr2_kp', 'plasmid_kp'  tags: AMR:blaTEM
    EC042 nodes: 'chr1_ec', 'chr2_ec'                  tags: housekeeping
    """
    engine = PanIndexEngine(pangenome_seed=b"test_namespace_seed_000000000")
    store = PanIndexStore()

    kp_nodes = [
        ("chr1_kp",    "ATGCATGCATGC", ["AMR:blaTEM", "resistance"]),
        ("chr2_kp",    "GCTAGCTAGCTA", ["AMR:blaTEM"]),
        ("plasmid_kp", "TTTAAACCCGGG", ["AMR:blaTEM", "mobile_element"]),
    ]
    ec_nodes = [
        ("chr1_ec",    "AAACCCGGGTTTT", ["housekeeping"]),
        ("chr2_ec",    "CCCGGGAAATTT",  ["housekeeping"]),
    ]

    parent = engine.root_hash
    for node_id, seq, tags in kp_nodes:
        addr = engine.derive_ratchet_address(parent, node_id)
        store.insert(node_id, addr, tags, {'seq': seq}, organism="KP001")
        parent = addr

    parent = engine.root_hash
    for node_id, seq, tags in ec_nodes:
        addr = engine.derive_ratchet_address(parent, node_id)
        store.insert(node_id, addr, tags, {'seq': seq}, organism="EC042")
        parent = addr

    return engine, store


# ======================================================================
# insert with organism
# ======================================================================

class TestOrganismInsert(unittest.TestCase):
    """Organism field must be stored and retrievable from node records."""

    def setUp(self):
        _, self.store = _build_multi_organism_store()

    def test_organism_stored_in_record(self):
        node = self.store.get_node("chr1_kp")
        self.assertIsNotNone(node)
        self.assertEqual(node['organism'], "KP001")

    def test_ec_organism_stored(self):
        node = self.store.get_node("chr1_ec")
        self.assertIsNotNone(node)
        self.assertEqual(node['organism'], "EC042")

    def test_default_organism_is_empty_string(self):
        engine = PanIndexEngine(pangenome_seed=b"default_org_test_seed_00000000")
        store = PanIndexStore()
        addr = engine.derive_ratchet_address(engine.root_hash, "anon")
        store.insert("anon", addr, [], {'seq': 'ATGC'})
        node = store.get_node("anon")
        self.assertEqual(node['organism'], '')


# ======================================================================
# get_node with organism filter
# ======================================================================

class TestOrganismGetNode(unittest.TestCase):
    """get_node must respect organism filter."""

    def setUp(self):
        _, self.store = _build_multi_organism_store()

    def test_correct_organism_returns_node(self):
        node = self.store.get_node("chr1_kp", organism="KP001")
        self.assertIsNotNone(node)

    def test_wrong_organism_returns_none(self):
        node = self.store.get_node("chr1_kp", organism="EC042")
        self.assertIsNone(node)

    def test_no_organism_filter_returns_any(self):
        node = self.store.get_node("chr1_kp", organism=None)
        self.assertIsNotNone(node)


# ======================================================================
# all_nodes with organism filter
# ======================================================================

class TestOrganismAllNodes(unittest.TestCase):
    """all_nodes must return only nodes matching the organism filter."""

    def setUp(self):
        _, self.store = _build_multi_organism_store()

    def test_all_nodes_unscoped_returns_all(self):
        nodes = self.store.all_nodes()
        self.assertEqual(len(nodes), 5)

    def test_all_nodes_scoped_kp(self):
        nodes = self.store.all_nodes(organism="KP001")
        self.assertEqual(len(nodes), 3)
        for nid in nodes:
            self.assertIn('kp', nid)

    def test_all_nodes_scoped_ec(self):
        nodes = self.store.all_nodes(organism="EC042")
        self.assertEqual(len(nodes), 2)
        for nid in nodes:
            self.assertIn('ec', nid)

    def test_unknown_organism_returns_empty(self):
        nodes = self.store.all_nodes(organism="UNKNOWN_ORG")
        self.assertEqual(nodes, [])


# ======================================================================
# lookup_by_tag with organism filter
# ======================================================================

class TestOrganismTagLookup(unittest.TestCase):
    """lookup_by_tag must scope correctly when organism is specified."""

    def setUp(self):
        _, self.store = _build_multi_organism_store()

    def test_global_amr_tag_returns_all_organisms(self):
        nodes = self.store.lookup_by_tag("AMR:blaTEM")
        self.assertGreaterEqual(len(nodes), 3)

    def test_scoped_amr_returns_only_kp(self):
        nodes = self.store.lookup_by_tag("AMR:blaTEM", organism="KP001")
        self.assertEqual(len(nodes), 3)
        for nid in nodes:
            node = self.store.get_node(nid)
            self.assertEqual(node['organism'], "KP001")

    def test_scoped_amr_ec_returns_empty(self):
        # EC042 has no AMR:blaTEM tag
        nodes = self.store.lookup_by_tag("AMR:blaTEM", organism="EC042")
        self.assertEqual(nodes, [])

    def test_housekeeping_global(self):
        nodes = self.store.lookup_by_tag("housekeeping")
        self.assertGreaterEqual(len(nodes), 2)

    def test_housekeeping_scoped_ec(self):
        nodes = self.store.lookup_by_tag("housekeeping", organism="EC042")
        self.assertEqual(len(nodes), 2)

    def test_housekeeping_scoped_kp_returns_empty(self):
        nodes = self.store.lookup_by_tag("housekeeping", organism="KP001")
        self.assertEqual(nodes, [])

    def test_nonexistent_tag_scoped_returns_empty(self):
        nodes = self.store.lookup_by_tag("NONEXISTENT", organism="KP001")
        self.assertEqual(nodes, [])

    def test_organism_none_same_as_default(self):
        # organism=None (default) returns all
        nodes_none = self.store.lookup_by_tag("AMR:blaTEM", organism=None)
        nodes_default = self.store.lookup_by_tag("AMR:blaTEM")
        self.assertEqual(set(nodes_none), set(nodes_default))


# ======================================================================
# Cross-organism AMR scan (the core use case)
# ======================================================================

class TestCrossOrganismAMRScan(unittest.TestCase):
    """
    The global AMR scan: find all organisms that carry a resistance gene.
    FRX's unique value proposition: O(1) hash lookup vs alignment-based scan.
    """

    def setUp(self):
        _, self.store = _build_multi_organism_store()

    def test_global_scan_finds_kp(self):
        carriers = self.store.lookup_by_tag("AMR:blaTEM")
        organisms = {
            self.store.get_node(nid)['organism']
            for nid in carriers
        }
        self.assertIn("KP001", organisms)

    def test_global_scan_does_not_include_ec(self):
        carriers = self.store.lookup_by_tag("AMR:blaTEM")
        organisms = {
            self.store.get_node(nid)['organism']
            for nid in carriers
        }
        self.assertNotIn("EC042", organisms)

    def test_mobile_element_scan(self):
        carriers = self.store.lookup_by_tag("mobile_element")
        self.assertGreater(len(carriers), 0)
        for nid in carriers:
            node = self.store.get_node(nid)
            self.assertIn('AMR:blaTEM', node['tags'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
