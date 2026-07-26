import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore
from bipartite import BipartiteGraph, reverse_complement

GFA_PATH = os.path.join(os.path.dirname(__file__), '..', 'test.gfa')


def make_graph():
    engine = PanIndexEngine(pangenome_seed=b"test_bipartite_seed_0000000000aa")
    store  = PanIndexStore()
    graph  = BipartiteGraph(engine, store)
    graph.parse_gfa(GFA_PATH)
    graph.compute_addresses()
    return graph, store


class TestReverseComplement(unittest.TestCase):

    def test_complement_mapping(self):
        self.assertEqual(reverse_complement("A"), "T")
        self.assertEqual(reverse_complement("T"), "A")
        self.assertEqual(reverse_complement("G"), "C")
        self.assertEqual(reverse_complement("C"), "G")

    def test_reverse_applied(self):
        # ACTG -> complement TGAC -> reversed CAGT
        self.assertEqual(reverse_complement("ACTG"), "CAGT")

    def test_double_rc_is_identity(self):
        seq = "ATGCGTCGTA"
        self.assertEqual(reverse_complement(reverse_complement(seq)), seq)


class TestBipartiteExpansion(unittest.TestCase):

    def setUp(self):
        self.graph, self.store = make_graph()

    def test_node_count_doubled(self):
        # 4 physical nodes -> 8 states
        self.assertEqual(len(self.graph.states), 8)

    def test_all_states_addressed(self):
        for key, node in self.graph.states.items():
            self.assertIsNotNone(node.address,
                                 f"State {key} has no address")
            self.assertEqual(len(node.address), 32)

    def test_forward_and_reverse_have_distinct_addresses(self):
        for nid in ['1', '2', '3', '4']:
            fwd_addr = self.graph.states[f"{nid}+"].address
            rev_addr = self.graph.states[f"{nid}-"].address
            self.assertNotEqual(fwd_addr, rev_addr,
                                 f"Node {nid}: forward and reverse share an address")

    def test_reverse_state_is_rc_of_forward(self):
        for nid in ['1', '2', '3', '4']:
            fwd_seq = self.graph.states[f"{nid}+"].seq
            rev_seq = self.graph.states[f"{nid}-"].seq
            self.assertEqual(rev_seq, reverse_complement(fwd_seq),
                             f"Node {nid}-: sequence is not reverse complement of {nid}+")

    def test_all_states_indexed(self):
        for key in self.graph.states:
            node = self.store.get_node(key)
            self.assertIsNotNone(node, f"State {key} not found in store")

    def test_strand_tags_correct(self):
        fwd_node = self.store.get_node("1+")
        rev_node  = self.store.get_node("1-")
        self.assertIn("strand:+", fwd_node['tags'])
        self.assertIn("strand:-", rev_node['tags'])

    def test_addresses_unique_across_all_states(self):
        addrs = [n.address for n in self.graph.states.values() if n.address]
        self.assertEqual(len(addrs), len(set(addrs)),
                         "Duplicate addresses found in bipartite expansion")

    def test_identity_layers_are_persisted(self):
        node = self.store.get_node("1+")
        self.assertEqual(node['metadata']['stable_id'], "1+")
        self.assertEqual(len(node['metadata']['content_id']), 64)
        self.assertEqual(len(node['metadata']['topology_id']), 64)

    def test_content_lookup_returns_all_matching_placements(self):
        first = self.store.get_node("1+")['content_id']
        self.assertEqual(self.store.lookup_by_content_id(first), ["1+"])


class TestCircularCycleBreaking(unittest.TestCase):
    """
    Simulate a circular graph (plasmid: A->B->A) and verify
    that cycle detection and canonical anchor election work.
    """

    def _build_circular_gfa(self, tmp_path):
        content = (
            "H\tVN:Z:1.0\n"
            "S\tP1\tATCG\n"
            "S\tP2\tGCTA\n"
            "S\tP3\tTTAA\n"
            "L\tP1\t+\tP2\t+\t0M\n"
            "L\tP2\t+\tP3\t+\t0M\n"
            "L\tP3\t+\tP1\t+\t0M\n"  # cycle back
        )
        with open(tmp_path, 'w') as f:
            f.write(content)

    def test_circular_graph_all_addressed(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_circular.gfa')
        self._build_circular_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_circular_seed_000000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store)
        graph.parse_gfa(tmp)
        graph.compute_addresses()

        os.remove(tmp)

        for key, node in graph.states.items():
            self.assertIsNotNone(node.address,
                                 f"Circular graph: state {key} unresolved")

    def test_circular_addresses_unique(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_circular2.gfa')
        self._build_circular_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_circular2_seed_00000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store)
        graph.parse_gfa(tmp)
        graph.compute_addresses()

        os.remove(tmp)

        addrs = [n.address for n in graph.states.values() if n.address]
        self.assertEqual(len(addrs), len(set(addrs)))


class TestSingleStrandedVirus(unittest.TestCase):
    """Viral single-stranded RNA: only N+ states created."""

    def _build_single_strand_gfa(self, path):
        content = (
            "H\tVN:Z:1.0\n"
            "S\tV1\tATCG\n"
            "S\tV2\tGCTA\n"
            "L\tV1\t+\tV2\t+\t0M\n"
        )
        with open(path, 'w') as f:
            f.write(content)

    def test_single_stranded_creates_only_forward_states(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_virus.gfa')
        self._build_single_strand_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_virus_seed_000000000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store, single_stranded=True)
        graph.parse_gfa(tmp)
        os.remove(tmp)

        # Only + states should exist
        self.assertIn("V1+", graph.states)
        self.assertNotIn("V1-", graph.states)
        self.assertIn("V2+", graph.states)
        self.assertNotIn("V2-", graph.states)

    def test_single_stranded_all_addressed(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_virus2.gfa')
        self._build_single_strand_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_virus2_seed_00000000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store, single_stranded=True)
        graph.parse_gfa(tmp)
        graph.compute_addresses()
        os.remove(tmp)

        for key, node in graph.states.items():
            self.assertIsNotNone(node.address)

    def test_single_stranded_tag_applied(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_virus3.gfa')
        self._build_single_strand_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_virus3_seed_00000000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store, single_stranded=True)
        graph.parse_gfa(tmp)
        graph.compute_addresses()
        os.remove(tmp)

        node = store.get_node("V1+")
        self.assertIn("genome_type:single_stranded", node['tags'])


class TestEukaryoticMultiChromosome(unittest.TestCase):
    """Multiple disconnected subgraphs -> separate component namespaces."""

    def _build_two_chromosome_gfa(self, path):
        # Chr1: A -> B  (disconnected from Chr2)
        # Chr2: C -> D
        content = (
            "H\tVN:Z:1.0\n"
            "S\tA\tAAAA\n"
            "S\tB\tCCCC\n"
            "S\tC\tGGGG\n"
            "S\tD\tTTTT\n"
            "L\tA\t+\tB\t+\t0M\n"
            "L\tC\t+\tD\t+\t0M\n"
        )
        with open(path, 'w') as f:
            f.write(content)

    def test_two_components_detected(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_euk.gfa')
        self._build_two_chromosome_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_euk_seed_0000000000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store)
        graph.parse_gfa(tmp)
        components = graph._find_connected_components()
        os.remove(tmp)

        # Bipartite: 4 segments -> 8 states.
        # A+ and A- are connected via complement edge, B+ and B- likewise.
        # A group and B group are connected via the A->B link.
        # C group and D group are similarly linked.
        # Result: 2 disconnected bipartite chromosome groups.
        # Each group has 4 states (A+,A-,B+,B-) and (C+,C-,D+,D-).
        self.assertEqual(len(components), 2)

    def test_different_components_have_different_namespaces(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_euk2.gfa')
        self._build_two_chromosome_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_euk2_seed_000000000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store)
        graph.parse_gfa(tmp)
        graph.compute_addresses()
        os.remove(tmp)

        node_a = store.get_node("A+")
        node_c = store.get_node("C+")

        comp_a = node_a['metadata']['component']
        comp_c = node_c['metadata']['component']
        self.assertNotEqual(comp_a, comp_c,
                            "Nodes in different chromosomes should have different component IDs")

    def test_all_chromosomes_addressed(self):
        tmp = os.path.join(os.path.dirname(__file__), 'tmp_euk3.gfa')
        self._build_two_chromosome_gfa(tmp)

        engine = PanIndexEngine(pangenome_seed=b"test_euk3_seed_000000000000000")
        store  = PanIndexStore()
        graph  = BipartiteGraph(engine, store)
        graph.parse_gfa(tmp)
        graph.compute_addresses()
        os.remove(tmp)

        for key, node in graph.states.items():
            self.assertIsNotNone(node.address,
                                 f"Multi-chromosome: state {key} unresolved")


if __name__ == "__main__":
    unittest.main(verbosity=2)

