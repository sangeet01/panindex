import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from meta_layer import PaninianRuleEngine, Rule
from default_rules import build_default_rule_engine


# ======================================================================
# build_default_rule_engine
# ======================================================================

class TestDefaultRuleEngineSetup(unittest.TestCase):
    """Rule engine must contain exactly five rules in correct precedence order."""

    def setUp(self):
        self.re = build_default_rule_engine()

    def test_has_five_rules(self):
        self.assertEqual(len(self.re.rules), 5)

    def test_rules_sorted_descending_by_precedence(self):
        precs = [r.precedence for r in self.re.rules]
        self.assertEqual(precs, sorted(precs, reverse=True))

    def test_highest_precedence_is_10(self):
        self.assertEqual(self.re.rules[0].precedence, 10)

    def test_lowest_precedence_is_1(self):
        self.assertEqual(self.re.rules[-1].precedence, 1)


# ======================================================================
# Rule resolution - individual rules
# ======================================================================

class TestRuleResolution(unittest.TestCase):
    """Each rule must fire correctly for matching node data."""

    def setUp(self):
        self.re = build_default_rule_engine()

    # --- genomic_segment (precedence 1, always fires) ---

    def test_empty_node_resolves_genomic(self):
        result = self.re.resolve({}, {})
        self.assertEqual(result, 'RES:genomic')

    def test_generic_node_resolves_genomic(self):
        result = self.re.resolve({'tags': ['core_chromosome'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:genomic')

    # --- resistance_candidate (precedence 3) ---

    def test_amr_tag_resolves_amr_candidate(self):
        result = self.re.resolve({'tags': ['AMR:blaTEM'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:amr_candidate')

    def test_bla_prefix_tag_resolves_amr_candidate(self):
        # blaTEM-1 has 'bla' prefix but no critical keyword (not KPC/NDM/carbapenem etc.)
        result = self.re.resolve({'tags': ['blaTEM-1'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:amr_candidate')

    def test_amr_lowercase_tag_resolves_candidate(self):
        result = self.re.resolve({'tags': ['amr_gene'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:amr_candidate')

    # --- mobile_element (precedence 5) ---

    def test_mobile_element_tag_resolves_mobile(self):
        result = self.re.resolve({'tags': ['mobile_element'], 'seq': 'ATGCATGC'}, {})
        self.assertEqual(result, 'RES:mobile')

    def test_plasmid_tag_resolves_mobile(self):
        result = self.re.resolve({'tags': ['plasmid_marker'], 'seq': 'ATGC'}, {})
        self.assertEqual(result, 'RES:mobile')

    def test_PLASMID_prefix_tag_resolves_mobile(self):
        result = self.re.resolve({'tags': ['PLASMID:pKP001'], 'seq': 'ATGC'}, {})
        self.assertEqual(result, 'RES:mobile')

    # --- amr_confirmed (precedence 8) ---

    def test_amr_and_long_seq_resolves_confirmed(self):
        long_seq = 'ATGC' * 30  # 120 bp >= 100
        result = self.re.resolve({'tags': ['AMR:blaTEM'], 'seq': long_seq}, {})
        self.assertEqual(result, 'RES:amr_confirmed')

    def test_amr_short_seq_resolves_only_candidate(self):
        short_seq = 'ATGCATGC'  # 8 bp < 100
        result = self.re.resolve({'tags': ['AMR:blaTEM'], 'seq': short_seq}, {})
        self.assertEqual(result, 'RES:amr_candidate')

    def test_long_seq_without_amr_resolves_genomic(self):
        long_seq = 'ATGC' * 30
        result = self.re.resolve({'tags': ['core_chromosome'], 'seq': long_seq}, {})
        self.assertEqual(result, 'RES:genomic')

    # --- critical_resistance (precedence 10) ---

    def test_carbapenem_resolves_critical(self):
        result = self.re.resolve({'tags': ['carbapenem_resistance'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:critical')

    def test_vancomycin_resolves_critical(self):
        result = self.re.resolve({'tags': ['vancomycin_resistance'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:critical')

    def test_NDM_tag_resolves_critical(self):
        result = self.re.resolve({'tags': ['NDM-1'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:critical')

    def test_KPC_tag_resolves_critical(self):
        result = self.re.resolve({'tags': ['KPC-3_carbapenemase'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:critical')


# ======================================================================
# Precedence ordering (Apavada beats Utsarga)
# ======================================================================

class TestRulePrecedence(unittest.TestCase):
    """Higher precedence (Apavada) must always win over lower (Utsarga)."""

    def setUp(self):
        self.re = build_default_rule_engine()

    def test_critical_beats_amr_confirmed(self):
        # carbapenem + AMR + long seq -> critical wins (prec 10 > 8)
        long_seq = 'ATGC' * 30
        result = self.re.resolve({'tags': ['AMR:blaTEM', 'carbapenem'], 'seq': long_seq}, {})
        self.assertEqual(result, 'RES:critical')

    def test_amr_confirmed_beats_mobile(self):
        # AMR + mobile_element + long seq -> amr_confirmed wins (prec 8 > 5)
        long_seq = 'ATGC' * 30
        result = self.re.resolve({'tags': ['AMR:blaTEM', 'mobile_element'], 'seq': long_seq}, {})
        self.assertEqual(result, 'RES:amr_confirmed')

    def test_mobile_beats_candidate(self):
        # AMR + mobile_element + short seq -> mobile wins (prec 5 > 3)
        result = self.re.resolve({'tags': ['AMR:blaTEM', 'mobile_element'], 'seq': 'ATG'}, {})
        self.assertEqual(result, 'RES:mobile')


# ======================================================================
# Integration with annotator
# ======================================================================

class TestAnnotatorRuleIntegration(unittest.TestCase):
    """Annotator must apply the rule engine and produce RES: tags."""

    def setUp(self):
        import tempfile
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _run_annotator(self, gfa_content: str) -> dict:
        """Run GFAAnnotator on in-memory GFA and return store."""
        import tempfile
        import os
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
        from annotator import GFAAnnotator

        fd_in, in_path = tempfile.mkstemp(suffix='.gfa', dir=self.tmp)
        fd_out, out_path = tempfile.mkstemp(suffix='_ann.gfa', dir=self.tmp)
        os.close(fd_in)
        os.close(fd_out)
        with open(in_path, 'w') as f:
            f.write(gfa_content)

        ann = GFAAnnotator(seed=b"test_rules_integration_seed_000")
        ann.annotate(in_path, out_path)
        return ann.store

    def test_plain_node_gets_res_genomic(self):
        gfa = "H\tVN:Z:1.0\nS\tnode1\tATGCATGC\n"
        store = self._run_annotator(gfa)
        node = store.get_node('node1')
        self.assertIsNotNone(node)
        tags = node['tags']
        self.assertIn('RES:genomic', tags)

    def test_amr_short_node_gets_res_amr_candidate(self):
        gfa = "H\tVN:Z:1.0\nS\tblatEM\tATGCATGC\tAN:Z:AMR:blaTEM\n"
        store = self._run_annotator(gfa)
        node = store.get_node('blatEM')
        self.assertIsNotNone(node)
        self.assertIn('RES:amr_candidate', node['tags'])

    def test_lookup_by_res_tag_works(self):
        gfa = "H\tVN:Z:1.0\nS\tgeneric\tGCTA\n"
        store = self._run_annotator(gfa)
        nodes = store.lookup_by_tag('RES:genomic')
        self.assertGreater(len(nodes), 0)

    def test_res_tag_in_streaming_annotator(self):
        import tempfile
        import os
        from streaming_annotator import StreamingGFAAnnotator

        gfa = "H\tVN:Z:1.0\nS\tstream_node\tATGCATGC\n"
        fd_in, in_path = tempfile.mkstemp(suffix='.gfa', dir=self.tmp)
        fd_out, out_path = tempfile.mkstemp(suffix='_stream.gfa', dir=self.tmp)
        os.close(fd_in)
        os.close(fd_out)
        with open(in_path, 'w') as f:
            f.write(gfa)

        ann = StreamingGFAAnnotator(seed=b"test_streaming_rule_seed_000000")
        ann.annotate(in_path, out_path)
        node = ann.store.get_node('stream_node')
        self.assertIsNotNone(node)
        self.assertIn('RES:genomic', node['tags'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
