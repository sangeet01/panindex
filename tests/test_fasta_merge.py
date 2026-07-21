import sys
import os
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore
from fasta_merge import FastaParser, FastaMerger, SubsequenceQuery


# ======================================================================
# Shared fixtures
# ======================================================================

FASTA_MULTILINE = """>seq1 description ignored
ATGC
ATGC
>seq2
GCTA
>seq3
TTAA
"""

# GFA with one * placeholder, one existing sequence, one missing from FASTA
GFA_WITH_STAR = (
    "H\tVN:Z:1.0\n"
    "S\tseq1\t*\n"           # * -> should be filled from FASTA
    "S\tseq2\tEXISTING\n"    # already has sequence -> keep
    "S\tmissing\t*\n"        # not in FASTA -> keep *
    "L\tseq1\t+\tseq2\t+\t0M\n"
)

GFA_SIMPLE = (
    "H\tVN:Z:1.0\n"
    "S\t1\tCAAATAAG\n"
    "S\t2\tA\n"
    "S\t3\tT\n"
    "L\t1\t+\t2\t+\t0M\n"
    "L\t1\t+\t3\t+\t0M\n"
)


def _write_temp(content: str, suffix: str) -> str:
    """Write content to a temp file and return its path."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, 'w', encoding='utf-8') as f:
        f.write(content)
    return path


def _build_store_from_gfa(gfa_content: str) -> tuple:
    """Build a minimal engine + store from a GFA string for query tests."""
    engine = PanIndexEngine(pangenome_seed=b"test_fasta_merge_seed_000000000")
    store = PanIndexStore()

    for line in gfa_content.splitlines():
        parts = line.split('\t')
        if parts[0] == 'S' and len(parts) >= 3 and parts[2] != '*':
            node_id = parts[1]
            seq = parts[2]
            addr = engine.derive_ratchet_address(engine.root_hash, node_id)
            store.insert(node_id, addr, [], {'seq': seq})

    return engine, store


# ======================================================================
# FastaParser
# ======================================================================

class TestFastaParser(unittest.TestCase):
    """FastaParser must correctly parse multi-line FASTA into name->seq dict."""

    def setUp(self):
        self.fasta_path = _write_temp(FASTA_MULTILINE, '.fasta')

    def tearDown(self):
        os.unlink(self.fasta_path)

    def test_record_count(self):
        seqs = FastaParser.parse(self.fasta_path)
        self.assertEqual(len(seqs), 3)

    def test_first_word_is_name(self):
        seqs = FastaParser.parse(self.fasta_path)
        self.assertIn('seq1', seqs)
        self.assertNotIn('seq1 description ignored', seqs)

    def test_multiline_concatenated(self):
        seqs = FastaParser.parse(self.fasta_path)
        # seq1 spans two lines: ATGC + ATGC = ATGCATGC
        self.assertEqual(seqs['seq1'], 'ATGCATGC')

    def test_uppercase_normalization(self):
        fasta = ">lower\natgcatgc\n"
        path = _write_temp(fasta, '.fasta')
        try:
            seqs = FastaParser.parse(path)
            self.assertEqual(seqs['lower'], 'ATGCATGC')
        finally:
            os.unlink(path)

    def test_empty_fasta_returns_empty_dict(self):
        path = _write_temp('', '.fasta')
        try:
            seqs = FastaParser.parse(path)
            self.assertEqual(seqs, {})
        finally:
            os.unlink(path)

    def test_stats_total_bp(self):
        seqs = FastaParser.parse(self.fasta_path)
        st = FastaParser.stats(seqs)
        # seq1=8, seq2=4, seq3=4 = 16
        self.assertEqual(st['total_bp'], 16)
        self.assertEqual(st['records'], 3)

    def test_stats_min_max(self):
        seqs = FastaParser.parse(self.fasta_path)
        st = FastaParser.stats(seqs)
        self.assertEqual(st['min_len'], 4)
        self.assertEqual(st['max_len'], 8)


# ======================================================================
# FastaMerger.merge
# ======================================================================

class TestFastaMergerMerge(unittest.TestCase):
    """Merger must fill * from FASTA, keep existing sequences, warn on missing."""

    def setUp(self):
        self.fasta_path = _write_temp(FASTA_MULTILINE, '.fasta')
        self.gfa_path   = _write_temp(GFA_WITH_STAR,   '.gfa')
        fd, self.out_path = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)

    def tearDown(self):
        for p in [self.fasta_path, self.gfa_path, self.out_path]:
            try:
                os.unlink(p)
            except FileNotFoundError:
                pass

    def _parse_slines(self, path):
        """Return dict of node_id -> sequence from S-lines in a GFA file."""
        result = {}
        with open(path, encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts[0] == 'S' and len(parts) >= 3:
                    result[parts[1]] = parts[2]
        return result

    def test_star_filled_from_fasta(self):
        merger = FastaMerger()
        merger.merge(self.fasta_path, self.gfa_path, self.out_path)
        slines = self._parse_slines(self.out_path)
        self.assertEqual(slines['seq1'], 'ATGCATGC')

    def test_existing_sequence_kept(self):
        merger = FastaMerger()
        merger.merge(self.fasta_path, self.gfa_path, self.out_path)
        slines = self._parse_slines(self.out_path)
        self.assertEqual(slines['seq2'], 'EXISTING')

    def test_missing_segment_stays_star(self):
        merger = FastaMerger()
        merger.merge(self.fasta_path, self.gfa_path, self.out_path)
        slines = self._parse_slines(self.out_path)
        self.assertEqual(slines['missing'], '*')

    def test_stats_filled_count(self):
        merger = FastaMerger()
        stats = merger.merge(self.fasta_path, self.gfa_path, self.out_path)
        self.assertEqual(stats['filled'], 1)
        self.assertEqual(stats['kept'], 1)
        self.assertEqual(stats['missing'], 1)

    def test_s_prefix_matching(self):
        fasta = ">42\nGCGCGCGC\n"
        gfa   = "H\tVN:Z:1.0\nS\ts42\t*\n"
        fasta_path = _write_temp(fasta, '.fasta')
        gfa_path   = _write_temp(gfa,   '.gfa')
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            merger = FastaMerger()
            stats = merger.merge(fasta_path, gfa_path, out)
            self.assertEqual(stats['filled'], 1)
        finally:
            for p in [fasta_path, gfa_path, out]:
                os.unlink(p)

    def test_link_lines_passthrough(self):
        merger = FastaMerger()
        merger.merge(self.fasta_path, self.gfa_path, self.out_path)
        with open(self.out_path, encoding='utf-8') as f:
            lines = f.readlines()
        link_lines = [l for l in lines if l.startswith('L\t')]
        self.assertEqual(len(link_lines), 1)


# ======================================================================
# FastaMerger.fasta_as_gfa
# ======================================================================

class TestFastaMergerFastaAsGfa(unittest.TestCase):
    """fasta_as_gfa must produce a valid GFA with one S-line per FASTA record."""

    def setUp(self):
        self.fasta_path = _write_temp(FASTA_MULTILINE, '.fasta')
        fd, self.out_path = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)

    def tearDown(self):
        for p in [self.fasta_path, self.out_path]:
            try:
                os.unlink(p)
            except FileNotFoundError:
                pass

    def _read_gfa(self, path):
        header = []
        segments = {}
        with open(path, encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                if parts[0] == 'H':
                    header.append(line)
                elif parts[0] == 'S' and len(parts) >= 3:
                    segments[parts[1]] = parts[2]
        return header, segments

    def test_header_line_present(self):
        merger = FastaMerger()
        merger.fasta_as_gfa(self.fasta_path, self.out_path)
        header, _ = self._read_gfa(self.out_path)
        self.assertGreater(len(header), 0)
        self.assertIn('VN:Z:1.0', header[0])

    def test_segment_count_matches_fasta(self):
        merger = FastaMerger()
        merger.fasta_as_gfa(self.fasta_path, self.out_path)
        _, segments = self._read_gfa(self.out_path)
        self.assertEqual(len(segments), 3)

    def test_sequences_correct(self):
        merger = FastaMerger()
        merger.fasta_as_gfa(self.fasta_path, self.out_path)
        _, segments = self._read_gfa(self.out_path)
        self.assertEqual(segments['seq1'], 'ATGCATGC')
        self.assertEqual(segments['seq2'], 'GCTA')

    def test_no_link_lines(self):
        merger = FastaMerger()
        merger.fasta_as_gfa(self.fasta_path, self.out_path)
        with open(self.out_path, encoding='utf-8') as f:
            lines = f.readlines()
        link_lines = [l for l in lines if l.startswith('L\t')]
        self.assertEqual(len(link_lines), 0)


# ======================================================================
# SubsequenceQuery - region extraction
# ======================================================================

class TestSubsequenceQueryRegion(unittest.TestCase):
    """Region extraction must return correct subsequences and clamp safely."""

    def setUp(self):
        self.engine, self.store = _build_store_from_gfa(GFA_SIMPLE)
        self.sq = SubsequenceQuery(self.engine, self.store)

    def test_valid_region_returns_correct_subseq(self):
        result = self.sq.extract_region('1', 0, 4)
        self.assertIsNotNone(result)
        self.assertEqual(result['subsequence'], 'CAAA')

    def test_full_region_returns_full_sequence(self):
        result = self.sq.extract_region('1', 0, 8)
        self.assertIsNotNone(result)
        self.assertEqual(result['subsequence'], 'CAAATAAG')

    def test_out_of_bounds_end_is_clamped(self):
        result = self.sq.extract_region('1', 0, 999)
        self.assertIsNotNone(result)
        self.assertEqual(result['subsequence'], 'CAAATAAG')

    def test_empty_range_returns_none(self):
        result = self.sq.extract_region('1', 5, 5)
        self.assertIsNone(result)

    def test_inverted_range_returns_none(self):
        result = self.sq.extract_region('1', 6, 3)
        self.assertIsNone(result)

    def test_unknown_segment_returns_none(self):
        result = self.sq.extract_region('nonexistent', 0, 10)
        self.assertIsNone(result)

    def test_result_length_field(self):
        result = self.sq.extract_region('1', 2, 6)
        self.assertEqual(result['length'], 4)
        self.assertEqual(len(result['subsequence']), 4)

    def test_region_address_is_32_bytes(self):
        result = self.sq.extract_region('1', 0, 4)
        self.assertIsNotNone(result)
        self.assertEqual(len(result['region_address']), 32)


# ======================================================================
# SubsequenceQuery - region string parsing
# ======================================================================

class TestSubsequenceQueryParse(unittest.TestCase):
    """Region string parsing must handle valid and invalid formats."""

    def test_valid_parse(self):
        seg, start, end = SubsequenceQuery._parse_region_string('seg1:100-200')
        self.assertEqual(seg, 'seg1')
        self.assertEqual(start, 100)
        self.assertEqual(end, 200)

    def test_segment_with_colon_in_name(self):
        # rsplit ensures only the last colon splits
        seg, start, end = SubsequenceQuery._parse_region_string('chr1:sub:0-50')
        self.assertEqual(seg, 'chr1:sub')
        self.assertEqual(start, 0)
        self.assertEqual(end, 50)

    def test_missing_colon_raises(self):
        with self.assertRaises(ValueError):
            SubsequenceQuery._parse_region_string('seg1_100-200')

    def test_missing_dash_raises(self):
        with self.assertRaises(ValueError):
            SubsequenceQuery._parse_region_string('seg1:100')

    def test_non_integer_raises(self):
        with self.assertRaises(ValueError):
            SubsequenceQuery._parse_region_string('seg1:abc-200')

    def test_query_method_delegates_correctly(self):
        engine, store = _build_store_from_gfa(GFA_SIMPLE)
        sq = SubsequenceQuery(engine, store)
        result = sq.query('1:0-4')
        self.assertIsNotNone(result)
        self.assertEqual(result['subsequence'], 'CAAA')


# ======================================================================
# SubsequenceQuery - pattern search
# ======================================================================

class TestSubsequenceQueryPattern(unittest.TestCase):
    """Pattern search must find all occurrences and return correct coordinates."""

    def setUp(self):
        self.engine, self.store = _build_store_from_gfa(GFA_SIMPLE)
        self.sq = SubsequenceQuery(self.engine, self.store)

    def test_exact_pattern_found(self):
        hits = self.sq.search_pattern('CAAATAAG')
        self.assertEqual(len(hits), 1)
        self.assertEqual(hits[0]['segment_id'], '1')
        self.assertEqual(hits[0]['start'], 0)
        self.assertEqual(hits[0]['end'], 8)

    def test_partial_pattern_found(self):
        hits = self.sq.search_pattern('CAAA')
        self.assertGreater(len(hits), 0)
        self.assertEqual(hits[0]['start'], 0)

    def test_absent_pattern_returns_empty(self):
        hits = self.sq.search_pattern('ZZZZ')
        self.assertEqual(hits, [])

    def test_case_insensitive_input(self):
        hits_upper = self.sq.search_pattern('CAAATAAG')
        hits_lower = self.sq.search_pattern('caaataag')
        self.assertEqual(len(hits_upper), len(hits_lower))

    def test_hit_address_is_32_bytes(self):
        hits = self.sq.search_pattern('CAAA')
        self.assertGreater(len(hits), 0)
        self.assertEqual(len(hits[0]['region_address']), 32)

    def test_two_different_patterns_different_addresses(self):
        hits_a = self.sq.search_pattern('CAAA')
        hits_b = self.sq.search_pattern('TAAG')
        if hits_a and hits_b and hits_a[0]['segment_id'] == hits_b[0]['segment_id']:
            self.assertNotEqual(
                hits_a[0]['region_address'],
                hits_b[0]['region_address']
            )

    def test_results_sorted_by_segment_then_position(self):
        hits = self.sq.search_pattern('A')
        keys = [(h['segment_id'], h['start']) for h in hits]
        self.assertEqual(keys, sorted(keys))


# ======================================================================
# SubsequenceQuery - address stability
# ======================================================================

class TestSubsequenceQueryAddressStability(unittest.TestCase):
    """Region addresses must be deterministic and stable across calls."""

    def setUp(self):
        self.engine, self.store = _build_store_from_gfa(GFA_SIMPLE)
        self.sq = SubsequenceQuery(self.engine, self.store)

    def test_same_region_same_address(self):
        r1 = self.sq.extract_region('1', 0, 4)
        r2 = self.sq.extract_region('1', 0, 4)
        self.assertEqual(r1['region_address'], r2['region_address'])

    def test_different_regions_different_addresses(self):
        r1 = self.sq.extract_region('1', 0, 4)
        r2 = self.sq.extract_region('1', 4, 8)
        self.assertNotEqual(r1['region_address'], r2['region_address'])

    def test_region_address_differs_from_parent(self):
        result = self.sq.extract_region('1', 0, 4)
        self.assertIsNotNone(result)
        parent_hex = result['parent_address']
        region_hex = result['region_address'].hex()
        self.assertNotEqual(parent_hex, region_hex)


if __name__ == '__main__':
    unittest.main(verbosity=2)
