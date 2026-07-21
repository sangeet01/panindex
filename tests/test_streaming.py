import sys
import os
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from annotator import GFAAnnotator
from streaming_annotator import StreamingGFAAnnotator, make_annotator


# ======================================================================
# Shared GFA fixture
# ======================================================================

SIMPLE_GFA = (
    "H\tVN:Z:1.0\n"
    "S\t1\tCAAATAAG\n"
    "S\t2\tA\n"
    "S\t3\tT\n"
    "L\t1\t+\t2\t+\t0M\n"
    "L\t1\t+\t3\t+\t0M\n"
)

STAR_GFA = (
    "H\tVN:Z:1.0\n"
    "S\tseg1\t*\n"
    "S\tseg2\tATCG\n"
    "L\tseg1\t+\tseg2\t+\t0M\n"
)

TAGGED_GFA = (
    "H\tVN:Z:1.0\n"
    "S\tgene1\tATGCATGC\tAN:Z:AMR:blaTEM\n"
    "S\tgene2\tGCTAGCTA\n"
    "L\tgene1\t+\tgene2\t+\t0M\n"
)


def _write_temp(content: str, suffix: str) -> str:
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, 'w', encoding='utf-8') as f:
        f.write(content)
    return path


def _read_slines(path: str) -> dict:
    """Return {node_id: (seq, extra_fields)} from S-lines in a GFA."""
    result = {}
    with open(path, encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == 'S' and len(parts) >= 3:
                result[parts[1]] = (parts[2], parts[3:])
    return result


# ======================================================================
# StreamingGFAAnnotator - correctness
# ======================================================================

class TestStreamingAnnotatorCorrectness(unittest.TestCase):
    """Streaming output must match standard in-memory annotator output."""

    def setUp(self):
        self.gfa_path = _write_temp(SIMPLE_GFA, '.gfa')

    def tearDown(self):
        os.unlink(self.gfa_path)
        for attr in ['std_out', 'stream_out']:
            p = getattr(self, attr, None)
            if p and os.path.exists(p):
                os.unlink(p)

    def _annotate_standard(self) -> str:
        fd, path = tempfile.mkstemp(suffix='_std.gfa')
        os.close(fd)
        self.std_out = path
        ann = GFAAnnotator(seed=b"test_streaming_seed_00000000000")
        ann.annotate(self.gfa_path, path)
        return path

    def _annotate_streaming(self) -> str:
        fd, path = tempfile.mkstemp(suffix='_stream.gfa')
        os.close(fd)
        self.stream_out = path
        ann = StreamingGFAAnnotator(seed=b"test_streaming_seed_00000000000")
        ann.annotate(self.gfa_path, path)
        return path

    def test_same_segment_count(self):
        std = _read_slines(self._annotate_standard())
        stream = _read_slines(self._annotate_streaming())
        self.assertEqual(set(std.keys()), set(stream.keys()))

    def test_sequences_preserved(self):
        std = _read_slines(self._annotate_standard())
        stream = _read_slines(self._annotate_streaming())
        for node_id in std:
            self.assertEqual(std[node_id][0], stream[node_id][0])

    def test_an_tag_injected(self):
        out = self._annotate_streaming()
        slines = _read_slines(out)
        for node_id, (seq, extras) in slines.items():
            an_fields = [f for f in extras if f.startswith('AN:Z:')]
            self.assertEqual(len(an_fields), 1,
                             f"Node {node_id} missing AN:Z: tag")

    def test_pa_tag_injected(self):
        out = self._annotate_streaming()
        slines = _read_slines(out)
        for node_id, (seq, extras) in slines.items():
            pa_fields = [f for f in extras if f.startswith('PA:Z:')]
            self.assertEqual(len(pa_fields), 1,
                             f"Node {node_id} missing PA:Z: tag")

    def test_af_tag_injected(self):
        out = self._annotate_streaming()
        slines = _read_slines(out)
        for node_id, (seq, extras) in slines.items():
            af_fields = [f for f in extras if f.startswith('AF:i:')]
            self.assertEqual(len(af_fields), 1,
                             f"Node {node_id} missing AF:i: tag")

    def test_link_lines_preserved(self):
        out = self._annotate_streaming()
        with open(out, encoding='utf-8') as f:
            lines = f.readlines()
        link_lines = [l for l in lines if l.startswith('L\t')]
        self.assertEqual(len(link_lines), 2)

    def test_header_line_preserved(self):
        out = self._annotate_streaming()
        with open(out, encoding='utf-8') as f:
            first = f.readline()
        self.assertTrue(first.startswith('H\t'))


# ======================================================================
# StreamingGFAAnnotator - store population
# ======================================================================

class TestStreamingAnnotatorStore(unittest.TestCase):
    """The in-memory store must be populated during streaming annotation."""

    def setUp(self):
        self.gfa_path = _write_temp(SIMPLE_GFA, '.gfa')
        fd, self.out_path = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)

    def tearDown(self):
        for p in [self.gfa_path, self.out_path]:
            try:
                os.unlink(p)
            except FileNotFoundError:
                pass

    def test_store_has_all_segments(self):
        ann = StreamingGFAAnnotator(seed=b"test_store_pop_seed_000000000000")
        ann.annotate(self.gfa_path, self.out_path)
        self.assertEqual(len(ann.store), 3)

    def test_store_node_ids_match_gfa(self):
        ann = StreamingGFAAnnotator(seed=b"test_store_pop_seed_000000000001")
        ann.annotate(self.gfa_path, self.out_path)
        self.assertIn('1', ann.store.all_nodes())
        self.assertIn('2', ann.store.all_nodes())
        self.assertIn('3', ann.store.all_nodes())

    def test_store_sequences_correct(self):
        ann = StreamingGFAAnnotator(seed=b"test_store_pop_seed_000000000002")
        ann.annotate(self.gfa_path, self.out_path)
        node = ann.store.get_node('1')
        self.assertIsNotNone(node)
        self.assertEqual(node['metadata']['seq'], 'CAAATAAG')

    def test_nodes_written_count(self):
        ann = StreamingGFAAnnotator(seed=b"test_store_pop_seed_000000000003")
        ann.annotate(self.gfa_path, self.out_path)
        self.assertEqual(ann.nodes_written, 3)

    def test_star_sequence_stored_as_empty_string(self):
        gfa_path = _write_temp(STAR_GFA, '.gfa')
        fd, out_path = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            ann = StreamingGFAAnnotator(seed=b"test_star_seq_seed_00000000000")
            ann.annotate(gfa_path, out_path)
            node = ann.store.get_node('seg1')
            self.assertIsNotNone(node)
            self.assertEqual(node['metadata']['seq'], '')
        finally:
            os.unlink(gfa_path)
            os.unlink(out_path)


# ======================================================================
# make_annotator auto-selection
# ======================================================================

class TestMakeAnnotator(unittest.TestCase):
    """make_annotator must select the right annotator based on file size."""

    def setUp(self):
        self.gfa_path = _write_temp(SIMPLE_GFA, '.gfa')

    def tearDown(self):
        os.unlink(self.gfa_path)

    def test_small_file_returns_standard(self):
        ann = make_annotator(
            self.gfa_path,
            seed=b"make_ann_test_seed_000000000000",
            mem_limit_bytes=10 * 1024 * 1024,  # 10 MB, file is tiny
        )
        from annotator import GFAAnnotator
        self.assertIsInstance(ann, GFAAnnotator)

    def test_force_streaming_overrides_size(self):
        ann = make_annotator(
            self.gfa_path,
            seed=b"make_ann_test_seed_000000000001",
            mem_limit_bytes=10 * 1024 * 1024,
            force_streaming=True,
        )
        self.assertIsInstance(ann, StreamingGFAAnnotator)

    def test_large_limit_returns_standard(self):
        ann = make_annotator(
            self.gfa_path,
            seed=b"make_ann_test_seed_000000000002",
            mem_limit_bytes=0,  # 0 bytes limit -> always stream
        )
        self.assertIsInstance(ann, StreamingGFAAnnotator)


if __name__ == '__main__':
    unittest.main(verbosity=2)
