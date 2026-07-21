import sys
import os
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from vg_frx import VGFRXNormalizer, VGNormalizationError, run_vg_import


# ======================================================================
# Shared GFA fixtures
# ======================================================================

# Standard GFA 1.0 output from vg (with vg-internal tags)
VG_GFA_10 = (
    "H\tVN:Z:1.0\n"
    "S\t1\tACGTACGTACGTACGT\tLN:i:16\tRC:i:12\n"
    "S\t2\tTGCATGCATGCATGCA\tLN:i:16\tRC:i:8\tFC:i:3\n"
    "S\t3\tGGCCGGCCGGCCGGCC\tLN:i:16\n"
    "L\t1\t+\t2\t+\t0M\n"
    "L\t2\t+\t3\t+\t0M\n"
    "P\tpath1\t1+,2+,3+\t*\n"
)

# GFA 1.1 output from newer vg (has W-lines)
VG_GFA_11 = (
    "H\tVN:Z:1.1\n"
    "S\t1\tACGT\tLN:i:4\n"
    "S\t2\tTGCA\tLN:i:4\n"
    "L\t1\t+\t2\t+\t0M\n"
    "W\tNA12878\t0\tchr1\t0\t4\t>1>2\n"
)

# Minimal valid GFA (no vg tags, just pure S/L)
PURE_GFA = (
    "H\tVN:Z:1.0\n"
    "S\tseg1\tATGCATGC\n"
    "S\tseg2\tGCTAGCTA\n"
    "L\tseg1\t+\tseg2\t+\t0M\n"
)

# GFA with extra S-line tags that are not vg-internal (should be preserved)
GFA_CUSTOM_TAGS = (
    "H\tVN:Z:1.0\n"
    "S\tseg1\tATGC\tAN:Z:my_annotation\tLN:i:4\n"
)

NOT_GFA = "This is not a GFA file at all.\n"


def _write_temp(content: str) -> str:
    fd, path = tempfile.mkstemp(suffix='.gfa')
    with os.fdopen(fd, 'w', encoding='utf-8') as f:
        f.write(content)
    return path


# ======================================================================
# VGFRXNormalizer - basic operation
# ======================================================================

class TestNormalizerBasic(unittest.TestCase):
    """Normalizer must produce valid GFA output and count lines correctly."""

    def setUp(self):
        self.norm = VGFRXNormalizer()
        self.gfa_path = _write_temp(VG_GFA_10)

    def tearDown(self):
        os.unlink(self.gfa_path)

    def test_returns_stats_dict(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            stats = self.norm.normalize(self.gfa_path, out)
            self.assertIsInstance(stats, dict)
        finally:
            os.unlink(out)

    def test_s_line_count(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            stats = self.norm.normalize(self.gfa_path, out)
            self.assertEqual(stats['s_lines'], 3)
        finally:
            os.unlink(out)

    def test_l_line_count(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            stats = self.norm.normalize(self.gfa_path, out)
            self.assertEqual(stats['l_lines'], 2)
        finally:
            os.unlink(out)

    def test_p_line_count(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            stats = self.norm.normalize(self.gfa_path, out)
            self.assertEqual(stats['p_lines'], 1)
        finally:
            os.unlink(out)

    def test_tags_stripped_count(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            stats = self.norm.normalize(self.gfa_path, out)
            # seg1 has 2 vg tags (LN:i:, RC:i:), seg2 has 3, seg3 has 1 -> 6 total
            self.assertEqual(stats['tags_stripped'], 6)
        finally:
            os.unlink(out)

    def test_output_has_no_ln_tag(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            self.norm.normalize(self.gfa_path, out)
            with open(out) as f:
                content = f.read()
            self.assertNotIn('LN:i:', content)
        finally:
            os.unlink(out)

    def test_output_has_no_rc_tag(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            self.norm.normalize(self.gfa_path, out)
            with open(out) as f:
                content = f.read()
            self.assertNotIn('RC:i:', content)
        finally:
            os.unlink(out)

    def test_output_sequences_preserved(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            self.norm.normalize(self.gfa_path, out)
            with open(out) as f:
                content = f.read()
            self.assertIn('ACGT', content)
            self.assertIn('TGCA', content)
            self.assertIn('GGCC', content)
        finally:
            os.unlink(out)


# ======================================================================
# VGFRXNormalizer - GFA 1.1 / W-line pass-through
# ======================================================================

class TestNormalizerGFA11(unittest.TestCase):
    """W-lines from GFA 1.1 must pass through unchanged."""

    def setUp(self):
        self.norm = VGFRXNormalizer()
        self.gfa_path = _write_temp(VG_GFA_11)

    def tearDown(self):
        os.unlink(self.gfa_path)

    def test_w_line_count_in_stats(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            stats = self.norm.normalize(self.gfa_path, out)
            self.assertEqual(stats['w_lines'], 1)
        finally:
            os.unlink(out)

    def test_w_line_preserved_in_output(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            self.norm.normalize(self.gfa_path, out)
            with open(out) as f:
                content = f.read()
            self.assertIn('W\tNA12878', content)
        finally:
            os.unlink(out)

    def test_gfa_version_detected_11(self):
        self.assertEqual(self.norm.detect_gfa_version(self.gfa_path), "1.1")


# ======================================================================
# VGFRXNormalizer - keep_vg_tags flag
# ======================================================================

class TestNormalizerKeepTags(unittest.TestCase):
    """keep_vg_tags=True must preserve all vg-internal tags."""

    def setUp(self):
        self.norm = VGFRXNormalizer(keep_vg_tags=True)
        self.gfa_path = _write_temp(VG_GFA_10)

    def tearDown(self):
        os.unlink(self.gfa_path)

    def test_ln_tag_preserved(self):
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            stats = self.norm.normalize(self.gfa_path, out)
            self.assertEqual(stats['tags_stripped'], 0)
            with open(out) as f:
                content = f.read()
            self.assertIn('LN:i:', content)
        finally:
            os.unlink(out)


# ======================================================================
# VGFRXNormalizer - custom non-vg tags preserved
# ======================================================================

class TestNormalizerCustomTags(unittest.TestCase):
    """Non-vg custom tags on S-lines must NOT be stripped."""

    def test_custom_an_tag_preserved(self):
        norm = VGFRXNormalizer()
        gfa_path = _write_temp(GFA_CUSTOM_TAGS)
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            norm.normalize(gfa_path, out)
            with open(out) as f:
                content = f.read()
            self.assertIn('AN:Z:my_annotation', content)
            self.assertNotIn('LN:i:', content)
        finally:
            os.unlink(gfa_path)
            os.unlink(out)


# ======================================================================
# VGFRXNormalizer - error cases
# ======================================================================

class TestNormalizerErrors(unittest.TestCase):
    """Normalizer must raise appropriate exceptions on bad input."""

    def setUp(self):
        self.norm = VGFRXNormalizer()

    def test_file_not_found_raises(self):
        with self.assertRaises(FileNotFoundError):
            self.norm.normalize('/tmp/does_not_exist_frx_test.gfa', '/tmp/out.gfa')

    def test_not_gfa_file_raises(self):
        bad_path = _write_temp(NOT_GFA)
        fd, out = tempfile.mkstemp(suffix='.gfa')
        os.close(fd)
        try:
            with self.assertRaises(VGNormalizationError):
                self.norm.normalize(bad_path, out)
        finally:
            os.unlink(bad_path)
            try:
                os.unlink(out)
            except FileNotFoundError:
                pass


# ======================================================================
# VGFRXNormalizer - normalize_to_tempfile
# ======================================================================

class TestNormalizerTempfile(unittest.TestCase):
    """normalize_to_tempfile must return a path that exists."""

    def test_tempfile_exists(self):
        norm = VGFRXNormalizer()
        gfa_path = _write_temp(PURE_GFA)
        try:
            tmp_path, stats = norm.normalize_to_tempfile(gfa_path)
            try:
                self.assertTrue(os.path.isfile(tmp_path))
            finally:
                os.unlink(tmp_path)
        finally:
            os.unlink(gfa_path)

    def test_tempfile_stats_correct(self):
        norm = VGFRXNormalizer()
        gfa_path = _write_temp(PURE_GFA)
        try:
            tmp_path, stats = norm.normalize_to_tempfile(gfa_path)
            os.unlink(tmp_path)
            self.assertEqual(stats['s_lines'], 2)
            self.assertEqual(stats['l_lines'], 1)
        finally:
            os.unlink(gfa_path)


# ======================================================================
# run_vg_import integration
# ======================================================================

class TestRunVGImport(unittest.TestCase):
    """run_vg_import must produce a valid .frx.db from vg GFA."""

    def setUp(self):
        self.gfa_path = _write_temp(VG_GFA_10)
        fd, self.db_path = tempfile.mkstemp(suffix='.frx.db')
        os.close(fd)
        os.unlink(self.db_path)  # run_vg_import creates it fresh

    def tearDown(self):
        os.unlink(self.gfa_path)
        try:
            os.unlink(self.db_path)
        except FileNotFoundError:
            pass

    def test_run_produces_db(self):
        result = run_vg_import(
            self.gfa_path, self.db_path,
            seed=b"test_vg_import_seed_000000000000",
            verbose=False,
        )
        self.assertTrue(os.path.isfile(self.db_path))

    def test_nodes_written_matches_s_lines(self):
        result = run_vg_import(
            self.gfa_path, self.db_path,
            seed=b"test_vg_import_seed_000000000001",
            verbose=False,
        )
        self.assertEqual(result['nodes_written'], 3)

    def test_kmer_index_not_built_by_default(self):
        result = run_vg_import(
            self.gfa_path, self.db_path,
            seed=b"test_vg_import_seed_000000000002",
            verbose=False,
        )
        self.assertEqual(result['kmer_entries'], 0)

    def test_kmer_index_built_when_requested(self):
        result = run_vg_import(
            self.gfa_path, self.db_path,
            seed=b"test_vg_import_seed_000000000003",
            build_kmer_index=True,
            verbose=False,
        )
        self.assertGreater(result['kmer_entries'], 0)

    def test_normalized_stats_returned(self):
        result = run_vg_import(
            self.gfa_path, self.db_path,
            seed=b"test_vg_import_seed_000000000004",
            verbose=False,
        )
        self.assertIn('normalized_stats', result)
        self.assertEqual(result['normalized_stats']['s_lines'], 3)


if __name__ == '__main__':
    unittest.main(verbosity=2)
