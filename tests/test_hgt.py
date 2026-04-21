import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from engine import PanIndexEngine
from index import PanIndexStore
from hgt_handler import OrganismGraph, HGTSimulation


PLASMID_SEQ = "ATGCGTCGTAGCTAGCTAGCTGATCGATCGATCGATCGAATTCGCTAGCTAGCTAGCATG"


class TestCanonicalCycleHashRotation(unittest.TestCase):
    """All rotations of a plasmid sequence must yield identical canonical hashes."""

    def setUp(self):
        self.engine = PanIndexEngine()

    def test_all_rotations_same_hash(self):
        seq = PLASMID_SEQ
        n = len(seq)
        hashes = set()
        for i in range(n):
            rotated = seq[i:] + seq[:i]
            hashes.add(self.engine.canonical_cycle_hash(rotated))
        self.assertEqual(len(hashes), 1,
                         f"Expected 1 unique hash across {n} rotations, got {len(hashes)}")


class TestHGTSymlink(unittest.TestCase):
    """After HGT, the recipient must point to the donor's exact canonical hash."""

    def setUp(self):
        engine = PanIndexEngine(pangenome_seed=b"test_hgt_symlink_seed_000000000")
        store = PanIndexStore()
        self.donor = OrganismGraph("Donor", engine, store)
        self.recipient = OrganismGraph("Recipient", engine, store)
        self.engine = engine

        self.donor.add_chromosome_node("d1", "AAAA", ["core"])
        self.blatem_hash = self.donor.register_plasmid(
            "pDonor-blaTEM", PLASMID_SEQ, ["AMR:blaTEM"]
        )

        self.recipient.add_chromosome_node("r1", "CCCC", ["core"])
        self.recipient.add_chromosome_node("r2", "GGGG", ["core"])

    def test_no_resistance_before_hgt(self):
        self.assertFalse(self.recipient.has_resistance(self.blatem_hash))

    def test_has_resistance_after_hgt(self):
        self.recipient.receive_hgt("HGT:blaTEM", self.blatem_hash, "Donor")
        self.assertTrue(self.recipient.has_resistance(self.blatem_hash))

    def test_symlink_points_to_donor_hash(self):
        self.recipient.receive_hgt("HGT:blaTEM", self.blatem_hash, "Donor")
        stored_hash = self.recipient.hgt_symlinks.get("HGT:blaTEM")
        self.assertEqual(stored_hash, self.blatem_hash)

    def test_donor_unaffected_by_hgt(self):
        self.recipient.receive_hgt("HGT:blaTEM", self.blatem_hash, "Donor")
        # Donor's plasmid registry must be unchanged
        self.assertIn(self.blatem_hash, self.donor.plasmids)


class TestNodeSplitOnHGTInsertion(unittest.TestCase):
    """
    When HGT inserts into a chromosome, the split halves must:
    1. Not share the same address.
    2. Each have a valid 32-byte address.
    3. The original node must be gone from the chromosome.
    """

    def setUp(self):
        engine = PanIndexEngine(pangenome_seed=b"test_node_split_seed_000000000a")
        store = PanIndexStore()
        self.recipient = OrganismGraph("Recv", engine, store)
        self.engine = engine

        self.recipient.add_chromosome_node("chr1", "ATCGATCGATCG", ["core"])
        self.recipient.add_chromosome_node("chr2", "GCTAGCTAGCTA", ["core"])

        donor_engine = PanIndexEngine()
        donor_store = PanIndexStore()
        donor = OrganismGraph("Donor", donor_engine, donor_store)
        self.blatem_hash = donor.register_plasmid(
            "pD-blaTEM", PLASMID_SEQ, ["AMR:blaTEM"]
        )

    def test_original_node_removed_after_split(self):
        self.recipient.receive_hgt("HGT:blaTEM", self.blatem_hash, "Donor",
                                   insertion_after="chr2")
        self.assertNotIn("chr2", self.recipient.chromosome)

    def test_split_halves_exist(self):
        self.recipient.receive_hgt("HGT:blaTEM", self.blatem_hash, "Donor",
                                   insertion_after="chr2")
        self.assertIn("chr2_L", self.recipient.chromosome)
        self.assertIn("chr2_R", self.recipient.chromosome)

    def test_split_halves_have_unique_addresses(self):
        self.recipient.receive_hgt("HGT:blaTEM", self.blatem_hash, "Donor",
                                   insertion_after="chr2")
        addr_l = self.recipient.chromosome["chr2_L"]["addr"]
        addr_r = self.recipient.chromosome["chr2_R"]["addr"]
        self.assertNotEqual(addr_l, addr_r)

    def test_split_halves_are_32_bytes(self):
        self.recipient.receive_hgt("HGT:blaTEM", self.blatem_hash, "Donor",
                                   insertion_after="chr2")
        self.assertEqual(len(self.recipient.chromosome["chr2_L"]["addr"]), 32)
        self.assertEqual(len(self.recipient.chromosome["chr2_R"]["addr"]), 32)


class TestGlobalAMRScan(unittest.TestCase):
    """Global tag scan must find all organisms carrying a resistance gene."""

    def test_multiple_recipients_detected(self):
        sim = HGTSimulation()
        blatem_hash, _ = sim.run()

        amr_nodes = sim.global_store.lookup_by_tag("AMR:blaTEM")
        self.assertTrue(len(amr_nodes) >= 1)

    def test_hash_is_same_content_same_organism(self):
        # Two organisms registering the exact same plasmid sequence
        # must get the exact same canonical hash
        engine = PanIndexEngine()
        h1 = engine.canonical_cycle_hash(PLASMID_SEQ)
        h2 = engine.canonical_cycle_hash(PLASMID_SEQ)
        self.assertEqual(h1, h2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
