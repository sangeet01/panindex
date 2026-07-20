import sys
import os
import hashlib
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore


class OrganismGraph:
    """
    Represents one bacterium's complete genome as a PanIndex graph.
    Holds a chromosome DAG plus a registry of plasmid symlinks.
    """

    def __init__(self, name: str, engine: PanIndexEngine, store: PanIndexStore):
        self.name = name
        self.engine = engine
        self.store = store
        # Chromosome nodes: node_id -> {seq, addr, tags}
        self.chromosome: Dict[str, dict] = {}
        # Plasmid registry: plasmid_canonical_hash -> metadata
        self.plasmids: Dict[bytes, dict] = {}
        # Topological symlinks: symlink_label -> target_hash
        self.hgt_symlinks: Dict[str, bytes] = {}

    def add_chromosome_node(self, node_id: str, seq: str, tags: List[str]):
        """Add a segment to the chromosome DAG and index it."""
        parent = (
            list(self.chromosome.values())[-1]['addr']
            if self.chromosome
            else self.engine.root_hash
        )
        derived = self.engine.derive_ratchet_address(parent, f"{self.name}/{node_id}")
        addr = self.engine.compute_node_address(seq, derived)
        self.chromosome[node_id] = {'seq': seq, 'addr': addr, 'tags': tags}
        self.store.insert(
            node_id, addr, tags,
            {'seq': seq, 'organism': self.name},
            organism=self.name,
        )

    def register_plasmid(self, plasmid_name: str, circular_sequence: str,
                         resistance_tags: List[str]) -> bytes:
        """
        Register a circular plasmid. Returns its canonical hash.
        Canonical Cycle Hash makes the hash rotation-invariant.
        """
        canonical_hash = self.engine.canonical_cycle_hash(circular_sequence)
        self.plasmids[canonical_hash] = {
            'name': plasmid_name,
            'sequence': circular_sequence,
            'tags': resistance_tags,
        }
        # Index the plasmid by its canonical hash and tags
        self.store.insert(
            node_id=f"PLASMID:{plasmid_name}",
            address=canonical_hash,
            tags=resistance_tags,
            metadata={'type': 'plasmid', 'organism': self.name, 'seq': circular_sequence},
            organism=self.name,
        )
        return canonical_hash

    def receive_hgt(self, symlink_label: str, donor_plasmid_hash: bytes,
                    donor_name: str, insertion_after: Optional[str] = None):
        """
        Model a Horizontal Gene Transfer event.

        Instead of copying the sequence, we create a Topological Symlink:
        a pointer from this organism's graph to the donor's plasmid hash.

        Args:
            symlink_label       : Human-readable label for this HGT event.
            donor_plasmid_hash  : The canonical hash of the donor plasmid.
            donor_name          : Name of the donor organism.
            insertion_after     : Chromosome node_id where the insertion occurs
                                  (triggers node split simulation).
        """
        self.hgt_symlinks[symlink_label] = donor_plasmid_hash
        print(f"  [{self.name}] HGT received: '{symlink_label}' -> "
              f"{donor_plasmid_hash.hex()[:12]}... (from {donor_name})")

        if insertion_after and insertion_after in self.chromosome:
            self._simulate_node_split(insertion_after, symlink_label, donor_plasmid_hash)

    def _simulate_node_split(self, node_id: str, hgt_label: str,
                             hgt_hash: bytes):
        """
        When HGT inserts into a chromosome segment, the Fractal Ratchet
        splits the node. Addresses of split halves remain independently stable.
        """
        node = self.chromosome[node_id]
        seq = node['seq']
        mid = len(seq) // 2

        seq_left = seq[:mid]
        seq_right = seq[mid:]

        # Derive addresses for the two halves
        parent_addr = node['addr']
        left_addr = self.engine.derive_ratchet_address(parent_addr, f"{node_id}_L")
        right_addr = self.engine.derive_ratchet_address(parent_addr, f"{node_id}_R")

        left_final = self.engine.compute_node_address(seq_left, left_addr)
        right_final = self.engine.compute_node_address(seq_right, right_addr)

        # Split in chromosome registry
        self.chromosome[f"{node_id}_L"] = {'seq': seq_left, 'addr': left_final, 'tags': node['tags']}
        self.chromosome[f"{node_id}_R"] = {'seq': seq_right, 'addr': right_final, 'tags': node['tags']}
        del self.chromosome[node_id]

        # Index split halves
        self.store.insert(f"{node_id}_L", left_final, node['tags'],
                          {'seq': seq_left, 'organism': self.name, 'split_from': node_id})
        self.store.insert(f"{node_id}_R", right_final, node['tags'],
                          {'seq': seq_right, 'organism': self.name, 'split_from': node_id})

        print(f"  [{self.name}] Node '{node_id}' split: "
              f"'{node_id}_L'={left_final.hex()[:10]}... | "
              f"HGT:{hgt_hash.hex()[:10]}... | "
              f"'{node_id}_R'={right_final.hex()[:10]}...")

    def has_resistance(self, resistance_hash: bytes) -> bool:
        """
        Check if this organism carries a known resistance gene.
        O(1) dictionary lookup - no alignment required.
        """
        return resistance_hash in self.plasmids or any(
            v == resistance_hash for v in self.hgt_symlinks.values()
        )


class HGTSimulation:
    """
    Runs a full HGT simulation between two bacteria sharing a resistance plasmid.
    """

    def __init__(self):
        # Shared index across all organisms (simulates a global AMR database)
        self.global_engine = PanIndexEngine(pangenome_seed=b"global_pangenome_seed_000000000")
        self.global_store = PanIndexStore()

    def run(self):
        print("=" * 60)
        print("PanIndex HGT Simulation: blaTEM-1 Resistance Transfer")
        print("=" * 60)

        # --- Donor: Klebsiella pneumoniae strain KP001 ---
        print("\n[Step 1] Building donor organism: K. pneumoniae KP001")
        kp001 = OrganismGraph("KP001", self.global_engine, self.global_store)
        kp001.add_chromosome_node("chrA", "ATGCATGCATGC", ["core_chromosome"])
        kp001.add_chromosome_node("chrB", "GCTAGCTAGCTA", ["core_chromosome"])
        kp001.add_chromosome_node("chrC", "TTAATTAATTAA", ["core_chromosome"])

        # Register a resistance plasmid (circular)
        # Note: sequence deliberately starts at an arbitrary rotation point
        plasmid_seq = "ATGCGTCGTAGCTAGCTAGCTGATCGATCGATCGATCGAATTCGCTAGCTAGCTAGCATG"
        blatem_hash = kp001.register_plasmid(
            "pKP001-blaTEM",
            plasmid_seq,
            ["AMR:blaTEM", "resistance:beta-lactam", "mobile_element"]
        )
        print(f"  [KP001] Plasmid registered. Canonical hash: {blatem_hash.hex()[:16]}...")

        # --- Recipient: Escherichia coli strain EC042 ---
        print("\n[Step 2] Building recipient organism: E. coli EC042 (no resistance)")
        ec042 = OrganismGraph("EC042", self.global_engine, self.global_store)
        ec042.add_chromosome_node("chr1", "CCCCGGGGAAAA", ["core_chromosome"])
        ec042.add_chromosome_node("chr2", "TTTTAAAACCCC", ["core_chromosome"])
        ec042.add_chromosome_node("chr3", "GGGGTTTTCCCC", ["core_chromosome"])

        print(f"  [EC042] Has blaTEM before HGT: {ec042.has_resistance(blatem_hash)}")

        # --- HGT Event ---
        print("\n[Step 3] Simulating HGT conjugation event: KP001 -> EC042")
        ec042.receive_hgt(
            symlink_label="HGT:pKP001-blaTEM",
            donor_plasmid_hash=blatem_hash,
            donor_name="KP001",
            insertion_after="chr2"  # Insertion into chromosome triggers split
        )

        # --- Detection ---
        print("\n[Step 4] Scanning EC042 for known resistance genes (O(1) lookup)")
        has_resistance = ec042.has_resistance(blatem_hash)
        print(f"  [EC042] Has blaTEM after HGT: {has_resistance}")

        # Global AMR scan: check all organisms in the global store
        print("\n[Step 5] Global AMR scan via tag index (O(1) per organism)")
        amr_nodes = self.global_store.lookup_by_tag("AMR:blaTEM")
        print(f"  Organisms/Nodes carrying AMR:blaTEM -> {amr_nodes}")

        print("\n[SUMMARY]")
        print(f"  Donor plasmid canonical hash : {blatem_hash.hex()[:24]}...")
        print(f"  HGT symlink in EC042         : {list(ec042.hgt_symlinks.keys())}")
        print(f"  Detection method             : O(1) hash dictionary lookup")
        print(f"  Total indexed nodes          : {len(self.global_store)}")
        print("\n  SUCCESS: Full HGT simulation complete.")
        print("=" * 60)
        return blatem_hash, ec042


if __name__ == "__main__":
    sim = HGTSimulation()
    sim.run()
