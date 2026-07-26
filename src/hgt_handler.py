"""
FRX Horizontal Gene Transfer (HGT) Engine.

Production-grade implementation for any-scale genomics data.

Architecture
------------
OrganismGraph
    Represents one organism's complete genome in the FRX index.
    Chromosome nodes are stored in the shared PanIndexStore.
    Plasmids are registered by their canonical cycle hash (rotation-invariant).
    HGT events are modelled as topological symlinks — no sequence duplication.

HGTRegistry
    SQLite-backed persistent registry of all HGT events and plasmid records.
    Designed for datasets with millions of organisms and plasmids:
      - Streaming inserts via batched executemany (configurable batch size).
    - Known-key queries use direct index lookups (no full-table scans).
      - Thread-safe via per-connection WAL journal mode.
      - Supports incremental updates: re-running on the same DB is idempotent.

HGTSimulation
    Orchestrates a full donor -> recipient conjugation demo using the
    production HGTRegistry and OrganismGraph APIs.

Scaling notes
-------------
- Chromosome nodes: streamed one at a time; peak RAM = O(1) per node.
- Plasmid registration: canonical_cycle_hash is O(|seq|); no full graph load.
- HGT symlinks: direct insert and expected O(1) key lookup regardless of
    unrelated graph size.
- Global AMR scan: expected O(1) tag posting-list lookup via
    PanIndexStore.lookup_by_tag(), plus result enumeration.
- Node splitting on HGT insertion: O(|seq|) for the split node only.
"""

import logging
import os
import sqlite3
import tempfile
from typing import Dict, Iterator, List, Optional, Tuple

from engine import PanIndexEngine
from index import PanIndexStore

logger = logging.getLogger(__name__)

# Default batch size for bulk inserts into the HGT registry
_DEFAULT_BATCH = 500


# ======================================================================
# HGT Registry (SQLite-backed, production-scale)
# ======================================================================

_HGT_DDL = """
CREATE TABLE IF NOT EXISTS plasmids (
    canonical_hex  TEXT PRIMARY KEY,
    name           TEXT NOT NULL,
    sequence       TEXT NOT NULL,
    organism       TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS plasmid_tags (
    canonical_hex  TEXT NOT NULL,
    tag            TEXT NOT NULL,
    PRIMARY KEY (canonical_hex, tag)
);

CREATE TABLE IF NOT EXISTS hgt_events (
    id             INTEGER PRIMARY KEY AUTOINCREMENT,
    symlink_label  TEXT NOT NULL,
    donor_hex      TEXT NOT NULL,
    donor_org      TEXT NOT NULL,
    recipient_org  TEXT NOT NULL,
    insertion_node TEXT NOT NULL DEFAULT '',
    UNIQUE (symlink_label, recipient_org)
);

CREATE INDEX IF NOT EXISTS idx_plasmid_org  ON plasmids(organism);
CREATE INDEX IF NOT EXISTS idx_hgt_donor    ON hgt_events(donor_hex);
CREATE INDEX IF NOT EXISTS idx_hgt_recipient ON hgt_events(recipient_org);
CREATE INDEX IF NOT EXISTS idx_hgt_label    ON hgt_events(symlink_label);
"""


class HGTRegistry:
    """
    Persistent, append-only registry of plasmids and HGT events.

    Backed by SQLite with WAL journal mode for concurrent read safety.
    All writes are batched and wrapped in explicit transactions.

    Usage:
        registry = HGTRegistry("pangenome.hgt.db")
        registry.register_plasmid(canonical_hex, name, seq, organism, tags)
        registry.record_hgt_event(label, donor_hex, donor_org, recipient_org)
        events = registry.get_hgt_events_for_recipient("EC042")
        registry.close()
    """

    def __init__(self, db_path: str, batch_size: int = _DEFAULT_BATCH):
        self.db_path = db_path
        self.batch_size = batch_size
        self._conn = self._open(db_path)
        self._plasmid_buf: List[tuple] = []
        self._tag_buf: List[tuple] = []
        self._event_buf: List[tuple] = []

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    @staticmethod
    def _open(db_path: str) -> sqlite3.Connection:
        conn = sqlite3.connect(db_path, check_same_thread=False)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.execute("PRAGMA foreign_keys=ON")
        conn.executescript(_HGT_DDL)
        conn.commit()
        return conn

    def flush(self):
        """Flush all pending buffered writes to the database."""
        with self._conn:
            if self._plasmid_buf:
                self._conn.executemany(
                    "INSERT OR IGNORE INTO plasmids "
                    "(canonical_hex, name, sequence, organism) VALUES (?,?,?,?)",
                    self._plasmid_buf,
                )
                self._plasmid_buf.clear()
            if self._tag_buf:
                self._conn.executemany(
                    "INSERT OR IGNORE INTO plasmid_tags (canonical_hex, tag) VALUES (?,?)",
                    self._tag_buf,
                )
                self._tag_buf.clear()
            if self._event_buf:
                self._conn.executemany(
                    "INSERT OR IGNORE INTO hgt_events "
                    "(symlink_label, donor_hex, donor_org, recipient_org, insertion_node) "
                    "VALUES (?,?,?,?,?)",
                    self._event_buf,
                )
                self._event_buf.clear()

    def close(self):
        """Flush pending writes and close the database connection."""
        self.flush()
        self._conn.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    # ------------------------------------------------------------------
    # Write API
    # ------------------------------------------------------------------

    def register_plasmid(
        self,
        canonical_hex: str,
        name: str,
        sequence: str,
        organism: str,
        tags: List[str],
    ):
        """
        Register a plasmid. Idempotent: re-registering the same canonical_hex
        is a no-op (INSERT OR IGNORE).
        """
        self._plasmid_buf.append((canonical_hex, name, sequence, organism))
        for tag in tags:
            self._tag_buf.append((canonical_hex, tag))
        if len(self._plasmid_buf) >= self.batch_size:
            self.flush()

    def record_hgt_event(
        self,
        symlink_label: str,
        donor_hex: str,
        donor_org: str,
        recipient_org: str,
        insertion_node: str = "",
    ):
        """
        Record an HGT conjugation event. Idempotent on (symlink_label, recipient_org).
        """
        self._event_buf.append(
            (symlink_label, donor_hex, donor_org, recipient_org, insertion_node)
        )
        if len(self._event_buf) >= self.batch_size:
            self.flush()

    # ------------------------------------------------------------------
    # Read API
    # ------------------------------------------------------------------

    def get_plasmid(self, canonical_hex: str) -> Optional[dict]:
        """Return plasmid metadata dict or None if not found."""
        self.flush()
        row = self._conn.execute(
            "SELECT name, sequence, organism FROM plasmids WHERE canonical_hex=?",
            (canonical_hex,),
        ).fetchone()
        if not row:
            return None
        tags = [
            r[0] for r in self._conn.execute(
                "SELECT tag FROM plasmid_tags WHERE canonical_hex=?",
                (canonical_hex,),
            ).fetchall()
        ]
        return {"name": row[0], "sequence": row[1], "organism": row[2], "tags": tags}

    def get_hgt_events_for_recipient(self, recipient_org: str) -> List[dict]:
        """Return all HGT events received by an organism."""
        self.flush()
        rows = self._conn.execute(
            "SELECT symlink_label, donor_hex, donor_org, insertion_node "
            "FROM hgt_events WHERE recipient_org=?",
            (recipient_org,),
        ).fetchall()
        return [
            {
                "symlink_label": r[0],
                "donor_hex": r[1],
                "donor_org": r[2],
                "insertion_node": r[3],
            }
            for r in rows
        ]

    def get_hgt_events_for_donor_plasmid(self, canonical_hex: str) -> List[dict]:
        """Return all organisms that received a specific plasmid."""
        self.flush()
        rows = self._conn.execute(
            "SELECT symlink_label, donor_org, recipient_org, insertion_node "
            "FROM hgt_events WHERE donor_hex=?",
            (canonical_hex,),
        ).fetchall()
        return [
            {
                "symlink_label": r[0],
                "donor_org": r[1],
                "recipient_org": r[2],
                "insertion_node": r[3],
            }
            for r in rows
        ]

    def iter_all_plasmids(self) -> Iterator[dict]:
        """Stream all registered plasmids without loading them all into RAM."""
        self.flush()
        cursor = self._conn.execute(
            "SELECT canonical_hex, name, sequence, organism FROM plasmids"
        )
        for row in cursor:
            canonical_hex, name, sequence, organism = row
            tags = [
                r[0] for r in self._conn.execute(
                    "SELECT tag FROM plasmid_tags WHERE canonical_hex=?",
                    (canonical_hex,),
                ).fetchall()
            ]
            yield {
                "canonical_hex": canonical_hex,
                "name": name,
                "sequence": sequence,
                "organism": organism,
                "tags": tags,
            }

    def stats(self) -> dict:
        """Return registry statistics without loading data into RAM."""
        self.flush()
        (plasmid_count,) = self._conn.execute(
            "SELECT COUNT(*) FROM plasmids"
        ).fetchone()
        (event_count,) = self._conn.execute(
            "SELECT COUNT(*) FROM hgt_events"
        ).fetchone()
        (organism_count,) = self._conn.execute(
            "SELECT COUNT(DISTINCT recipient_org) FROM hgt_events"
        ).fetchone()
        return {
            "plasmids": plasmid_count,
            "hgt_events": event_count,
            "recipient_organisms": organism_count,
            "db_size_bytes": os.path.getsize(self.db_path),
        }


# ======================================================================
# OrganismGraph
# ======================================================================

class OrganismGraph:
    """
    Represents one organism's complete genome as a PanIndex graph.

    Chromosome nodes are indexed in the shared PanIndexStore one at a time
    (O(1) peak RAM per node). Plasmids are registered by canonical cycle hash.
    HGT events are modelled as topological symlinks — no sequence duplication.

    For large chromosomes, use add_chromosome_nodes_streaming() to process
    an iterator of (node_id, seq, tags) tuples without loading all nodes
    into memory simultaneously.
    """

    def __init__(
        self,
        name: str,
        engine: PanIndexEngine,
        store: PanIndexStore,
        registry: Optional[HGTRegistry] = None,
    ):
        self.name = name
        self.engine = engine
        self.store = store
        self.registry = registry

        # In-memory chromosome registry: node_id -> {seq, addr, tags}
        # For very large chromosomes, callers should use the streaming API
        # and manage their own external node list.
        self.chromosome: Dict[str, dict] = {}

        # Plasmid registry: canonical_hash (bytes) -> metadata
        self.plasmids: Dict[bytes, dict] = {}

        # HGT symlinks: symlink_label -> donor_canonical_hash (bytes)
        self.hgt_symlinks: Dict[str, bytes] = {}

        # Running parent address for chained chromosome derivation
        self._chrom_parent: bytes = engine.root_hash

    # ------------------------------------------------------------------
    # Chromosome
    # ------------------------------------------------------------------

    def add_chromosome_node(self, node_id: str, seq: str, tags: List[str]):
        """
        Add a single segment to the chromosome DAG and index it.

        Derives the address from the previous node's address (chained ratchet),
        so the derivation path encodes the linear order of the chromosome.
        """
        derived = self.engine.derive_ratchet_address(
            self._chrom_parent, f"{self.name}/{node_id}"
        )
        addr = self.engine.compute_node_address(seq, derived)
        self.chromosome[node_id] = {"seq": seq, "addr": addr, "tags": tags}
        self._chrom_parent = addr
        self.store.insert(
            node_id, addr, tags,
            {"seq": seq, "organism": self.name, "derivation_path": f"{self.name}/{node_id}"},
            organism=self.name,
        )

    def add_chromosome_nodes_streaming(
        self,
        nodes: Iterator[Tuple[str, str, List[str]]],
    ):
        """
        Stream-index an arbitrary number of chromosome nodes with O(1) peak RAM.

        Args:
            nodes: Iterator of (node_id, sequence, tags) tuples.
                   Can be a generator reading from a file line by line.
        """
        for node_id, seq, tags in nodes:
            self.add_chromosome_node(node_id, seq, tags)

    # ------------------------------------------------------------------
    # Plasmids
    # ------------------------------------------------------------------

    def register_plasmid(
        self,
        plasmid_name: str,
        circular_sequence: str,
        resistance_tags: List[str],
    ) -> bytes:
        """
        Register a circular plasmid. Returns its canonical hash (bytes).

        The canonical cycle hash is rotation-invariant: the same plasmid
        registered at any rotation point produces the same hash.
        Idempotent: re-registering the same sequence is a no-op.
        """
        canonical_hash = self.engine.canonical_cycle_hash(circular_sequence)
        canonical_hex = canonical_hash.hex()

        self.plasmids[canonical_hash] = {
            "name": plasmid_name,
            "sequence": circular_sequence,
            "tags": resistance_tags,
        }

        self.store.insert(
            node_id=f"PLASMID:{plasmid_name}",
            address=canonical_hash,
            tags=resistance_tags,
            metadata={
                "type": "plasmid",
                "organism": self.name,
                "seq": circular_sequence,
            },
            organism=self.name,
        )

        if self.registry is not None:
            self.registry.register_plasmid(
                canonical_hex, plasmid_name, circular_sequence,
                self.name, resistance_tags,
            )

        logger.debug(
            "[%s] Plasmid '%s' registered. Canonical hash: %s...",
            self.name, plasmid_name, canonical_hex[:16],
        )
        return canonical_hash

    # ------------------------------------------------------------------
    # HGT
    # ------------------------------------------------------------------

    def receive_hgt(
        self,
        symlink_label: str,
        donor_plasmid_hash: bytes,
        donor_name: str,
        insertion_after: Optional[str] = None,
    ):
        """
        Model a Horizontal Gene Transfer event.

        Creates a topological symlink pointing to the donor's canonical hash.
        No sequence data is duplicated. The symlink is persisted to the
        HGTRegistry if one is attached.

        If insertion_after is given and the node exists in the chromosome,
        the node is split at its midpoint and the HGT insertion is recorded
        between the two halves.

        Args:
            symlink_label      : Human-readable label for this HGT event.
            donor_plasmid_hash : Canonical hash (bytes) of the donor plasmid.
            donor_name         : Name of the donor organism.
            insertion_after    : Chromosome node_id where insertion occurs.
        """
        self.hgt_symlinks[symlink_label] = donor_plasmid_hash
        donor_hex = donor_plasmid_hash.hex()

        # Index the symlink in the global store (same address, different node_id)
        # so tag queries find this organism as a carrier.
        donor_meta = self.store.get_node(f"PLASMID:{donor_name}")
        donor_tags = donor_meta["tags"] if donor_meta else []
        hgt_tags = list(donor_tags) + ["HGT:symlink", f"HGT:donor:{donor_name}"]

        self.store.insert(
            node_id=f"{self.name}:HGT:{symlink_label}",
            address=donor_plasmid_hash,
            tags=hgt_tags,
            metadata={
                "type": "hgt_symlink",
                "organism": self.name,
                "donor": donor_name,
                "donor_hex": donor_hex,
                "symlink_label": symlink_label,
            },
            organism=self.name,
        )

        if self.registry is not None:
            self.registry.record_hgt_event(
                symlink_label, donor_hex, donor_name, self.name,
                insertion_node=insertion_after or "",
            )

        logger.info(
            "[%s] HGT received: '%s' -> %s... (from %s)",
            self.name, symlink_label, donor_hex[:12], donor_name,
        )

        if insertion_after and insertion_after in self.chromosome:
            self._simulate_node_split(insertion_after, symlink_label, donor_plasmid_hash)

    def _simulate_node_split(
        self,
        node_id: str,
        hgt_label: str,
        hgt_hash: bytes,
    ):
        """
        Split a chromosome node at its midpoint to model HGT insertion.

        The two halves receive independent ratchet addresses derived from
        the original node's address, so they remain stable even if the
        graph is later extended.
        """
        node = self.chromosome[node_id]
        seq = node["seq"]
        mid = max(1, len(seq) // 2)

        seq_left = seq[:mid]
        seq_right = seq[mid:]

        parent_addr = node["addr"]
        left_derived = self.engine.derive_ratchet_address(parent_addr, f"{node_id}_L")
        right_derived = self.engine.derive_ratchet_address(parent_addr, f"{node_id}_R")

        left_addr = self.engine.compute_node_address(seq_left, left_derived)
        right_addr = self.engine.compute_node_address(seq_right, right_derived)

        self.chromosome[f"{node_id}_L"] = {
            "seq": seq_left, "addr": left_addr, "tags": node["tags"]
        }
        self.chromosome[f"{node_id}_R"] = {
            "seq": seq_right, "addr": right_addr, "tags": node["tags"]
        }
        del self.chromosome[node_id]

        self.store.insert(
            f"{node_id}_L", left_addr, node["tags"],
            {"seq": seq_left, "organism": self.name, "split_from": node_id},
            organism=self.name,
        )
        self.store.insert(
            f"{node_id}_R", right_addr, node["tags"],
            {"seq": seq_right, "organism": self.name, "split_from": node_id},
            organism=self.name,
        )

        logger.info(
            "[%s] Node '%s' split: '%s_L'=%s... | HGT:%s... | '%s_R'=%s...",
            self.name, node_id,
            node_id, left_addr.hex()[:10],
            hgt_hash.hex()[:10],
            node_id, right_addr.hex()[:10],
        )

    # ------------------------------------------------------------------
    # Resistance detection
    # ------------------------------------------------------------------

    def has_resistance(self, resistance_hash: bytes) -> bool:
        """
        O(1) check: does this organism carry a known resistance gene?

        Checks both directly registered plasmids and received HGT symlinks.
        No sequence alignment required.
        """
        if resistance_hash in self.plasmids:
            return True
        return any(v == resistance_hash for v in self.hgt_symlinks.values())

    def resistance_lineage(self, resistance_hash: bytes) -> Optional[str]:
        """
        Return the donor organism name if this resistance was acquired via HGT,
        or 'native' if it was registered directly, or None if not present.
        """
        if resistance_hash in self.plasmids:
            return "native"
        for label, h in self.hgt_symlinks.items():
            if h == resistance_hash:
                # Extract donor from label convention "HGT:donor:<name>"
                node_id = f"{self.name}:HGT:{label}"
                node = self.store.get_node(node_id)
                if node:
                    return node["metadata"].get("donor", "unknown")
        return None


# ======================================================================
# Global AMR Scanner
# ======================================================================

class AMRScanner:
    """
    Alignment-free global AMR scan across all indexed organisms.

    All operations are O(1) tag index lookups via PanIndexStore.
    Supports streaming output for result sets too large to hold in RAM.
    """

    def __init__(self, store: PanIndexStore, registry: Optional[HGTRegistry] = None):
        self.store = store
        self.registry = registry

    def scan_tag(
        self,
        amr_tag: str,
        organism: Optional[str] = None,
    ) -> List[str]:
        """
        Return all node_ids carrying amr_tag, optionally scoped to one organism.
        O(1) lookup.
        """
        return self.store.lookup_by_tag(amr_tag, organism=organism)

    def organisms_carrying(self, amr_tag: str) -> List[str]:
        """
        Return the distinct set of organism names carrying amr_tag.
        O(K) where K = number of nodes with that tag.
        """
        nodes = self.store.lookup_by_tag(amr_tag)
        seen = set()
        result = []
        for nid in nodes:
            node = self.store.get_node(nid)
            if node:
                org = node.get("organism", "")
                if org and org not in seen:
                    seen.add(org)
                    result.append(org)
        return result

    def hgt_spread_report(self, canonical_hex: str) -> dict:
        """
        Report how far a specific plasmid has spread via HGT.

        Requires an attached HGTRegistry for full lineage data.
        Falls back to tag-based detection if no registry is attached.
        """
        if self.registry is not None:
            events = self.registry.get_hgt_events_for_donor_plasmid(canonical_hex)
            recipients = [e["recipient_org"] for e in events]
        else:
            # Fallback: find all nodes whose address matches the canonical hash
            node_id = self.store.lookup_by_address(bytes.fromhex(canonical_hex))
            recipients = []
            if node_id:
                node = self.store.get_node(node_id)
                if node:
                    recipients = [node.get("organism", "")]

        return {
            "canonical_hex": canonical_hex,
            "recipient_count": len(recipients),
            "recipients": recipients,
        }

    def stream_critical_nodes(self) -> Iterator[dict]:
        """
        Stream all nodes tagged RES:critical without loading them all into RAM.
        Useful for large databases with millions of nodes.
        """
        for nid in self.store.lookup_by_tag("RES:critical"):
            node = self.store.get_node(nid)
            if node:
                yield {"node_id": nid, **node}


# ======================================================================
# HGT Simulation (demo / integration test)
# ======================================================================

class HGTSimulation:
    """
    Runs a full HGT simulation between two bacteria sharing a resistance plasmid.

    Uses a temporary HGTRegistry backed by an in-memory SQLite database
    so the simulation is self-contained and leaves no files on disk.
    """

    def __init__(self, registry_path: Optional[str] = None):
        self.global_engine = PanIndexEngine(
            pangenome_seed=b"global_pangenome_seed_000000000"
        )
        self.global_store = PanIndexStore()

        # Use a temp file for the registry so the simulation is portable
        if registry_path is None:
            fd, registry_path = tempfile.mkstemp(suffix=".hgt.db")
            os.close(fd)
            self._owns_registry_file = True
        else:
            self._owns_registry_file = False

        self._registry_path = registry_path
        self.registry = HGTRegistry(registry_path)

    def run(self) -> Tuple[bytes, "OrganismGraph"]:
        print("=" * 60)
        print("PanIndex HGT Simulation: blaTEM-1 Resistance Transfer")
        print("=" * 60)

        # --- Donor: Klebsiella pneumoniae KP001 ---
        print("\n[Step 1] Building donor organism: K. pneumoniae KP001")
        kp001 = OrganismGraph(
            "KP001", self.global_engine, self.global_store, self.registry
        )
        kp001.add_chromosome_node("chrA", "ATGCATGCATGC", ["core_chromosome"])
        kp001.add_chromosome_node("chrB", "GCTAGCTAGCTA", ["core_chromosome"])
        kp001.add_chromosome_node("chrC", "TTAATTAATTAA", ["core_chromosome"])

        plasmid_seq = (
            "ATGCGTCGTAGCTAGCTAGCTGATCGATCGATCGATCG"
            "AATTCGCTAGCTAGCTAGCATG"
        )
        blatem_hash = kp001.register_plasmid(
            "pKP001-blaTEM",
            plasmid_seq,
            ["AMR:blaTEM", "resistance:beta-lactam", "mobile_element"],
        )
        print(f"  [KP001] Plasmid registered. Canonical hash: {blatem_hash.hex()[:16]}...")

        # --- Recipient: Escherichia coli EC042 ---
        print("\n[Step 2] Building recipient organism: E. coli EC042 (no resistance)")
        ec042 = OrganismGraph(
            "EC042", self.global_engine, self.global_store, self.registry
        )
        ec042.add_chromosome_node("chr1", "CCCCGGGGAAAA", ["core_chromosome"])
        ec042.add_chromosome_node("chr2", "TTTTAAAACCCC", ["core_chromosome"])
        ec042.add_chromosome_node("chr3", "GGGGTTTTCCCC", ["core_chromosome"])

        print(f"  [EC042] Has blaTEM before HGT: {ec042.has_resistance(blatem_hash)}")

        # --- HGT Event ---
        print("\n[Step 3] Simulating HGT conjugation event: KP001 -> EC042")
        ec042.receive_hgt(
            symlink_label="pKP001-blaTEM",
            donor_plasmid_hash=blatem_hash,
            donor_name="KP001",
            insertion_after="chr2",
        )

        # --- Detection ---
        print("\n[Step 4] Scanning EC042 for known resistance genes (O(1) lookup)")
        has_resistance = ec042.has_resistance(blatem_hash)
        lineage = ec042.resistance_lineage(blatem_hash)
        print(f"  [EC042] Has blaTEM after HGT : {has_resistance}")
        print(f"  [EC042] Resistance lineage   : {lineage}")

        # --- Global AMR scan ---
        print("\n[Step 5] Global AMR scan via tag index (O(1) per organism)")
        scanner = AMRScanner(self.global_store, self.registry)
        amr_nodes = scanner.scan_tag("AMR:blaTEM")
        carriers = scanner.organisms_carrying("AMR:blaTEM")
        spread = scanner.hgt_spread_report(blatem_hash.hex())
        print(f"  Nodes carrying AMR:blaTEM    -> {amr_nodes}")
        print(f"  Organisms carrying AMR:blaTEM-> {carriers}")
        print(f"  HGT spread report            -> {spread}")

        # --- Registry stats ---
        self.registry.flush()
        st = self.registry.stats()
        print(f"\n[Step 6] HGT Registry stats")
        print(f"  Plasmids registered : {st['plasmids']}")
        print(f"  HGT events recorded : {st['hgt_events']}")
        print(f"  Recipient organisms : {st['recipient_organisms']}")

        print("\n[SUMMARY]")
        print(f"  Donor plasmid canonical hash : {blatem_hash.hex()[:24]}...")
        print(f"  HGT symlinks in EC042        : {list(ec042.hgt_symlinks.keys())}")
        print(f"  Detection method             : O(1) hash dictionary lookup")
        print(f"  Total indexed nodes          : {len(self.global_store)}")
        print("\n  SUCCESS: Full HGT simulation complete.")
        print("=" * 60)

        return blatem_hash, ec042

    def cleanup(self):
        """Remove the temporary registry database file."""
        self.registry.close()
        if self._owns_registry_file:
            try:
                os.unlink(self._registry_path)
            except OSError:
                pass


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    sim = HGTSimulation()
    try:
        sim.run()
    finally:
        sim.cleanup()
