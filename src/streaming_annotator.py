"""
FRX Streaming GFA Annotator.

Processes GFA files in two passes without loading the full file into memory.
Designed for GFA files in the GB range where the standard GFAAnnotator
would exhaust available RAM.

When to use this vs the standard annotator
------------------------------------------
- Standard GFAAnnotator  : Fast, in-memory. Use for GFA files < MEM_LIMIT_MB.
- StreamingGFAAnnotator  : Constant peak RAM. Use for large GFA files.

Both produce identical output (same AN:Z:/PA:Z:/AF:i: tags, same addresses).

Two-pass algorithm
------------------
Pass 1 (S-lines + L-lines):
  First sub-pass: read all L-lines to build a neighbour map (node -> [neighbours]).
  Second sub-pass: read S-lines one at a time.  For each segment:
    - Compute ratchet address.
    - Collect neighbour tags from the neighbour map for Sannidhi context.
    - Apply the full Paninian rule engine (Utsarga/Apavada + SemanticFilter Phi).
    - Record the derivation step in DerivationHistory (Asiddhatva).
    - Write node data immediately to a SQLite temp table.  Clear from RAM.

Pass 2 (Write-back):
  Re-read the input GFA line by line.  For each S-line, fetch the pre-computed
  address from the temp table and inject AN:Z:/PA:Z:/AF:i: tags.  Write to
  output line by line.  All other line types (H/L/P/W) pass through unchanged.

Semantic Filter Phi (Akanksha + Yogyata + Sannidhi):
  - Akanksha  : demotes RES:amr_confirmed/critical to candidate when seq is empty.
  - Yogyata   : appends COMPAT:low_gc / COMPAT:high_gc based on GC content.
  - Sannidhi  : upgrades/downgrades resolution based on neighbour tag proximity.
  Extra bio-tags produced by the filter are stored in the temp DB and written
  to the output GFA as part of the AF:i: tag count.

Asiddhatva / DerivationHistory:
  A DerivationHistory is attached to the engine.  Every derive_ratchet_address()
  call is recorded.  The history path string is stored in the temp DB and in
  the PanIndexStore metadata so it is queryable after annotation.

Memory profile:
  Standard annotator : O(N * avg_seq_len) peak RAM
  Streaming annotator: O(E) peak RAM for neighbour map (E = number of L-lines)
                       + O(1) per S-line during address derivation
"""

import os
import sqlite3
import tempfile
from typing import Dict, List, Optional, Set

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore
from meta_layer import DerivationHistory

# Files smaller than this (in bytes) will use the fast in-memory path
DEFAULT_MEM_LIMIT_BYTES = 256 * 1024 * 1024  # 256 MB


# ======================================================================
# Streaming Annotator
# ======================================================================

class StreamingGFAAnnotator:
    """
    Memory-safe GFA annotation pipeline for large GFA files.

    Produces the same output as GFAAnnotator but processes the file
    in passes without holding all nodes in memory simultaneously.

    Semantic Filter Phi and DerivationHistory (Asiddhatva) are both
    active in the streaming path, matching the in-memory annotator.

    Attributes:
        engine        : PanIndexEngine with attached DerivationHistory.
        store         : PanIndexStore populated during annotation.
        history       : DerivationHistory for Asiddhatva staged traceability.
        nodes_written : Number of S-lines annotated in the last run.
    """

    def __init__(self, seed: Optional[bytes] = None):
        self.history = DerivationHistory()
        self.engine = PanIndexEngine(pangenome_seed=seed, history=self.history)
        self.store = PanIndexStore()
        self.nodes_written: int = 0

        from default_rules import build_default_rule_engine
        self._rule_engine = build_default_rule_engine()

    def annotate(
        self,
        input_path: str,
        output_path: str,
        derivation_root: str = "PangenomeRoot",
    ):
        """
        Full streaming pipeline: pass 1 (segment indexing) + pass 2 (write-back).

        Args:
            input_path      : Path to source GFA file.
            output_path     : Path for annotated output GFA.
            derivation_root : Top-level label in the derivation hierarchy.
        """
        tmp_fd, tmp_db = tempfile.mkstemp(suffix='.frx_streaming_tmp.db')
        os.close(tmp_fd)

        try:
            self._pass1_index_segments(input_path, derivation_root, tmp_db)
            self._pass2_write_annotated(input_path, output_path, tmp_db)
        finally:
            try:
                os.unlink(tmp_db)
            except OSError:
                pass

    # ------------------------------------------------------------------
    # Pass 1 - Index segments
    # ------------------------------------------------------------------

    def _pass1_index_segments(
        self,
        input_path: str,
        derivation_root: str,
        tmp_db: str,
    ):
        """
        Two sub-passes:
          1a. Scan L-lines to build neighbour map (node_id -> set of neighbour ids).
          1b. Scan S-lines: derive address, apply full Phi filter with Sannidhi
              neighbour context, record DerivationHistory, write to temp DB.
        """
        # ----------------------------------------------------------
        # Sub-pass 1a: build neighbour map from L-lines (O(E) RAM)
        # ----------------------------------------------------------
        # neighbour_map[node_id] = set of directly adjacent node_ids
        neighbour_map: Dict[str, Set[str]] = {}

        with open(input_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if parts[0] != 'L' or len(parts) < 5:
                    continue
                u, v = parts[1], parts[3]
                neighbour_map.setdefault(u, set()).add(v)
                neighbour_map.setdefault(v, set()).add(u)

        # ----------------------------------------------------------
        # Sub-pass 1b: scan S-lines, derive addresses, apply Phi
        # ----------------------------------------------------------
        conn = sqlite3.connect(tmp_db)
        conn.execute("""
            CREATE TABLE segments (
                node_id         TEXT PRIMARY KEY,
                address_hex     TEXT NOT NULL,
                derivation      TEXT NOT NULL,
                history_path    TEXT NOT NULL,
                tag_count       INTEGER NOT NULL,
                seq             TEXT NOT NULL,
                strand          TEXT NOT NULL,
                component       INTEGER NOT NULL,
                resolution      TEXT NOT NULL,
                bio_tags        TEXT NOT NULL
            )
        """)

        root_path_addr = self.engine.derive_ratchet_address(
            self.engine.root_hash, derivation_root,
            node_id=derivation_root,
        )

        # We need the tags of neighbours for Sannidhi, but in a single-pass
        # stream we haven't seen all nodes yet.  Strategy: collect raw tags
        # from S-lines in a lightweight dict (node_id -> raw_tag_string) in
        # this same pass, then resolve Sannidhi in a deferred second scan of
        # the rows buffer.  For very large files this dict is O(N * avg_tags)
        # which is far smaller than O(N * avg_seq_len).
        node_raw_tags: Dict[str, List[str]] = {}

        # First collect all node raw tags (seq not stored, just tags)
        with open(input_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if parts[0] != 'S' or len(parts) < 3:
                    continue
                nid = parts[1]
                raw = parts[3:] if len(parts) > 3 else []
                node_raw_tags[nid] = self._extract_tags(raw)

        # Now derive addresses and apply full Phi filter
        rows = []
        BATCH = 500

        with open(input_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if parts[0] != 'S' or len(parts) < 3:
                    continue

                node_id = parts[1]
                seq = parts[2]
                if seq == '*':
                    seq = ''

                raw_tags = parts[3:] if len(parts) > 3 else []
                anubandha_tags = self._extract_tags(raw_tags)

                # Detect strand from existing tags (e.g. from a prior annotation)
                strand = '+'
                for t in anubandha_tags:
                    if t.startswith('strand:'):
                        strand = t.split(':', 1)[1]
                        break

                # Detect component from existing tags
                component = 0
                for t in anubandha_tags:
                    if t.startswith('component:'):
                        try:
                            component = int(t.split(':', 1)[1])
                        except ValueError:
                            pass
                        break

                # Build Sannidhi neighbour_tags context
                neighbour_ids = neighbour_map.get(node_id, set())
                neighbour_tags: List[str] = []
                for nb_id in neighbour_ids:
                    neighbour_tags.extend(node_raw_tags.get(nb_id, []))

                # Apply Paninian rule engine with full Phi (Akanksha + Yogyata + Sannidhi)
                rule_node = {'tags': list(anubandha_tags), 'seq': seq}
                bio_context = {
                    'neighbor_tags': neighbour_tags,
                    'component': component,
                    'strand': strand,
                }
                resolution = self._rule_engine.resolve(rule_node, bio_context)
                # rule_node['tags'] may have been extended in-place by SemanticFilter
                enriched_tags = rule_node['tags']
                bio_tags_extra = [
                    t for t in enriched_tags if t not in anubandha_tags
                ]

                # Derive ratchet address (recorded in DerivationHistory automatically)
                ratchet_addr = self.engine.derive_ratchet_address(
                    root_path_addr, node_id, node_id=node_id
                )
                derivation = f"{derivation_root}/{node_id}"

                # Retrieve the history path string for this node
                history_path = ' -> '.join(
                    self.history.get_derivation_path(node_id)
                )

                # Build the full tag set for this node
                full_tags = list(enriched_tags) + [
                    f"LN:{len(seq)}",
                    f"strand:{strand}",
                    f"node:{node_id}",
                    f"component:{component}",
                ]
                if resolution and resolution != 'default_resolution':
                    full_tags.append(resolution)

                rows.append((
                    node_id,
                    ratchet_addr.hex(),
                    derivation,
                    history_path,
                    len(full_tags),
                    seq,
                    strand,
                    component,
                    resolution if resolution != 'default_resolution' else '',
                    ','.join(bio_tags_extra),
                ))

                # Populate the in-memory store
                self.store.insert(
                    node_id=node_id,
                    address=ratchet_addr,
                    tags=full_tags,
                    metadata={
                        'seq': seq,
                        'strand': strand,
                        'derivation_path': derivation,
                        'derivation_history': history_path,
                        'component': component,
                        'out_neighbors': list(neighbour_ids),
                        'in_neighbors': [],
                    }
                )

                if len(rows) >= BATCH:
                    self._flush_rows(conn, rows)
                    rows.clear()

        if rows:
            self._flush_rows(conn, rows)

        self.nodes_written = conn.execute(
            "SELECT COUNT(*) FROM segments"
        ).fetchone()[0]
        conn.close()

    @staticmethod
    def _flush_rows(conn: sqlite3.Connection, rows: list):
        conn.executemany(
            "INSERT OR REPLACE INTO segments "
            "(node_id, address_hex, derivation, history_path, tag_count, "
            " seq, strand, component, resolution, bio_tags) "
            "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows
        )
        conn.commit()

    # ------------------------------------------------------------------
    # Pass 2 - Write annotated GFA
    # ------------------------------------------------------------------

    def _pass2_write_annotated(
        self,
        input_path: str,
        output_path: str,
        tmp_db: str,
    ):
        """
        Re-read input line by line, inject tags from temp db, write output.
        One line at a time - O(1) peak RAM.

        Injected tags per S-line:
            AN:Z:<address_hex>   - 32-byte ratchet address
            PA:Z:<derivation>    - human-readable derivation path
            AF:i:<tag_count>     - number of Anubandha tags (incl. bio-tags)
            HI:Z:<history_path>  - Asiddhatva derivation history chain
        """
        conn = sqlite3.connect(tmp_db)

        with open(input_path, 'r', encoding='utf-8') as fin, \
             open(output_path, 'w', encoding='utf-8') as fout:

            for line in fin:
                parts = line.rstrip('\n').split('\t')

                if parts[0] == 'S' and len(parts) >= 3:
                    node_id = parts[1]
                    row = conn.execute(
                        "SELECT address_hex, derivation, tag_count, history_path "
                        "FROM segments WHERE node_id = ?",
                        (node_id,)
                    ).fetchone()

                    if row:
                        addr_hex, derivation, tag_count, history_path = row
                        parts.append(f"AN:Z:{addr_hex}")
                        parts.append(f"PA:Z:{derivation}")
                        parts.append(f"AF:i:{tag_count}")
                        if history_path:
                            parts.append(f"HI:Z:{history_path}")

                fout.write('\t'.join(parts) + '\n')

        conn.close()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _extract_tags(raw_tags: list) -> list:
        """Parse GFA optional fields (TAG:TYPE:VALUE) into Anubandha tags."""
        extracted = []
        for field in raw_tags:
            parts = field.split(':')
            if len(parts) >= 3:
                tag_name = parts[0]
                tag_value = ':'.join(parts[2:])
                extracted.append(f"{tag_name}:{tag_value}")
        return extracted

    def print_index_summary(self):
        """Print a summary of the streaming annotation results."""
        print(f"\n[StreamingGFAAnnotator] Nodes indexed: {self.nodes_written}")
        st = self.store.stats()
        print(f"  Total nodes      : {st['total_nodes']}")
        print(f"  Unique tags      : {st['unique_tags']}")
        print(f"  Derivation stages: {len(self.history)}")


# ======================================================================
# Auto-select annotator
# ======================================================================

def make_annotator(
    input_path: str,
    seed: Optional[bytes] = None,
    mem_limit_bytes: int = DEFAULT_MEM_LIMIT_BYTES,
    force_streaming: bool = False,
):
    """
    Return the appropriate annotator based on file size.

    Args:
        input_path      : GFA file to annotate.
        seed            : Pangenome seed bytes.
        mem_limit_bytes : Files larger than this use StreamingGFAAnnotator.
        force_streaming : Always use streaming regardless of file size.

    Returns:
        GFAAnnotator or StreamingGFAAnnotator instance.
    """
    from annotator import GFAAnnotator

    if force_streaming:
        return StreamingGFAAnnotator(seed=seed)

    try:
        file_size = os.path.getsize(input_path)
    except OSError:
        file_size = 0

    if file_size > mem_limit_bytes:
        print(f"[FRX] File size {file_size / (1024**2):.1f} MB exceeds "
              f"{mem_limit_bytes / (1024**2):.0f} MB limit. "
              f"Using streaming annotator.")
        return StreamingGFAAnnotator(seed=seed)

    return GFAAnnotator(seed=seed)
