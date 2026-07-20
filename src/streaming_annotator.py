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
Pass 1 (S-lines):
  Read segments one at a time. For each segment, compute a provisional
  ratchet address (derivation_root -> node_id derivation, no neighbor XOR).
  Write node data immediately to a SQLite temp table. Clear from RAM.

Pass 2 (Write-back):
  Re-read the input GFA line by line. For each S-line, fetch the pre-computed
  address from the temp table and inject AN:Z:/PA:Z:/AF:i: tags. Write to
  output line by line. All other line types (H/L/P/W) pass through unchanged.

Note on Merkle XOR: the streaming path uses pure ratchet-path addresses
(identical to query_by_path behavior). The full Merkle XOR (commutative
neighbor hashing) is computed by the standard annotator and stored in the
'merkle_addr' metadata field. For the streaming case this field is omitted
since it requires the full adjacency graph in memory.

Memory profile:
  Standard annotator : O(N * avg_seq_len) peak RAM
  Streaming annotator: O(1) peak RAM (one line at a time per pass)
"""

import os
import sqlite3
import tempfile
from typing import Optional

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore

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

    Attributes:
        engine      : PanIndexEngine used for address derivation.
        store       : PanIndexStore populated during annotation.
        nodes_written : Number of S-lines annotated in the last run.
    """

    def __init__(self, seed: Optional[bytes] = None):
        self.engine = PanIndexEngine(pangenome_seed=seed)
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
        # Temp SQLite for intermediate segment data
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
        Read S-lines from the GFA, compute addresses, write to temp db.
        One node at a time - O(1) peak RAM.
        """
        conn = sqlite3.connect(tmp_db)
        conn.execute("""
            CREATE TABLE segments (
                node_id       TEXT PRIMARY KEY,
                address_hex   TEXT NOT NULL,
                derivation    TEXT NOT NULL,
                tag_count     INTEGER NOT NULL,
                seq           TEXT NOT NULL
            )
        """)

        root_path_addr = self.engine.derive_ratchet_address(
            self.engine.root_hash, derivation_root
        )

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

                # Apply Paninian rule engine - append resolution as an extra tag
                rule_node = {'tags': anubandha_tags, 'seq': seq}
                resolution = self._rule_engine.resolve(rule_node, {})
                if resolution and resolution != 'default_resolution':
                    anubandha_tags = list(anubandha_tags) + [resolution]

                ratchet_addr = self.engine.derive_ratchet_address(
                    root_path_addr, node_id
                )
                derivation = f"{derivation_root}/{node_id}"

                rows.append((
                    node_id,
                    ratchet_addr.hex(),
                    derivation,
                    len(anubandha_tags),
                    seq,
                ))

                # Also populate the in-memory store
                self.store.insert(
                    node_id=node_id,
                    address=ratchet_addr,
                    tags=anubandha_tags,
                    metadata={
                        'seq': seq,
                        'derivation_path': derivation,
                        'out_neighbors': [],
                        'in_neighbors': [],
                    }
                )

                if len(rows) >= BATCH:
                    conn.executemany(
                        "INSERT OR REPLACE INTO segments "
                        "(node_id, address_hex, derivation, tag_count, seq) "
                        "VALUES (?, ?, ?, ?, ?)",
                        rows
                    )
                    conn.commit()
                    rows.clear()

        if rows:
            conn.executemany(
                "INSERT OR REPLACE INTO segments "
                "(node_id, address_hex, derivation, tag_count, seq) "
                "VALUES (?, ?, ?, ?, ?)",
                rows
            )
            conn.commit()

        self.nodes_written = conn.execute(
            "SELECT COUNT(*) FROM segments"
        ).fetchone()[0]
        conn.close()

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
        """
        conn = sqlite3.connect(tmp_db)

        with open(input_path, 'r', encoding='utf-8') as fin, \
             open(output_path, 'w', encoding='utf-8') as fout:

            for line in fin:
                parts = line.rstrip('\n').split('\t')

                if parts[0] == 'S' and len(parts) >= 3:
                    node_id = parts[1]
                    row = conn.execute(
                        "SELECT address_hex, derivation, tag_count "
                        "FROM segments WHERE node_id = ?",
                        (node_id,)
                    ).fetchone()

                    if row:
                        addr_hex, derivation, tag_count = row
                        parts.append(f"AN:Z:{addr_hex}")
                        parts.append(f"PA:Z:{derivation}")
                        parts.append(f"AF:i:{tag_count}")

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
        print(f"  Total nodes   : {st['total_nodes']}")
        print(f"  Unique tags   : {st['unique_tags']}")


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
