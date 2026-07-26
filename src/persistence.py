"""
FRX Index Persistence - SQLite-backed save and load for PanIndexStore.

Uses Python's built-in sqlite3 module (zero extra dependencies).

Schema
------
nodes table   : one row per indexed node (node_id, address_hex, metadata_json)
node_tags table: one row per (tag, node_id) pair

Indexes
-------
idx_address : on nodes.address_hex  -> fast address lookup
idx_tag     : on node_tags.tag      -> fast tag lookup

Round-trip guarantee
--------------------
save_store(store, path) -> load_store(path) produces a PanIndexStore
that returns identical results for all three query modes as the original.

Crash safety
------------
save_store writes to a temp file in the same directory, then atomically
renames it over the destination. A crash mid-write never corrupts an
existing database. WAL journal mode ensures durability on the temp file.
"""

import json
import os
import sqlite3
import tempfile
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from index import PanIndexStore


# ======================================================================
# DDL
# ======================================================================

_DDL = """
CREATE TABLE IF NOT EXISTS nodes (
    node_id       TEXT PRIMARY KEY,
    address_hex   TEXT NOT NULL,
    merkle_hex    TEXT NOT NULL DEFAULT '',
    content_hex   TEXT NOT NULL DEFAULT '',
    metadata_json TEXT NOT NULL,
    organism      TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS node_tags (
    tag     TEXT NOT NULL,
    node_id TEXT NOT NULL,
    PRIMARY KEY (tag, node_id)
);

CREATE INDEX IF NOT EXISTS idx_address  ON nodes(address_hex);
CREATE INDEX IF NOT EXISTS idx_merkle   ON nodes(merkle_hex);
CREATE INDEX IF NOT EXISTS idx_tag      ON node_tags(tag);
CREATE INDEX IF NOT EXISTS idx_organism ON nodes(organism);
"""


# ======================================================================
# Save
# ======================================================================

def save_store(store: "PanIndexStore", path: str) -> int:
    """
    Serialize a PanIndexStore to a SQLite database file.

    Writes atomically: data goes to a sibling temp file first, then
    os.replace() swaps it into place. A crash mid-write never leaves
    the destination in a corrupt state.

    Args:
        store : Populated PanIndexStore instance.
        path  : Filesystem path for the output .db file.
                Created or overwritten atomically.

    Returns:
        Number of nodes written.
    """
    dir_name = os.path.dirname(os.path.abspath(path)) or '.'
    fd, tmp_path = tempfile.mkstemp(dir=dir_name, suffix='.frx_tmp')
    os.close(fd)

    try:
        conn = sqlite3.connect(tmp_path)
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.executescript(_DDL)

        node_rows = []
        tag_rows = []

        for node_id, record in store._node_store.items():
            if record.get('metadata', {}).get('bipartite_state'):
                continue

            node_rows.append((
                node_id,
                record['address'],
                record.get('merkle_addr', ''),
                record.get('content_id', ''),
                json.dumps(record['metadata'], ensure_ascii=True),
                record.get('organism', ''),
            ))
            for tag in record['tags']:
                tag_rows.append((tag, node_id))

        with conn:
            conn.executemany(
                "INSERT INTO nodes "
                "(node_id, address_hex, merkle_hex, content_hex, metadata_json, organism) "
                "VALUES (?, ?, ?, ?, ?, ?)",
                node_rows
            )
            conn.executemany(
                "INSERT OR IGNORE INTO node_tags (tag, node_id) VALUES (?, ?)",
                tag_rows
            )

        conn.close()
        os.replace(tmp_path, path)
        return len(node_rows)

    except Exception:
        try:
            os.unlink(tmp_path)
        except OSError:
            pass
        raise


# ======================================================================
# Load
# ======================================================================

def load_store(path: str) -> "PanIndexStore":
    """
    Deserialize a PanIndexStore from a SQLite database file.

    Rebuilds the in-memory bisect address index and tag dict from the
    stored rows so all existing query paths work identically to a
    freshly built store.

    Args:
        path: Filesystem path of a .db file written by save_store().

    Returns:
        Populated PanIndexStore instance.

    Raises:
        FileNotFoundError if path does not exist.
        sqlite3.DatabaseError if the file is not a valid FRX index.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"FRX index not found: {path}")

    from index import PanIndexStore

    conn = sqlite3.connect(path)
    store = PanIndexStore()

    try:
        # Handle old .db files that may lack content_hex, merkle_hex, or organism columns
        try:
            cursor = conn.execute(
                "SELECT node_id, address_hex, merkle_hex, content_hex, metadata_json, organism "
                "FROM nodes ORDER BY address_hex"
            )
            node_records = [(r[0], r[1], r[2], r[3], r[4], r[5]) for r in cursor.fetchall()]
        except sqlite3.OperationalError:
            try:
                cursor = conn.execute(
                    "SELECT node_id, address_hex, metadata_json, organism "
                    "FROM nodes ORDER BY address_hex"
                )
                node_records = [(r[0], r[1], '', '', r[2], r[3]) for r in cursor.fetchall()]
            except sqlite3.OperationalError:
                cursor = conn.execute(
                    "SELECT node_id, address_hex, metadata_json "
                    "FROM nodes ORDER BY address_hex"
                )
                node_records = [(r[0], r[1], '', '', r[2], '') for r in cursor.fetchall()]

        tag_cursor = conn.execute("SELECT node_id, tag FROM node_tags")
        tags_by_node = {}
        for node_id, tag in tag_cursor.fetchall():
            tags_by_node.setdefault(node_id, []).append(tag)

        for node_id, address_hex, merkle_hex, content_hex, metadata_json, organism in node_records:
            tags = tags_by_node.get(node_id, [])
            metadata = json.loads(metadata_json)
            address_bytes = bytes.fromhex(address_hex)
            content_hex = content_hex or metadata.get('content_id', '')

            # Strip is_alias from loaded metadata: a node restored from disk
            # is a first-class queryable entry, not a bipartite alias.
            # Keeping is_alias=True would cause PanIndexStore.insert() to skip
            # tag indexing, leaving lookup_by_tag() returning empty results.
            metadata.pop('is_alias', None)

            store.insert(
                node_id, address_bytes, tags, metadata,
                organism=organism,
                merkle_addr=merkle_hex or '',
                content_id=content_hex,
            )

    finally:
        conn.close()

    return store


# ======================================================================
# Stats
# ======================================================================

def db_stats(path: str) -> dict:
    """
    Return lightweight statistics about a saved FRX index without
    loading the full store into memory.

    Returns:
        Dict with keys: nodes, tags, unique_tags, db_size_bytes.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"FRX index not found: {path}")

    conn = sqlite3.connect(path)
    try:
        (node_count,) = conn.execute("SELECT COUNT(*) FROM nodes").fetchone()
        (tag_count,)  = conn.execute("SELECT COUNT(*) FROM node_tags").fetchone()
        (unique_tags,) = conn.execute(
            "SELECT COUNT(DISTINCT tag) FROM node_tags"
        ).fetchone()
    finally:
        conn.close()

    return {
        'nodes': node_count,
        'tags': tag_count,
        'unique_tags': unique_tags,
        'db_size_bytes': os.path.getsize(path),
    }
