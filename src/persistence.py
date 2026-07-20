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
"""

import json
import sqlite3
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
    metadata_json TEXT NOT NULL,
    organism      TEXT NOT NULL DEFAULT ''
);

CREATE TABLE IF NOT EXISTS node_tags (
    tag     TEXT NOT NULL,
    node_id TEXT NOT NULL,
    PRIMARY KEY (tag, node_id)
);

CREATE INDEX IF NOT EXISTS idx_address  ON nodes(address_hex);
CREATE INDEX IF NOT EXISTS idx_tag      ON node_tags(tag);
CREATE INDEX IF NOT EXISTS idx_organism ON nodes(organism);
"""


# ======================================================================
# Save
# ======================================================================

def save_store(store: "PanIndexStore", path: str) -> int:
    """
    Serialize a PanIndexStore to a SQLite database file.

    Args:
        store : Populated PanIndexStore instance.
        path  : Filesystem path for the output .db file.
                Will be created or overwritten.

    Returns:
        Number of nodes written.
    """
    conn = sqlite3.connect(path)
    try:
        conn.executescript(_DDL)

        # Clear existing data (overwrite semantics)
        conn.execute("DELETE FROM node_tags")
        conn.execute("DELETE FROM nodes")

        node_rows = []
        tag_rows = []

        for node_id, record in store._node_store.items():
            node_rows.append((
                node_id,
                record['address'],
                json.dumps(record['metadata'], ensure_ascii=True),
                record.get('organism', ''),
            ))
            for tag in record['tags']:
                tag_rows.append((tag, node_id))

        conn.executemany(
            "INSERT INTO nodes (node_id, address_hex, metadata_json, organism) "
            "VALUES (?, ?, ?, ?)",
            node_rows
        )
        conn.executemany(
            "INSERT OR IGNORE INTO node_tags (tag, node_id) VALUES (?, ?)",
            tag_rows
        )

        conn.commit()
        return len(node_rows)

    finally:
        conn.close()


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
    import os
    if not os.path.isfile(path):
        raise FileNotFoundError(f"FRX index not found: {path}")

    # Import here to avoid circular dependency at module load time
    from index import PanIndexStore

    conn = sqlite3.connect(path)
    store = PanIndexStore()

    try:
        # Load all nodes - handle old .db files that lack the organism column
        try:
            cursor = conn.execute(
                "SELECT node_id, address_hex, metadata_json, organism "
                "FROM nodes ORDER BY address_hex"
            )
            node_records = [(r[0], r[1], r[2], r[3]) for r in cursor.fetchall()]
        except sqlite3.OperationalError:
            # Old .db without organism column - default to empty string
            cursor = conn.execute(
                "SELECT node_id, address_hex, metadata_json FROM nodes ORDER BY address_hex"
            )
            node_records = [(r[0], r[1], r[2], '') for r in cursor.fetchall()]

        # Load all tags keyed by node_id
        tag_cursor = conn.execute(
            "SELECT node_id, tag FROM node_tags"
        )
        tags_by_node = {}
        for node_id, tag in tag_cursor.fetchall():
            tags_by_node.setdefault(node_id, []).append(tag)

        for node_id, address_hex, metadata_json, organism in node_records:
            tags = tags_by_node.get(node_id, [])
            metadata = json.loads(metadata_json)
            address_bytes = bytes.fromhex(address_hex)
            store.insert(node_id, address_bytes, tags, metadata, organism=organism)

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
    import os
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
