# PanIndex CLI & API Reference

## Command Line Interface (`frx`)

The `frx` tool provides the main interface for interacting with the FractalIndex.

### `frx vg-import`
Imports a `vg`-generated GFA file, normalizes it, and builds a `.frx.db` database.
- `--gfa`: Input GFA file path (required)
- `--out`: Output SQLite db path (required)
- `--build-kmer-index`: Extracts k-mers for `frx pattern` compatibility.
- `--streaming`: Uses memory-efficient streaming parsing for gigabyte-scale GFAs.

### `frx pattern`
Finds nodes containing an exact sequence using the k-mer index.
- `--index`: The database path
- `--pattern`: DNA sequence string

### `frx variant`
Compares two physical nodes using MinHash LSH to detect sequence homology.
- `--index`: The database path
- `--node1`, `--node2`: Physical node IDs to compare
- `--threshold`: Similarity threshold (default: 0.5)

### `frx query`
Address calculation and structural lookup based on tags or structural paths.
- `--path`: Derivation path (`PangenomeRoot/Chromosome_1/Node_14`)
- `--tag`: Look for anubandha tag (e.g., `AMR:blaTEM`)

### `frx hgt-sim`
Runs an interactive, built-in demonstration of the Horizontal Gene Transfer algorithms.

---

## REST API Reference

The web visualization relies on a Python Flask backend (`viz.server`). 

**Authentication**: 
Configure the environment variable `FRX_API_TOKEN`. Requests must include:
`Authorization: Bearer <TOKEN>`
If `FRX_API_TOKEN` is unset, the server runs in DEV Mode (no auth).

### `GET /health`
Returns the status of the server. No authentication required.
```json
{
  "status": "ok",
  "nodes": 451,
  "auth": true
}
```

### `GET /graph`
Returns the full unrolled Bipartite Graph formatted for Cytoscape.js (`.elements`).
Requires Auth. Rate-limited to 100/min per IP.

### `GET /query?path=...` (or `?tag=...` or `?similarity=...`)
Executes an engine search against the active pangenome and returns matching physical nodes along with their structural hashes.
Requires Auth. Rate-limited.
```json
{
  "matched_nodes": ["state_14_fwd", "state_14_rev"],
  "derived_address": "8a32b0c9f13..."
}
```
