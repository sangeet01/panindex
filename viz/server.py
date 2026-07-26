import sys
import os
import time
from collections import defaultdict
from functools import wraps

from flask import Flask, jsonify, request
from flask_cors import CORS

# Add src to path to import PanIndex components
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

from engine import PanIndexEngine
from index import PanIndexStore
from query import PanIndexQuery
from bipartite import BipartiteGraph

app = Flask(__name__)
CORS(app) # Enable CORS for frontend communication

# ----------------------------------------------------------------------
# Auth & rate-limit configuration
#
# _API_TOKEN: when non-empty, /graph and /query require
#             "Authorization: Bearer <token>". Empty string = dev mode
#             (no auth), matching the "Ensure you connect using your API
#             token if security is enabled" note in the README.
# _RATE_LIMIT / _RATE_WINDOW_SECONDS: simple fixed-window limiter applied
#             per client IP to /graph and /query. /health is exempt.
#
# These are module-level (not captured in closures) so tests can flip
# them per-case via the imported module object.
# ----------------------------------------------------------------------
_API_TOKEN = os.environ.get("FRX_API_TOKEN", "")
_RATE_LIMIT = int(os.environ.get("FRX_RATE_LIMIT", "100"))
_RATE_WINDOW_SECONDS = 60
_rate_registry = defaultdict(list)


def _is_authorized() -> bool:
    if not _API_TOKEN:
        return True  # dev mode: auth disabled
    auth_header = request.headers.get("Authorization", "")
    if not auth_header.startswith("Bearer "):
        return False
    token = auth_header[len("Bearer "):].strip()
    if not token:
        return False
    return token == _API_TOKEN


def _rate_limited(key: str) -> bool:
    """Returns True if `key` has exceeded the current rate limit."""
    now = time.time()
    window_start = now - _RATE_WINDOW_SECONDS
    timestamps = _rate_registry[key]
    while timestamps and timestamps[0] < window_start:
        timestamps.pop(0)
    if len(timestamps) >= _RATE_LIMIT:
        return True
    timestamps.append(now)
    return False


def require_auth_and_rate_limit(fn):
    """Decorator applied to protected endpoints (/graph, /query)."""
    @wraps(fn)
    def wrapper(*args, **kwargs):
        if not _is_authorized():
            return jsonify({"error": "Unauthorized"}), 401

        client_key = request.remote_addr or "unknown"
        if _rate_limited(client_key):
            resp = jsonify({"error": "TooManyRequests"})
            resp.headers["Retry-After"] = str(_RATE_WINDOW_SECONDS)
            return resp, 429

        return fn(*args, **kwargs)
    return wrapper


@app.route('/health', methods=['GET'])
def health():
    """Always public, regardless of auth configuration."""
    return jsonify({"status": "ok", "auth": bool(_API_TOKEN)}), 200


# Initialize with real data
engine = PanIndexEngine(pangenome_seed=b"panindex_viz_real_seed_v1")
store = PanIndexStore()
graph = BipartiteGraph(engine, store)

# Use the real converted Kaggle GFA if present; otherwise fall back to the
# small bundled sample so `python -m viz.server` works out of the box.
input_gfa = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'kp_pangenome.gfa'))
sample_gfa = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'sample_kp_pangenome.gfa'))

if os.path.exists(input_gfa):
    active_gfa = input_gfa
elif os.path.exists(sample_gfa):
    print(f"Note: {input_gfa} not found; falling back to bundled sample dataset.")
    active_gfa = sample_gfa
else:
    active_gfa = None

if active_gfa:
    print(f"Loading pangenome and building BipartiteGraph: {active_gfa}...")
    graph.parse_gfa(active_gfa)
    graph.compute_addresses()
else:
    print("Warning: no pangenome GFA found under data/. Serving an empty graph.")

query_engine = PanIndexQuery(engine, store)

@app.route('/graph', methods=['GET'])
@require_auth_and_rate_limit
def get_graph():
    """
    Exposes the expanded Bipartite graph for Cytoscape.js
    """
    elements = []
    
    # Add Bipartite State Nodes
    for state_key, state in graph.states.items():
        node_info = store.get_node(state_key)
        if not node_info: continue
        
        elements.append({
            'data': {
                'id': state_key,
                'label': f"{state.node_id}{state.strand}",
                'physical_id': state.node_id,
                'strand': state.strand,
                'address': node_info['address'],
                'component': node_info['metadata'].get('component', 0),
                'derivation_path': node_info['metadata'].get('derivation_path', ''),
                'seq_preview': node_info['metadata'].get('seq', '')[:30] + "...",
                'type': 'state'
            }
        })
    
    # Add Bipartite Edges (strictly directed)
    for state_key, node in graph.states.items():
        for parent in node.parents:
            elements.append({
                'data': {
                    'id': f"edge_{parent}_{state_key}",
                    'source': parent,
                    'target': state_key,
                    'type': 'bipartite_link'
                }
            })
            
    return jsonify(elements)

@app.route('/query', methods=['GET'])
@require_auth_and_rate_limit
def run_query():
    path = request.args.get('path')
    tag = request.args.get('tag')
    similarity = request.args.get('similarity')
    
    result = None
    if path:
        result = query_engine.query_by_path(path)
    elif tag:
        result = query_engine.query_by_tag(tag)
    elif similarity:
        result = query_engine.query_by_similarity(similarity)
        
    if result:
        return jsonify({
            'matched_nodes': result.matched_nodes,
            'derived_address': result.derived_address.hex() if result.derived_address else None
        })
    return jsonify({'error': 'No query parameters provided'}), 400

def main():
    """Entry point for the `frx-api` console script."""
    debug = os.environ.get("FRX_DEBUG", "0") == "1"
    port = int(os.environ.get("FRX_PORT", "5000"))
    if not _API_TOKEN:
        print("Warning: FRX_API_TOKEN is not set — running in dev mode with auth disabled.")
    app.run(debug=debug, port=port)


if __name__ == '__main__':
    main()
