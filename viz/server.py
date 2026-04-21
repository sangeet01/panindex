import sys
import os
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

# Initialize with real data
engine = PanIndexEngine(pangenome_seed=b"panindex_viz_real_seed_v1")
store = PanIndexStore()
graph = BipartiteGraph(engine, store)

# Use the real converted Kaggle GFA
input_gfa = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data', 'kp_pangenome.gfa'))

print(f"Loading real pangenome and building BipartiteGraph: {input_gfa}...")
graph.parse_gfa(input_gfa)
graph.compute_addresses()

query_engine = PanIndexQuery(engine, store)

@app.route('/graph', methods=['GET'])
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

if __name__ == '__main__':
    # Run server on port 5000
    app.run(debug=True, port=5000)
