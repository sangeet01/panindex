import sys
import os
from engine import PanIndexEngine

class GFAParser:
    """
    Normalization Layer: GFA to PanIndex Internal Graph.
    """
    def __init__(self, engine: PanIndexEngine):
        self.engine = engine
        self.nodes = {} # id -> {seq, tags, neighbors, panindex_addr}
        self.edges = []

    def parse_file(self, filepath: str):
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if not parts: continue
                
                type = parts[0]
                if type == 'S': # Segment
                    node_id = parts[1]
                    seq = parts[2]
                    tags = parts[3:] if len(parts) > 3 else []
                    self.nodes[node_id] = {
                        'seq': seq,
                        'tags': tags,
                        'out_neighbors': [],
                        'in_neighbors': [],
                        'panindex_addr': None
                    }
                elif type == 'L': # Link
                    u, v = parts[1], parts[3]
                    self.edges.append((u, v))

        # Build adjacency
        for u, v in self.edges:
            if u in self.nodes and v in self.nodes:
                self.nodes[u]['out_neighbors'].append(v)
                self.nodes[v]['in_neighbors'].append(u)

    def compute_all_addresses(self, root_context: str = "Root"):
        """
        Compute PanIndex addresses for all nodes in the topological layer.
        """
        # For a true Merkle graph, we'd need a topological sort (if DAG).
        # For a pangenome, we can use the neighborhood context.
        # Step 1: Base sequence hashes
        # Step 2: Neighbor context (Commutative XOR)
        
        # We'll use a multi-pass approach or topological sort for a DAG
        sorted_nodes = self._topological_sort()
        
        for node_id in sorted_nodes:
            node = self.nodes[node_id]
            
            # Use Ratchet for hierarchical context (simplified here as parent -> child)
            # In a real system, we'd use the hierarchical tags (Chr, Gene)
            parent_addr = self.engine.root_hash # Default to root if no parents
            
            if node['in_neighbors']:
                # Commutative XOR of parent addresses
                parent_hashes = [self.nodes[p]['panindex_addr'] for p in node['in_neighbors'] if self.nodes[p]['panindex_addr']]
                if parent_hashes:
                    parent_addr = self.engine.compute_commutative_hash(parent_hashes)
            
            # Derivation Path component (Simplified structural label)
            # In production, this would be derived from metadata tags (AN:Z:...)
            structural_ctx = node_id 
            derived_addr = self.engine.derive_ratchet_address(parent_addr, structural_ctx)
            
            # Content component
            node['panindex_addr'] = self.engine.compute_node_address(node['seq'], derived_addr)

    def _topological_sort(self):
        # Kahn's algorithm
        in_degree = {node_id: len(node['in_neighbors']) for node_id, node in self.nodes.items()}
        queue = [node_id for node_id, degree in in_degree.items() if degree == 0]
        sorted_nodes = []
        
        while queue:
            u = queue.pop(0)
            sorted_nodes.append(u)
            for v in self.nodes[u]['out_neighbors']:
                in_degree[v] -= 1
                if in_degree[v] == 0:
                    queue.append(v)
        return sorted_nodes

if __name__ == "__main__":
    engine = PanIndexEngine()
    parser = GFAParser(engine)
    parser.parse_file("test.gfa")
    parser.compute_all_addresses()
    
    for nid, data in parser.nodes.items():
        addr_hex = data['panindex_addr'].hex()
        print(f"Node {nid}: {addr_hex[:16]}...")
