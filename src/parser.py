import sys
import os
import logging
from typing import Dict, List, Set
from engine import PanIndexEngine

logger = logging.getLogger(__name__)


class GFAParser:
    """
    Normalization Layer: GFA to PanIndex Internal Graph.

    Handles both DAG and cyclic (plasmid/inversion) topologies.
    Cyclic nodes are resolved via lexicographic canonical anchor election
    (same strategy as BipartiteGraph) rather than being silently dropped.
    """

    def __init__(self, engine: PanIndexEngine):
        self.engine = engine
        self.nodes: Dict[str, dict] = {}
        self.edges: List[tuple] = []

    def parse_file(self, filepath: str):
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if not parts:
                    continue
                record_type = parts[0]
                if record_type == 'S' and len(parts) >= 3:
                    node_id = parts[1]
                    seq = parts[2] if parts[2] != '*' else ''
                    tags = parts[3:] if len(parts) > 3 else []
                    self.nodes[node_id] = {
                        'seq': seq,
                        'tags': tags,
                        'out_neighbors': [],
                        'in_neighbors': [],
                        'panindex_addr': None,
                    }
                elif record_type == 'L' and len(parts) >= 4:
                    u, v = parts[1], parts[3]
                    self.edges.append((u, v))

        for u, v in self.edges:
            if u in self.nodes and v in self.nodes:
                self.nodes[u]['out_neighbors'].append(v)
                self.nodes[v]['in_neighbors'].append(u)

    def compute_all_addresses(self, root_context: str = "Root"):
        """
        Compute PanIndex addresses for all nodes.

        DAG nodes are processed in topological order.
        Cyclic nodes are resolved by electing the lexicographically smallest
        sequence in each cycle as the canonical anchor (Local Root), then
        deriving the remaining cycle members from it — identical to the
        strategy used in BipartiteGraph.compute_addresses().
        """
        dag_nodes, cycle_groups = self._classify_nodes()

        # Build a forward adjacency map for efficient child lookup
        children: Dict[str, List[str]] = {nid: [] for nid in self.nodes}
        for nid, node in self.nodes.items():
            for child in node['out_neighbors']:
                children[nid].append(child)

        # Derive root path address (mirrors annotator / query_by_path)
        root_path_addr = self.engine.derive_ratchet_address(
            self.engine.root_hash, root_context
        )

        # --- Process DAG nodes in topological order ---
        for node_id in dag_nodes:
            self._derive_address(node_id, root_path_addr)

        # --- Process each cycle group ---
        for cycle_group in cycle_groups:
            # Elect canonical anchor: lexicographically smallest sequence
            anchor = min(cycle_group, key=lambda nid: self.nodes[nid]['seq'])
            logger.debug(
                "Cycle detected (%d nodes). Canonical anchor: '%s'",
                len(cycle_group), anchor,
            )

            # Anchor derives from root
            self._derive_address(anchor, root_path_addr)

            # BFS from anchor to address remaining cycle members
            visited: Set[str] = {anchor}
            queue = list(self.nodes[anchor]['out_neighbors'])
            while queue:
                nid = queue.pop(0)
                if nid not in cycle_group or nid in visited:
                    continue
                visited.add(nid)
                # Use anchor's address as parent for all cycle members
                # (stable: anchor address is fixed regardless of traversal order)
                anchor_addr = self.nodes[anchor]['panindex_addr']
                derived = self.engine.derive_ratchet_address(anchor_addr, nid)
                self.nodes[nid]['panindex_addr'] = self.engine.compute_node_address(
                    self.nodes[nid]['seq'], derived
                )
                queue.extend(self.nodes[nid]['out_neighbors'])

        unresolved = [
            nid for nid, n in self.nodes.items() if n['panindex_addr'] is None
        ]
        if unresolved:
            logger.warning(
                "%d node(s) could not be addressed (isolated or malformed): %s",
                len(unresolved), unresolved[:10],
            )
            # Fallback: derive from root so they are never None
            for nid in unresolved:
                self._derive_address(nid, root_path_addr)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _derive_address(self, node_id: str, root_path_addr: bytes):
        """Derive and store the PanIndex address for a single node."""
        node = self.nodes[node_id]
        if node['panindex_addr'] is not None:
            return

        if node['in_neighbors']:
            resolved_parents = [
                self.nodes[p]['panindex_addr']
                for p in node['in_neighbors']
                if self.nodes[p]['panindex_addr'] is not None
            ]
            parent_addr = (
                self.engine.compute_commutative_hash(resolved_parents)
                if resolved_parents
                else root_path_addr
            )
        else:
            parent_addr = root_path_addr

        derived = self.engine.derive_ratchet_address(parent_addr, node_id)
        node['panindex_addr'] = self.engine.compute_node_address(
            node['seq'], derived
        )

    def _classify_nodes(self):
        """
        Separate nodes into DAG nodes (topological order) and cycle groups.

        Returns:
            (dag_order: List[str], cycle_groups: List[List[str]])
        """
        in_degree = {
            nid: len(node['in_neighbors'])
            for nid, node in self.nodes.items()
        }

        # Kahn's algorithm
        queue = [nid for nid, deg in in_degree.items() if deg == 0]
        dag_order: List[str] = []

        while queue:
            u = queue.pop(0)
            dag_order.append(u)
            for v in self.nodes[u]['out_neighbors']:
                in_degree[v] -= 1
                if in_degree[v] == 0:
                    queue.append(v)

        dag_set = set(dag_order)
        cycle_nodes = [nid for nid in self.nodes if nid not in dag_set]

        if not cycle_nodes:
            return dag_order, []

        # Group cycle nodes into connected components via undirected adjacency
        cycle_set = set(cycle_nodes)
        visited: Set[str] = set()
        cycle_groups: List[List[str]] = []

        for start in cycle_nodes:
            if start in visited:
                continue
            component: List[str] = []
            stack = [start]
            while stack:
                nid = stack.pop()
                if nid in visited or nid not in cycle_set:
                    continue
                visited.add(nid)
                component.append(nid)
                for nb in self.nodes[nid]['out_neighbors'] + self.nodes[nid]['in_neighbors']:
                    if nb not in visited and nb in cycle_set:
                        stack.append(nb)
            if component:
                cycle_groups.append(component)

        return dag_order, cycle_groups

if __name__ == "__main__":
    import sys as _sys
    logging.basicConfig(level=logging.DEBUG)
    _engine = PanIndexEngine()
    _parser = GFAParser(_engine)
    _parser.parse_file("test.gfa")
    _parser.compute_all_addresses()
    for nid, data in _parser.nodes.items():
        addr_hex = data['panindex_addr'].hex()
        print(f"Node {nid}: {addr_hex[:16]}...")
