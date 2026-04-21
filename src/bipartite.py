"""
Bipartite Expansion Engine for Bidirected Variation Graphs.

Solves the core problem: a standard GFA encodes bidirected edges
(L 1 + 2 - 0M) which introduce cycles that break topological sort
and destroy the Fractal Ratchet's O(1) derivation.

Solution (from vga_circular.txt):
    For every physical node N, create two mathematical states:
        N+ : forward state, hashes raw sequence
        N- : reverse state, hashes reverse complement

    Links are then strictly directed between states, eliminating
    bidirectionality. Remaining cycles (plasmids, inversions) are
    broken by Lexicographical Symmetry Breaking - the node with the
    alphabetically lowest sequence becomes the Local Root.

    The ratchet then derives addresses on this expanded, strictly
    directed graph with no modification.
"""
import hashlib
from typing import Dict, List, Optional, Tuple
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore


COMPLEMENT = str.maketrans('ATGCatgc', 'TACGtacg')


def reverse_complement(seq: str) -> str:
    """Biological reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


class BipartiteNode:
    """
    One state of a physical GFA segment.
    key = "node_id+" or "node_id-"
    """
    __slots__ = ('key', 'node_id', 'strand', 'seq', 'parents', 'address')

    def __init__(self, node_id: str, strand: str, seq: str):
        self.key = f"{node_id}{strand}"
        self.node_id = node_id
        self.strand = strand
        self.seq = seq
        self.parents: List[str] = []
        self.address: Optional[bytes] = None


class BipartiteGraph:
    """
    Expanded bidirected GFA graph.

    Handles all biological DNA/RNA topologies:
        - Antiparallel double helices  : N+ / N- bipartite expansion
        - Prokaryotic circular DNA     : lexicographic canonical cycle anchor
        - Eukaryotic multi-chromosome  : auto-detected disconnected components,
                                         each gets its own HKDF namespace root
        - Viral single-stranded        : single_stranded=True skips N- expansion
        - Viral segmented              : each disconnected subgraph is a segment
                                         under a single Viral ROOT_HASH
    """

    def __init__(self, engine: PanIndexEngine, store: PanIndexStore,
                 single_stranded: bool = False):
        self.engine = engine
        self.store = store
        self.single_stranded = single_stranded
        self.states: Dict[str, BipartiteNode] = {}

    # ------------------------------------------------------------------
    # Parse
    # ------------------------------------------------------------------

    def parse_gfa(self, filepath: str):
        """
        Parse a GFA 1.0 file into the bipartite expanded representation.

        Segment lines create two state nodes (N+, N-).
        Link lines create directed state-to-state edges, plus the
        biological reverse complement edge.
        """
        segments = {}

        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if not parts:
                    continue

                if parts[0] == 'S' and len(parts) >= 3:
                    nid = parts[1]
                    seq = parts[2]

                    fwd = BipartiteNode(nid, '+', seq)
                    self.states[fwd.key] = fwd

                    if not self.single_stranded:
                        rev = BipartiteNode(nid, '-', reverse_complement(seq))
                        self.states[rev.key] = rev

                elif parts[0] == 'L' and len(parts) >= 5:
                    parent_id  = parts[1]
                    parent_dir = parts[2]
                    child_id   = parts[3]
                    child_dir  = parts[4]

                    parent_key = f"{parent_id}{parent_dir}"
                    child_key  = f"{child_id}{child_dir}"

                    # Primary directed edge: parent state -> child state
                    if child_key in self.states and parent_key in self.states:
                        self.states[child_key].parents.append(parent_key)

                    # Biological complement edge (double-stranded only):
                    # If A+ -> B-, then by double-helix symmetry B+ -> A-
                    if not self.single_stranded:
                        rev_parent_key = f"{child_id}{'+' if child_dir == '-' else '-'}"
                        rev_child_key  = f"{parent_id}{'+' if parent_dir == '-' else '-'}"
                        if rev_child_key in self.states and rev_parent_key in self.states:
                            self.states[rev_child_key].parents.append(rev_parent_key)


    # ------------------------------------------------------------------
    # Connected component detection (Eukaryotic multi-chromosome support)
    # ------------------------------------------------------------------

    def _find_connected_components(self) -> List[List[str]]:
        """
        Identify disconnected subgraphs (chromosomes, viral segments).
        Returns list of components, each a list of state keys.

        Used to assign per-component namespace roots:
            Address_Chr1 = HKDF(ROOT_HASH, "component_0")
            Address_Chr2 = HKDF(ROOT_HASH, "component_1")
        """
        visited = set()
        components = []

        # Build undirected adjacency for component detection
        adjacency: Dict[str, set] = {k: set() for k in self.states}
        for key, node in self.states.items():
            # Link states by biological edges (parents)
            for parent in node.parents:
                adjacency[key].add(parent)
                adjacency.setdefault(parent, set()).add(key)
            
            # Explicitly link N+ and N- for component identity
            rev_strand = '-' if node.strand == '+' else '+'
            rev_key = f"{node.node_id}{rev_strand}"
            if rev_key in self.states:
                adjacency[key].add(rev_key)
                adjacency[rev_key].add(key)

        for start in self.states:
            if start in visited:
                continue
            component = []
            stack = [start]
            while stack:
                node = stack.pop()
                if node in visited:
                    continue
                visited.add(node)
                component.append(node)
                for neighbor in adjacency.get(node, []):
                    if neighbor not in visited:
                        stack.append(neighbor)
            components.append(component)

        return components

    # ------------------------------------------------------------------
    # Cycle detection and canonical anchor selection
    # ------------------------------------------------------------------

    def _find_cycle_roots(self) -> Dict[str, str]:
        """
        Detect nodes that are part of a cycle (in-degree > 0 but no path
        to a true root). For each cycle, elect the node with the
        lexicographically smallest sequence as the Local Root.

        Returns dict: state_key -> elected_root_key for cycle members.
        """
        # Build in-degree map
        in_degree = {k: len(v.parents) for k, v in self.states.items()}

        # Kahn's algorithm to identify nodes in cycles
        queue = [k for k, deg in in_degree.items() if deg == 0]
        visited = set()

        while queue:
            current = queue.pop(0)
            visited.add(current)
            for k, node in self.states.items():
                if current in node.parents:
                    in_degree[k] -= 1
                    if in_degree[k] == 0:
                        queue.append(k)

        # Nodes not visited are in cycles
        cycle_nodes = [k for k in self.states if k not in visited]

        if not cycle_nodes:
            return {}

        # Group connected cycle components (simple union-find via adjacency)
        cycle_roots = {}
        components = []
        seen = set()

        for start in cycle_nodes:
            if start in seen:
                continue
            component = []
            stack = [start]
            while stack:
                node = stack.pop()
                if node in seen or node not in self.states:
                    continue
                seen.add(node)
                component.append(node)
                for parent in self.states[node].parents:
                    if parent in cycle_nodes and parent not in seen:
                        stack.append(parent)
            if component:
                components.append(component)

        # Elect canonical root per component: lowest sequence lexicographically
        for component in components:
            elected = min(component, key=lambda k: self.states[k].seq)
            for k in component:
                cycle_roots[k] = elected

        return cycle_roots

    # ------------------------------------------------------------------
    # Address derivation
    # ------------------------------------------------------------------

    def compute_addresses(self, derivation_root: str = "PangenomeRoot"):
        """
        Derive ratchet addresses for all states in the expanded graph.

        - DAG nodes processed in topological order.
        - Cycle nodes resolved by canonical anchor election.
        - Disconnected components (chromosomes/segments) each receive
          their own HKDF namespace:  HKDF(ROOT, derivation_root/component_N)
        """
        # Detect connected components for multi-chromosome / segmented genome support
        components = self._find_connected_components()
        multi_component = len(components) > 1

        # Build per-component namespace roots
        component_roots: Dict[str, bytes] = {}
        base_root = self.engine.derive_ratchet_address(
            self.engine.root_hash, derivation_root
        )

        for comp_idx, component in enumerate(components):
            if multi_component:
                # Each chromosome/segment gets its own namespace
                comp_namespace = f"component_{comp_idx}"
                comp_root = self.engine.derive_ratchet_address(
                    base_root, comp_namespace
                )
            else:
                comp_root = base_root

            for key in component:
                component_roots[key] = comp_root

        cycle_roots = self._find_cycle_roots()

        # Kahn's topological sort on expanded graph
        in_degree = {k: len(v.parents) for k, v in self.states.items()}
        queue = sorted(
            [k for k, deg in in_degree.items() if deg == 0],
            key=lambda k: self.states[k].seq
        )

        for cycle_key, root_key in cycle_roots.items():
            if root_key not in queue:
                queue.append(root_key)

        processed = set()

        while queue:
            key = queue.pop(0)
            if key in processed:
                continue

            node = self.states[key]
            comp_root_addr = component_roots.get(key, base_root)

            if not node.parents or key == cycle_roots.get(key):
                # True root or elected cycle root: use component namespace root
                node.address = self.engine.derive_ratchet_address(
                    comp_root_addr, key
                )
            else:
                resolved_parents = [
                    p for p in node.parents if self.states.get(p) and
                    self.states[p].address is not None
                ]
                if not resolved_parents:
                    queue.append(key)
                    continue

                canonical_parent = min(resolved_parents,
                                       key=lambda p: self.states[p].seq)
                node.address = self.engine.derive_ratchet_address(
                    self.states[canonical_parent].address, key
                )

            processed.add(key)

            for candidate_key, candidate_node in self.states.items():
                if key in candidate_node.parents and candidate_key not in processed:
                    queue.append(candidate_key)

        # Index
        for key, node in self.states.items():
            if node.address:
                comp_idx = next(
                    (i for i, comp in enumerate(components) if key in comp), 0
                )
                tags = [
                    f"strand:{node.strand}",
                    f"node:{node.node_id}",
                    f"component:{comp_idx}",
                ]
                if self.single_stranded:
                    tags.append("genome_type:single_stranded")

                self.store.insert(
                    node_id=key,
                    address=node.address,
                    tags=tags,
                    metadata={
                        'seq': node.seq,
                        'strand': node.strand,
                        'physical_node': node.node_id,
                        'parents': node.parents,
                        'component': comp_idx,
                        'derivation_path': f"{derivation_root}/component_{comp_idx}/{key}",
                    }
                )

    def summary(self):
        addressed = sum(1 for n in self.states.values() if n.address)
        print(f"BipartiteGraph: {len(self.states)} states "
              f"({len(self.states)//2} physical nodes), "
              f"{addressed} addressed")
        for key in sorted(self.states.keys()):
            node = self.states[key]
            addr = node.address.hex()[:12] + '...' if node.address else 'UNRESOLVED'
            rc_flag = ' [RC]' if node.strand == '-' else ''
            print(f"  {key:<8} seq={node.seq[:12]!r:<16} "
                  f"addr={addr}{rc_flag}")


if __name__ == "__main__":
    engine = PanIndexEngine(pangenome_seed=b"bipartite_demo_seed_000000000000")
    store  = PanIndexStore()
    graph  = BipartiteGraph(engine, store)

    gfa_path = os.path.join(os.path.dirname(__file__), '..', 'test.gfa')
    print(f"Parsing: {gfa_path}")
    graph.parse_gfa(gfa_path)
    graph.compute_addresses()
    graph.summary()

    # Verify: forward and reverse states of the same node have different addresses
    node1_fwd = store.get_node("1+")
    node1_rev = store.get_node("1-")
    assert node1_fwd and node1_rev
    assert node1_fwd['address'] != node1_rev['address']
    print("\nVerification: 1+ and 1- have distinct addresses.")

    # Verify: reverse state sequence is the reverse complement of the forward
    fwd_seq = graph.states["1+"].seq
    rev_seq = graph.states["1-"].seq
    assert rev_seq == reverse_complement(fwd_seq)
    print("Verification: 1- sequence is reverse complement of 1+.")
    print("\nBipartite Expansion Engine: OK")
