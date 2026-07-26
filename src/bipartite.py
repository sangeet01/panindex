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

    The ratchet then derives stable placement addresses on this expanded,
    strictly directed graph with no modification.

Asiddhatva / DerivationHistory:
    When a DerivationHistory is attached to the engine, every
    derive_ratchet_address() call is automatically recorded.  The history
    is stored on the engine and can be inspected after compute_addresses()
    for full staged traceability and lift-over support.

Semantic Filter (Sannidhi):
    _populate_store() builds a neighbor_tags context for each node by
    collecting the Anubandha tags of all parent states.  This context is
    passed to the rule engine so the SemanticFilter can apply proximity
    (Sannidhi) upgrades/downgrades based on adjacent biological signals.
"""
import hashlib
from typing import Dict, List, Optional, Protocol, Tuple
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from engine import PanIndexEngine
from index import PanIndexStore


class RuleEngineProtocol(Protocol):
    def resolve(self, node_data: dict, context: dict) -> str:
        ...


COMPLEMENT = str.maketrans('ATGCatgc', 'TACGtacg')


def reverse_complement(seq: str) -> str:
    """Biological reverse complement of a DNA sequence."""
    return seq.translate(COMPLEMENT)[::-1]


def extract_anubandha_tags(raw_tags: List[str]) -> List[str]:
    """Parse GFA optional fields into Anubandha tag strings."""
    extracted: List[str] = []
    for field in raw_tags:
        parts = field.split(':')
        if len(parts) >= 3:
            tag_name = parts[0]
            tag_value = ':'.join(parts[2:])
            extracted.append(f"{tag_name}:{tag_value}")
    return extracted


class BipartiteNode:
    """
    One state of a physical GFA segment.
    key = "node_id+" or "node_id-"
    """
    __slots__ = ('key', 'node_id', 'strand', 'seq', 'parents', 'address',
                 'content_id', 'topology_id', 'tags')

    def __init__(self, node_id: str, strand: str, seq: str, tags: Optional[List[str]] = None):
        self.key = f"{node_id}{strand}"
        self.node_id = node_id
        self.strand = strand
        self.seq = seq
        self.parents: List[str] = []
        self.address: Optional[bytes] = None
        self.content_id: Optional[bytes] = None
        self.topology_id: Optional[bytes] = None
        self.tags: List[str] = tags or []


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

    Asiddhatva / DerivationHistory:
        Attach a DerivationHistory to the engine before calling
        compute_addresses().  Every ratchet step is then recorded
        automatically via the engine.

    Semantic Filter (Sannidhi):
        _populate_store() collects neighbor tags for each node and passes
        them as 'neighbor_tags' in the BioContext so the SemanticFilter
        can apply proximity-based upgrades and downgrades.
    """

    def __init__(self, engine: PanIndexEngine, store: PanIndexStore,
                 single_stranded: bool = False, rule_engine: Optional[RuleEngineProtocol] = None):
        self.engine = engine
        self.store = store
        self.single_stranded = single_stranded
        self.rule_engine = rule_engine
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
                    seq = parts[2] if parts[2] != '*' else ''
                    raw_tags = parts[3:] if len(parts) > 3 else []
                    anubandha_tags = extract_anubandha_tags(raw_tags)

                    fwd = BipartiteNode(nid, '+', seq, tags=anubandha_tags)
                    self.states[fwd.key] = fwd

                    if not self.single_stranded:
                        rev = BipartiteNode(nid, '-', reverse_complement(seq), tags=anubandha_tags)
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
        Detect nodes that are part of a cycle and elect a canonical root per cycle.

        Uses a forward adjacency dict built once (O(E)) so the full method
        is O(N + E) rather than the previous O(N²).

        Returns dict: state_key -> elected_root_key for cycle members.
        """
        # Build forward adjacency once: parent -> set of children
        children: Dict[str, List[str]] = {k: [] for k in self.states}
        for key, node in self.states.items():
            for parent in node.parents:
                children[parent].append(key)

        # In-degree map
        in_degree = {k: len(v.parents) for k, v in self.states.items()}

        # Kahn’s algorithm — O(N + E)
        queue = [k for k, deg in in_degree.items() if deg == 0]
        visited: set = set()

        while queue:
            current = queue.pop(0)
            visited.add(current)
            for child in children[current]:
                in_degree[child] -= 1
                if in_degree[child] == 0:
                    queue.append(child)

        cycle_nodes = [k for k in self.states if k not in visited]
        if not cycle_nodes:
            return {}

        # Group cycle nodes into connected components via undirected adjacency
        # Build undirected adjacency restricted to cycle nodes only
        cycle_set = set(cycle_nodes)
        undirected: Dict[str, List[str]] = {k: [] for k in cycle_nodes}
        for key in cycle_nodes:
            for parent in self.states[key].parents:
                if parent in cycle_set:
                    undirected[key].append(parent)
                    undirected[parent].append(key)
            for child in children.get(key, []):
                if child in cycle_set:
                    undirected[key].append(child)
                    undirected[child].append(key)

        seen: set = set()
        cycle_roots: Dict[str, str] = {}

        for start in cycle_nodes:
            if start in seen:
                continue
            component: List[str] = []
            stack = [start]
            while stack:
                node = stack.pop()
                if node in seen:
                    continue
                seen.add(node)
                component.append(node)
                for nb in undirected.get(node, []):
                    if nb not in seen:
                        stack.append(nb)
            if component:
                elected = min(component, key=lambda k: self.states[k].seq)
                for k in component:
                    cycle_roots[k] = elected

        return cycle_roots

    # ------------------------------------------------------------------
    # Address derivation
    # ------------------------------------------------------------------

    def compute_addresses(self, derivation_root: str = "PangenomeRoot",
                          populate_store: bool = True):
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
                parent_address = self.states[canonical_parent].address
                assert parent_address is not None
                node.address = self.engine.derive_ratchet_address(
                    parent_address, key
                )

            processed.add(key)

            for candidate_key, candidate_node in self.states.items():
                if key in candidate_node.parents and candidate_key not in processed:
                    queue.append(candidate_key)

        if populate_store:
            self._populate_store(components, derivation_root)

    def _populate_store(self, components: List[List[str]], derivation_root: str):
        """
        Populate the attached store from computed state addresses.

        Sannidhi (proximity) context: for each node, collect the Anubandha
        tags of all parent states and pass them as 'neighbor_tags' in the
        BioContext dict so the SemanticFilter can apply proximity-based
        upgrades and downgrades.

        Asiddhatva: derivation history is recorded automatically by the
        engine during derive_ratchet_address() calls above; here we store
        the derivation_path in metadata so it is queryable from the store.
        """
        root_path_addr = self.engine.derive_ratchet_address(
            self.engine.root_hash, derivation_root
        )

        for key, node in self.states.items():
            if node.address:
                node.content_id = self.engine.compute_content_address(
                    node.seq, node.strand
                )
                parent_addresses = []
                neighbor_tags: List[str] = []
                for parent in node.parents:
                    parent_node = self.states.get(parent)
                    if parent_node and parent_node.address is not None:
                        parent_addresses.append(parent_node.address)
                        neighbor_tags.extend(parent_node.tags)
                node.topology_id = self.engine.compute_topology_address(
                    node.seq, parent_addresses, node.strand
                )
                comp_idx = next(
                    (i for i, comp in enumerate(components) if key in comp), 0
                )

                final_tags = list(node.tags)
                if self.rule_engine is not None:
                    # Build BioContext with Sannidhi neighbor information
                    bio_context = {
                        'neighbor_tags': neighbor_tags,
                        'component': comp_idx,
                        'strand': node.strand,
                    }
                    rule_node = {'tags': final_tags, 'seq': node.seq}
                    resolution = self.rule_engine.resolve(rule_node, bio_context)
                    if resolution and resolution != 'default_resolution':
                        # rule_node['tags'] may have been extended in-place by
                        # SemanticFilter (extra bio-tags); use that enriched list
                        final_tags = list(rule_node['tags']) + [resolution]

                tags = final_tags + [
                    f"strand:{node.strand}",
                    f"node:{node.node_id}",
                    f"component:{comp_idx}",
                ]
                if self.single_stranded:
                    tags.append("genome_type:single_stranded")

                # Retrieve derivation history path for this node if available
                history_path = ''
                if self.engine.history is not None:
                    path_entries = self.engine.history.get_derivation_path(key)
                    history_path = ' -> '.join(path_entries)

                self.store.insert(
                    node_id=key,
                    address=node.address,
                    tags=tags,
                    metadata={
                        'seq': node.seq,
                        'strand': node.strand,
                        'physical_node': node.node_id,
                        'stable_id': node.key,
                        'content_id': node.content_id.hex(),
                        'topology_id': node.topology_id.hex(),
                        'parents': node.parents,
                        'component': comp_idx,
                        'derivation_path': f"{derivation_root}/component_{comp_idx}/{key}",
                        'derivation_history': history_path,
                        'bipartite_state': True,
                    },
                    merkle_addr=node.topology_id.hex(),
                    content_id=node.content_id.hex(),
                )

        # Register physical forward-node aliases so original S-line IDs resolve
        for key, node in self.states.items():
            if node.strand != '+':
                continue

            if node.address and node.content_id is not None and node.topology_id is not None:
                comp_idx = next((i for i, comp in enumerate(components) if key in comp), 0)

                # Collect neighbor tags for Sannidhi context
                neighbor_tags = []
                for parent in node.parents:
                    parent_node = self.states.get(parent)
                    if parent_node:
                        neighbor_tags.extend(parent_node.tags)

                final_tags = list(node.tags)
                if self.rule_engine is not None:
                    bio_context = {
                        'neighbor_tags': neighbor_tags,
                        'component': comp_idx,
                        'strand': node.strand,
                    }
                    rule_node = {'tags': final_tags, 'seq': node.seq}
                    resolution = self.rule_engine.resolve(rule_node, bio_context)
                    if resolution and resolution != 'default_resolution':
                        final_tags = list(rule_node['tags']) + [resolution]

                tags = final_tags + [
                    f"strand:{node.strand}",
                    f"node:{node.node_id}",
                    f"component:{comp_idx}",
                ]
                if self.single_stranded:
                    tags.append("genome_type:single_stranded")

                history_path = ''
                if self.engine.history is not None:
                    path_entries = self.engine.history.get_derivation_path(node.node_id)
                    history_path = ' -> '.join(path_entries)

                path_address = self.engine.derive_ratchet_address(
                    root_path_addr, node.node_id
                )

                self.store.insert(
                    node_id=node.node_id,
                    address=path_address,
                    tags=tags,
                    metadata={
                        'seq': node.seq,
                        'strand': node.strand,
                        'physical_node': node.node_id,
                        'stable_id': node.node_id,
                        'content_id': node.content_id.hex(),
                        'topology_id': node.topology_id.hex(),
                        'parents': node.parents,
                        'component': comp_idx,
                        'derivation_path': f"{derivation_root}/{node.node_id}",
                        'derivation_history': history_path,
                        'is_alias': True,
                    },
                    merkle_addr=node.topology_id.hex(),
                    content_id='',
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
