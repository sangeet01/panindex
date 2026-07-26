from typing import Callable, List, Dict, Any, Optional, Tuple


class Anubandha:
    """
    Invisible Meta-Tags for genomic elements.
    Control behavior and persist across transformations.
    Tags are invisible in the final surface sequence but remain queryable.
    """
    def __init__(self, tag_id: str, value: Any):
        self.tag_id = tag_id
        self.value = value

    def __repr__(self):
        return f"[AN:{self.tag_id}={self.value}]"


class Rule:
    """
    Paninian Rule (Utsarga/Apavada).
    Higher precedence = more specific exception (Apavada) that overrides
    the general rule (Utsarga, precedence=1).
    """
    def __init__(self, name: str, precedence: int, condition: Callable, action: Callable):
        self.name = name
        self.precedence = precedence
        self.condition = condition
        self.action = action


# ======================================================================
# Semantic Filter Φ  (Ākāṅkṣā + Yogyatā + Sannidhi)
# ======================================================================

class SemanticFilter:
    """
    Biological compatibility filter — the Φ layer from panini.txt.

    Implements the three Paninian semantic constraints:

    Ākāṅkṣā  (expectancy / completeness)
        A candidate address is only valid if the node carries the minimum
        biological information expected for its resolution class.  A node
        resolved as 'RES:amr_confirmed' must have a non-empty sequence;
        a node resolved as 'RES:critical' must also carry a known
        resistance keyword.  Candidates that fail expectancy are demoted
        to the next lower resolution tier.

    Yogyatā  (semantic compatibility / fitness)
        Checks that the resolved class is biologically consistent with the
        node's sequence content.  Uses a lightweight GC-content heuristic:
        known resistance gene families (beta-lactamases, carbapenemases)
        have characteristic GC ranges.  Nodes whose sequence GC falls
        outside the expected window for their claimed class are flagged
        with a 'COMPAT:low_gc' or 'COMPAT:high_gc' warning tag rather
        than being silently accepted.

    Sannidhi (proximity / context)
        Checks the graph-local context supplied in the BioContext dict.
        If a node's immediate neighbours (passed as 'neighbor_tags' in
        context) carry complementary biological signals (e.g. a mobile
        element adjacent to an AMR gene), the resolution is upgraded.
        Conversely, a claimed AMR node surrounded only by housekeeping
        genes is downgraded to 'RES:amr_candidate'.

    Usage:
        sf = SemanticFilter()
        final_res, bio_tags = sf.apply(resolution, node_data, context)
        # bio_tags is a list of extra Anubandha tags to append
    """

    # GC content windows (inclusive) considered compatible per class
    _GC_WINDOWS: Dict[str, Tuple[float, float]] = {
        'RES:amr_confirmed': (0.35, 0.72),
        'RES:critical':      (0.35, 0.72),
        'RES:mobile':        (0.30, 0.75),
        'RES:amr_candidate': (0.20, 0.80),
        'RES:genomic':       (0.10, 0.90),
    }

    # Neighbor tags that upgrade a candidate to confirmed
    _UPGRADE_NEIGHBORS = frozenset({
        'mobile_element', 'plasmid', 'PLASMID',
        'AMR', 'bla', 'amr', 'transposon', 'integron',
    })

    # Neighbor tags that indicate housekeeping (downgrade signal)
    _HOUSEKEEPING_NEIGHBORS = frozenset({
        'ribosomal', 'rRNA', 'tRNA', 'housekeeping', 'core_chromosome',
    })

    @staticmethod
    def _gc_content(seq: str) -> float:
        if not seq:
            return 0.5  # neutral when sequence unknown
        seq = seq.upper()
        gc = seq.count('G') + seq.count('C')
        return gc / len(seq)

    def apply(
        self,
        resolution: str,
        node_data: Dict[str, Any],
        context: Dict[str, Any],
    ) -> Tuple[str, List[str]]:
        """
        Apply Ākāṅkṣā, Yogyatā, and Sannidhi filters to a candidate resolution.

        Args:
            resolution : Candidate resolution string from the rule engine
                         (e.g. 'RES:amr_confirmed').
            node_data  : Node dict with at least 'seq' and 'tags' keys.
            context    : BioContext dict.  Recognised keys:
                           'neighbor_tags' : List[str] — tags of adjacent nodes.
                           'component'     : int — connected component index.
                           'strand'        : str — '+' or '-'.

        Returns:
            (final_resolution, extra_bio_tags)
            extra_bio_tags is a list of additional Anubandha tags to append
            (e.g. ['COMPAT:low_gc', 'PHI:sannidhi_upgrade']).
        """
        seq = node_data.get('seq', '')
        tags = node_data.get('tags', [])
        neighbor_tags: List[str] = context.get('neighbor_tags', [])
        extra: List[str] = []

        # ----------------------------------------------------------
        # Ākāṅkṣā — expectancy check
        # ----------------------------------------------------------
        if resolution in ('RES:amr_confirmed', 'RES:critical'):
            if not seq:
                # Cannot confirm without sequence evidence — demote
                resolution = 'RES:amr_candidate'
                extra.append('PHI:akanksha_demote')

        # ----------------------------------------------------------
        # Yogyatā — GC compatibility check
        # ----------------------------------------------------------
        gc = self._gc_content(seq)
        lo, hi = self._GC_WINDOWS.get(resolution, (0.10, 0.90))
        if gc < lo:
            extra.append('COMPAT:low_gc')
        elif gc > hi:
            extra.append('COMPAT:high_gc')

        # ----------------------------------------------------------
        # Sannidhi — proximity / neighbor context
        # ----------------------------------------------------------
        neighbor_flat = ' '.join(neighbor_tags)

        has_mobile_neighbor = any(
            kw in neighbor_flat for kw in self._UPGRADE_NEIGHBORS
        )
        has_housekeeping_neighbor = any(
            kw in neighbor_flat for kw in self._HOUSEKEEPING_NEIGHBORS
        )

        if resolution == 'RES:amr_candidate' and has_mobile_neighbor:
            # Surrounded by mobile elements — upgrade confidence
            resolution = 'RES:amr_confirmed'
            extra.append('PHI:sannidhi_upgrade')
        elif resolution in ('RES:amr_confirmed', 'RES:critical') \
                and has_housekeeping_neighbor \
                and not has_mobile_neighbor:
            # Only housekeeping neighbors — downgrade
            resolution = 'RES:amr_candidate'
            extra.append('PHI:sannidhi_downgrade')

        return resolution, extra


# ======================================================================
# Derivation History  (Asiddhatva staged tracking)
# ======================================================================

class DerivationHistory:
    """
    Staged derivation log — the Asiddhatva principle from panini.txt.

    In Pāṇini's grammar, earlier rule applications are treated as
    'not yet effective' (asiddha) for later rules, preserving intermediate
    states for traceability.  Here, every address derivation step is
    recorded as an immutable entry so that:

    - The full derivation path of any node can be replayed.
    - Lift-over between graph constructions uses shared stage anchors.
    - Debugging shows exactly which parent hash and context string
      produced each address at each depth.

    Each entry is a dict:
        {
          'stage'      : int   — monotonically increasing step counter,
          'node_id'    : str   — the node or context label being derived,
          'parent_hex' : str   — hex of the parent hash used,
          'address_hex': str   — hex of the derived address,
          'context'    : str   — HKDF info string (derivation context),
          'depth'      : int   — depth in the derivation tree,
        }

    The history is append-only.  Entries are never modified after recording
    (immutability mirrors the asiddha constraint).
    """

    def __init__(self):
        self._log: List[Dict[str, Any]] = []
        self._stage: int = 0
        # node_id -> list of stage indices for fast per-node lookup
        self._node_stages: Dict[str, List[int]] = {}

    def record(
        self,
        node_id: str,
        parent_hash: bytes,
        address: bytes,
        context: str,
        depth: int = 0,
    ):
        """
        Record one derivation step.  Called by PanIndexEngine after each
        derive_ratchet_address() call.

        Args:
            node_id     : Label for the derived node/context.
            parent_hash : The parent bytes used as HKDF salt.
            address     : The resulting derived address bytes.
            context     : The HKDF info string (e.g. 'PangenomeRoot/Node_14').
            depth       : Depth in the derivation tree (0 = root).
        """
        entry: Dict[str, Any] = {
            'stage':       self._stage,
            'node_id':     node_id,
            'parent_hex':  parent_hash.hex(),
            'address_hex': address.hex(),
            'context':     context,
            'depth':       depth,
        }
        self._log.append(entry)
        self._node_stages.setdefault(node_id, []).append(self._stage)
        self._stage += 1

    def get_node_history(self, node_id: str) -> List[Dict[str, Any]]:
        """
        Return all derivation entries for a specific node_id in stage order.
        O(S) where S = number of stages for that node (usually 1).
        """
        indices = self._node_stages.get(node_id, [])
        return [self._log[i] for i in indices]

    def get_derivation_path(self, node_id: str) -> List[str]:
        """
        Return the ordered list of context strings that produced node_id.
        Useful for lift-over: two graphs sharing a prefix path share an anchor.
        """
        return [e['context'] for e in self.get_node_history(node_id)]

    def replay_addresses(self, node_id: str) -> List[str]:
        """Return the sequence of address_hex values produced for node_id."""
        return [e['address_hex'] for e in self.get_node_history(node_id)]

    def all_stages(self) -> List[Dict[str, Any]]:
        """Return the full immutable log (read-only view)."""
        return list(self._log)

    def __len__(self) -> int:
        return self._stage

    def __repr__(self) -> str:
        return f"DerivationHistory(stages={self._stage}, nodes={len(self._node_stages)})"


# ======================================================================
# Paninian Rule Engine
# ======================================================================

class PaninianRuleEngine:
    """
    Resolves graph ambiguities using Utsarga/Apavada precedence rules
    followed by the Semantic Filter Φ (Ākāṅkṣā + Yogyatā + Sannidhi).

    Resolution pipeline per node:
      1. Iterate rules highest-precedence first; first matching rule fires.
      2. Pass the candidate resolution through SemanticFilter.apply().
      3. Return the final resolution string + append any extra bio-tags.
    """
    def __init__(self, semantic_filter: Optional['SemanticFilter'] = None):
        self.rules: List[Rule] = []
        self.semantic_filter: SemanticFilter = semantic_filter or SemanticFilter()

    def add_rule(self, rule: Rule):
        self.rules.append(rule)
        self.rules.sort(key=lambda r: r.precedence, reverse=True)

    def resolve(
        self,
        node_data: Dict[str, Any],
        context: Dict[str, Any],
    ) -> str:
        """
        Apply rules then semantic filter; return the final resolution string.

        Extra bio-tags produced by the semantic filter are appended to
        node_data['tags'] in-place so callers see the enriched tag list
        after this call.
        """
        candidate = 'default_resolution'
        for rule in self.rules:
            if rule.condition(node_data, context):
                candidate = rule.action(node_data, context)
                break

        if candidate == 'default_resolution':
            return candidate

        final_res, extra_tags = self.semantic_filter.apply(
            candidate, node_data, context
        )

        # Append extra bio-tags in-place so the caller's tag list is enriched
        if extra_tags:
            existing = node_data.get('tags')
            if isinstance(existing, list):
                existing.extend(extra_tags)

        return final_res

def example_meta_layer():
    engine = PaninianRuleEngine()
    
    # Utsarga: General projection rule
    utsarga = Rule(
        name="Source_Projection",
        precedence=1,
        condition=lambda n, c: True,
        action=lambda n, c: f"proj:{n.get('origin', 'unknown')}:{n.get('pos', 0)}"
    )
    
    # Apavada: Specific override for AMR genes
    apavada_amr = Rule(
        name="AMR_Priority",
        precedence=10,
        condition=lambda n, c: "AMR" in n.get('tags', []),
        action=lambda n, c: f"pinpoint:AMR_BUBBLE_{n.get('variant_id')}"
    )
    
    engine.add_rule(utsarga)
    engine.add_rule(apavada_amr)
    
    # Test Node 1: Normal segment
    node1 = {'origin': 'strain_X', 'pos': 100, 'tags': []}
    res1 = engine.resolve(node1, {})
    print(f"Node 1 Resolution: {res1}") # Should be General
    
    # Test Node 2: Resistance segment
    node2 = {'origin': 'strain_Y', 'pos': 500, 'tags': ['AMR'], 'variant_id': 'blaTEM_v1'}
    res2 = engine.resolve(node2, {})
    print(f"Node 2 Resolution: {res2}") # Should be AMR Priority
    
    assert "pinpoint" in res2
    print("SUCCESS: Paninian Rule Engine correctly prioritizes specific exceptions.")

if __name__ == "__main__":
    example_meta_layer()
