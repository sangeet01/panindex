"""
FRX Default Paninian Rules.

Provides a pre-configured PaninianRuleEngine with five curated rules that
are applied to every node during annotation. Rules follow Utsarga/Apavada
precedence: higher precedence wins when multiple rules fire.

Rule hierarchy
--------------
 1  genomic_segment     : General rule (Utsarga). Fires for every node.
 3  resistance_candidate: AMR tag detected in any Anubandha tag string.
 5  mobile_element      : Node carries a mobile_element or plasmid tag.
 8  amr_confirmed       : AMR tag present AND sequence length >= 100 bp.
10  critical_resistance : Carbapenem or vancomycin resistance detected.

The resolution string (e.g. "RES:amr_confirmed") is appended as an
additional Anubandha tag on the annotated node. It is written to the
GFA output and stored in the PanIndexStore for tag-based queries.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from meta_layer import PaninianRuleEngine, Rule


# ======================================================================
# Rule conditions
# ======================================================================

def _has_amr_tag(node_data: dict, context: dict) -> bool:
    tags = node_data.get('tags', [])
    return any('AMR' in t or 'bla' in t or 'amr' in t for t in tags)


def _has_mobile_tag(node_data: dict, context: dict) -> bool:
    tags = node_data.get('tags', [])
    return any(
        'mobile_element' in t or 'plasmid' in t or 'PLASMID' in t
        for t in tags
    )


def _amr_and_long_seq(node_data: dict, context: dict) -> bool:
    return _has_amr_tag(node_data, context) and len(node_data.get('seq', '')) >= 100


def _has_critical_tag(node_data: dict, context: dict) -> bool:
    tags = node_data.get('tags', [])
    keywords = ('carbapenem', 'vancomycin', 'mcr', 'NDM', 'OXA-48', 'KPC')
    return any(any(kw in t for kw in keywords) for t in tags)


# ======================================================================
# Factory
# ======================================================================

def build_default_rule_engine() -> PaninianRuleEngine:
    """
    Build and return a PaninianRuleEngine with all five default FRX rules.

    Rules are sorted internally by precedence (highest wins). Callers
    should create one instance at annotator startup and reuse it.

    Returns:
        Configured PaninianRuleEngine ready for resolve() calls.
    """
    engine = PaninianRuleEngine()

    # Precedence 1: General (Utsarga) - applies to everything
    engine.add_rule(Rule(
        name="genomic_segment",
        precedence=1,
        condition=lambda n, c: True,
        action=lambda n, c: "RES:genomic",
    ))

    # Precedence 3: AMR candidate detected via tag string
    engine.add_rule(Rule(
        name="resistance_candidate",
        precedence=3,
        condition=_has_amr_tag,
        action=lambda n, c: "RES:amr_candidate",
    ))

    # Precedence 5: Mobile genetic element (plasmid/transposon tag)
    engine.add_rule(Rule(
        name="mobile_element",
        precedence=5,
        condition=_has_mobile_tag,
        action=lambda n, c: "RES:mobile",
    ))

    # Precedence 8: AMR confirmed - tag present AND sequence >= 100 bp
    engine.add_rule(Rule(
        name="amr_confirmed",
        precedence=8,
        condition=_amr_and_long_seq,
        action=lambda n, c: "RES:amr_confirmed",
    ))

    # Precedence 10: Critical resistance (Apavada exception)
    engine.add_rule(Rule(
        name="critical_resistance",
        precedence=10,
        condition=_has_critical_tag,
        action=lambda n, c: "RES:critical",
    ))

    return engine
