from typing import List, Dict, Any, Optional

class Anubandha:
    """
    Invisible Meta-Tags for genomic elements.
    Control behavior and persist across transformations.
    """
    def __init__(self, tag_id: str, value: Any):
        self.tag_id = tag_id
        self.value = value

    def __repr__(self):
        return f"[AN:{self.tag_id}={self.value}]"

class Rule:
    """
    Paninian Rule (Utsarga/Apavada).
    """
    def __init__(self, name: str, precedence: int, condition: callable, action: callable):
        self.name = name
        self.precedence = precedence # Higher = more specific (Apavada)
        self.condition = condition
        self.action = action

class PaninianRuleEngine:
    """
    Resolves graph ambiguities using precedence and semantic filters.
    """
    def __init__(self):
        self.rules: List[Rule] = []

    def add_rule(self, rule: Rule):
        self.rules.append(rule)
        # Keep rules sorted by precedence
        self.rules.sort(key=lambda r: r.precedence, reverse=True)

    def resolve(self, node_data: Dict[str, Any], context: Dict[str, Any]) -> str:
        """
        Apply rules based on precedence to resolve a canonical address.
        """
        for rule in self.rules:
            if rule.condition(node_data, context):
                return rule.action(node_data, context)
        return "default_resolution"

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
