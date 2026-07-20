import sys
import os
from typing import Optional

# Allow running from project root
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from engine import PanIndexEngine
from parser import GFAParser
from index import PanIndexStore


class GFAAnnotator:
    """
    Phase 2 - Normalization Layer Write-Back.

    Reads a standard GFA file, computes a PanIndex address for every node,
    and writes an annotated GFA file with PanIndex optional tags injected:

        AN:Z:<address_hex>  - the 32-byte PanIndex ratchet address (hex)
        PA:Z:<derivation>   - the human-readable derivation path
        AF:i:<tag_count>    - number of Anubandha tags attached

    The output file is valid GFA1 - standard tools (vg, odgi, Bandage)
    ignore unrecognised tags silently.
    """

    def __init__(self, seed: Optional[bytes] = None):
        self.engine = PanIndexEngine(pangenome_seed=seed)
        self.store = PanIndexStore()
        self._parser: Optional[GFAParser] = None

        from default_rules import build_default_rule_engine
        self._rule_engine = build_default_rule_engine()

    def annotate(self, input_path: str, output_path: str,
                 derivation_root: str = "PangenomeRoot"):
        """
        Full pipeline: parse -> address -> index -> write.

        Args:
            input_path      : Path to the source GFA file.
            output_path     : Path where the annotated GFA will be written.
            derivation_root : Label for the top of the derivation path.
        """
        # 1. Parse and compute topology-based Merkle addresses
        self._parser = GFAParser(self.engine)
        self._parser.parse_file(input_path)
        self._parser.compute_all_addresses(root_context=derivation_root)

        # Pre-compute the root ratchet path node used by query_by_path
        root_path_addr = self.engine.derive_ratchet_address(
            self.engine.root_hash, derivation_root
        )

        # 2. Index all nodes
        for node_id, data in self._parser.nodes.items():
            merkle_addr = data.get('panindex_addr')
            if merkle_addr is None:
                continue

            # Pure HKDF path address - mirrors what query_by_path derives
            # Path: root_hash -> derive(derivation_root) -> derive(node_id)
            ratchet_path_addr = self.engine.derive_ratchet_address(
                root_path_addr, node_id
            )

            anubandha_tags = self._extract_anubandha_tags(data.get('tags', []))
            derivation_path = f"{derivation_root}/{node_id}"

            # Apply Paninian rule engine - append resolution as an extra tag
            rule_node = {'tags': anubandha_tags, 'seq': data.get('seq', '')}
            resolution = self._rule_engine.resolve(rule_node, {})
            if resolution and resolution != 'default_resolution':
                anubandha_tags = list(anubandha_tags) + [resolution]

            # Primary registration: by ratchet path address (enables query_by_path)
            self.store.insert(
                node_id=node_id,
                address=ratchet_path_addr,
                tags=anubandha_tags,
                metadata={
                    'seq': data['seq'],
                    'derivation_path': derivation_path,
                    'merkle_addr': merkle_addr.hex(),
                    'out_neighbors': data['out_neighbors'],
                    'in_neighbors': data['in_neighbors'],
                }
            )

        # 3. Write annotated GFA
        self._write_annotated_gfa(input_path, output_path)

    def _extract_anubandha_tags(self, raw_tags: list) -> list:
        """
        Parse GFA optional fields (TAG:TYPE:VALUE) into Anubandha tag strings.
        Example: 'AN:Z:variant_A' -> 'variant_A'
        """
        extracted = []
        for field in raw_tags:
            parts = field.split(':')
            if len(parts) >= 3:
                tag_name = parts[0]
                tag_value = ':'.join(parts[2:])
                extracted.append(f"{tag_name}:{tag_value}")
        return extracted

    def _write_annotated_gfa(self, input_path: str, output_path: str):
        """
        Re-write the GFA file, injecting PanIndex tags onto S-lines.
        All other line types (L, P, W, H) are passed through unchanged.
        """
        with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
            for line in fin:
                parts = line.rstrip('\n').split('\t')

                if parts[0] == 'S' and len(parts) >= 3:
                    node_id = parts[1]
                    node_info = self.store.get_node(node_id)

                    if node_info:
                        addr_hex = node_info['address']
                        derivation = node_info['metadata'].get('derivation_path', '')
                        tag_count = len(node_info['tags'])

                        # Append PanIndex optional tags
                        parts.append(f"AN:Z:{addr_hex}")
                        parts.append(f"PA:Z:{derivation}")
                        parts.append(f"AF:i:{tag_count}")

                fout.write('\t'.join(parts) + '\n')

    def print_index_summary(self):
        print(f"\nIndex Summary: {self.store}")
        print(f"{'Node':<6} {'Address (first 16 chars)':<36} {'Tags'}")
        print("-" * 70)
        for node_id in sorted(self.store.all_nodes()):
            node = self.store.get_node(node_id)
            addr_short = node['address'][:16] + "..."
            tags = ', '.join(node['tags']) if node['tags'] else 'none'
            print(f"{node_id:<6} {addr_short:<36} {tags}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="PanIndex GFA Annotator")
    parser.add_argument('input', help="Input GFA file")
    parser.add_argument('output', help="Output annotated GFA file")
    parser.add_argument('--seed', help="Pangenome seed (string)", default="panindex_default_seed")
    args = parser.parse_args()

    annotator = GFAAnnotator(seed=args.seed.encode())

    print(f"Annotating '{args.input}' -> '{args.output}'")
    annotator.annotate(args.input, args.output)
    annotator.print_index_summary()
    print("\nSUCCESS: GFA annotated and PanIndex tags injected.")
