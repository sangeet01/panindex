import hashlib
import os
import binascii
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.kdf.hkdf import HKDF
from typing import TYPE_CHECKING, List, Optional, Union

if TYPE_CHECKING:
    from meta_layer import DerivationHistory


class PanIndexEngine:
    """
    PanIndex Core Hashing Engine.
    Implements Commutative Merkle Addressing and Fractal Ratchet HKDF Derivation.

    Optionally accepts a DerivationHistory instance (Asiddhatva staged tracking).
    When attached, every derive_ratchet_address() call is recorded so the full
    derivation path of any node can be replayed or used for lift-over.
    """

    def __init__(
        self,
        pangenome_seed: Optional[bytes] = None,
        history: Optional['DerivationHistory'] = None,
    ):
        if pangenome_seed is None:
            self.root_seed = os.urandom(32)
        else:
            self.root_seed = pangenome_seed

        # Pangenome Root Hash (H_root)
        self.root_hash = hashlib.sha256(self.root_seed).digest()

        # Optional Asiddhatva derivation history log
        self.history: Optional['DerivationHistory'] = history

        # Running depth counter for history entries
        self._depth: int = 0

    def derive_ratchet_address(
        self,
        parent_hash: bytes,
        context: str,
        node_id: Optional[str] = None,
    ) -> bytes:
        """
        Fractal Ratchet - HKDF derivation for a known path component.

        One component is derived in constant time; a path with depth ``d``
        costs O(d) before the resulting address is looked up in an index.

        If a DerivationHistory is attached, the step is recorded for
        Asiddhatva staged traceability.  node_id defaults to context when
        not supplied.
        """
        hkdf = HKDF(
            algorithm=hashes.SHA256(),
            length=32,
            salt=parent_hash,
            info=context.encode(),
        )
        address = hkdf.derive(parent_hash)

        if self.history is not None:
            self.history.record(
                node_id=node_id if node_id is not None else context,
                parent_hash=parent_hash,
                address=address,
                context=context,
                depth=self._depth,
            )
            self._depth += 1

        return address

    def compute_content_address(self, sequence: str, strand: str = '+') -> bytes:
        """Return a sequence identity independent of graph placement."""
        normalized = sequence.upper().encode()
        return hashlib.sha256(normalized + b'|' + strand.encode()).digest()

    def compute_topology_address(
        self,
        sequence: str,
        parent_addresses: List[bytes],
        strand: str = '+',
    ) -> bytes:
        """Return an order-independent identity for sequence placement."""
        content_id = self.compute_content_address(sequence, strand)
        neighborhood = self.compute_commutative_hash(parent_addresses)
        return hashlib.sha256(
            b'FRX-TOPOLOGY\x00' + content_id + neighborhood
        ).digest()

    def compute_commutative_hash(self, neighbor_hashes: List[bytes]) -> bytes:
        """
        Commutative neighborhood hash using XOR (Zobrist-inspired).
        Invariant to the order of neighbors.
        """
        if not neighbor_hashes:
            return b'\x00' * 32
        
        # XOR all neighbor hashes together in a high-dimensional space
        # Here we use 256-bit hashes as our vectors
        res = int.from_bytes(neighbor_hashes[0], 'big')
        for h in neighbor_hashes[1:]:
            res ^= int.from_bytes(h, 'big')
            
        return res.to_bytes(32, 'big')

    def compute_node_address(self, sequence: str, neighborhood_context_hash: bytes) -> bytes:
        """
        Combined Node Address: content + structural context.
        """
        seq_hash = hashlib.sha256(sequence.encode()).digest()
        # Mix sequence and neighborhood context
        return hashlib.sha256(seq_hash + neighborhood_context_hash).digest()

    def canonical_cycle_hash(self, circular_sequence: str) -> bytes:
        """
        Canonical Cycle Hashing for plasmids.
        Unrolls circular DNA by finding the lexicographical minimum k-mer.
        """
        # Lexicographical minimum cyclic shift
        n = len(circular_sequence)
        s = circular_sequence + circular_sequence
        min_s = circular_sequence
        for i in range(1, n):
            current_s = s[i:i+n]
            if current_s < min_s:
                min_s = current_s
        
        return hashlib.sha256(min_s.encode()).digest()

def example_usage():
    engine = PanIndexEngine(pangenome_seed=b"fixed_seed_for_demo_01234567890")
    print(f"Pangenome Root Hash: {binascii.hexlify(engine.root_hash).decode()}")

    # 1. Fractal Ratchet Derivation for a known path
    # Deriving Species -> Chromosome -> Gene
    chr4_addr = engine.derive_ratchet_address(engine.root_hash, "Chr4")
    brca1_addr = engine.derive_ratchet_address(chr4_addr, "BRCA1")
    
    print(f"Chr4 Address:   {binascii.hexlify(chr4_addr).decode()}")
    print(f"BRCA1 Address:  {binascii.hexlify(brca1_addr).decode()}")

    # 2. Split Stability Verification
    # Initial situation: Node N has neighbors A and B
    hash_a = hashlib.sha256(b"neighbor_a").digest()
    hash_b = hashlib.sha256(b"neighbor_b").digest()
    
    initial_neighborhood = engine.compute_commutative_hash([hash_a, hash_b])
    print(f"Initial Neighborhood Hash: {binascii.hexlify(initial_neighborhood).decode()}")
    
    # Split event: Node A splits into A1 and A2
    # In a commutative XOR group, we can represent this if we define the group operations carefully.
    # For simplicity in this POC, we show that sequence order doesn't matter:
    alt_neighborhood = engine.compute_commutative_hash([hash_b, hash_a])
    print(f"Order Invariant Hash:      {binascii.hexlify(alt_neighborhood).decode()}")
    
    assert initial_neighborhood == alt_neighborhood
    print("Verification: Neighborhood hash is order-invariant.")

if __name__ == "__main__":
    example_usage()
