# PanIndex Technical Specification
## The Pangenomic Fractal Indexing Architecture

PanIndex is a content-addressable, high-performance meta-layer for pangenome graphs. It solves the problem of $O(1)$ coordinate derivation in non-linear, bidirectional, and circular genomic topologies.

### 1. Core Architecture: The Fractal Ratchet
The heart of PanIndex is a hierarchical address derivation system built on **HKDF (HMAC-based Key Derivation Function)**.

- **Deterministic Paths**: A path like `PangenomeRoot/Chromosome_1/Node_14` is converted into a 32-byte address without walking the graph.
- **Hierarchical Chaining**: Each level's address is derived from its parent's address and its own identity:
  $$Addr_{child} = HKDF(Addr_{parent}, ID_{child})$$
- **O(1) Access**: This allows researchers to "compute" the address of any gene or mutation in any strain instantly, bypassing expensive $O(N)$ graph alignments.

### 2. Bipartite Expansion: Solving Bidirectionality
Genomic graphs (GFA) are bidirected because DNA is an antiparallel double helix. This creates topological complexity that breaks standard DAG algorithms.

- **State Doubling**: Every physical segment $N$ is expanded into two mathematical states: $N+$ (forward) and $N-$ (reverse complement).
- **Strand-Aware Addressing**: The ratchet derives distinct addresses for each strand. A search hit automatically identifies the strand orientation because the address itself is strand-specific.

### 3. Canonical Cycle Breaker
Bacteria and viruses often have circular genomes (plasmids), which introduce infinite loops.

- **Symmetry Breaking**: When a cycle is detected ($A \to B \to C \to A$), PanIndex elects a **Canonical Anchor**—the node with the lexicographically smallest DNA sequence.
- **Unrolling**: This node acts as the "Local Root" for that specific cycle, allowing the ratchet to "unroll" the circular structure into a deterministic derivation path.

### 4. Commutative Merkle Graph
To ensure neighborhood stability, nodes are addressed using a combination of their content and their topological context.

- **Zobrist/XOR Hashing**: Neighbor addresses are combined using XOR to ensure that the order in which neighbors are discovered does not change the final node hash (commutativity).
- **Mutation Pinpointing**: A single SNP change only changes the addresses of the affected node and its immediate neighbors, preventing a global "hash landslide."

### 5. Paninian Rule Engine
Inspired by the 2500-year-old grammar of Pāṇini, this engine resolves ambiguity in genomic merges.

- **Utsarga (General) & Apavada (Exception)**: Defines a hierarchy of rules for how sequences should be concatenated or tagged when multiple valid paths exist.
- **Axiomatic Stability**: Ensures that two researchers indexing the same graph will always arrive at the exact same addresses.
