# PanIndex: Getting Started Guide

This guide covers how to index, query, and visualize pangenomic data using the PanIndex stack.

### 1. Data Ingestion
PanIndex works with standard **GFA 1.0 (Graphical Fragment Assembly)** files. To begin, you must annotate your graph with PanIndex tags.

```bash
# Annotate a standard GFA file
python src/annotator.py data/input.gfa data/annotated.gfa
```

This command computes addresses for every node and injects two new tags:
- `AN:Z:<hex>`: The 32-byte PanIndex ratchet address.
- `PA:Z:<path>`: The human-readable derivation path.

### 2. Querying the Pangenome
The `query.py` interface provides three distinct search modes:

#### A. Ratchet Search (Path-based)
Find a node by its categorical path. Best for known chromosomal positions.
```bash
python src/query.py --gfa data/annotated.gfa --path "PangenomeRoot/Node14"
```

#### B. Tag Search (Keyword-based)
Find nodes by Anubandha tags (e.g., AMR gene names, variant IDs).
```bash
python src/query.py --gfa data/annotated.gfa --tag "AN:OXA-48"
```

#### C. Similarity Search (Sequence-based)
Find nodes that contain a specific DNA sequence, even with mutations, using MinHash/LSH.
```bash
python src/query.py --gfa data/annotated.gfa --similarity "ATGCG..."
```

### 3. Handling Real Biological Quirks
PanIndex handles complex biology natively:
- **Plasmids**: Circular DNA is automatically "unrolled" using lexicographic anchoring.
- **HGT**: Horizontal Gene Transfer is tracked via topological "symlinks" between strain subgraphs.
- **Double Strands**: Every query is strand-aware. Matching the reverse complement will return a different, strand-specific address.

### 4. Visualization
To see the graph and the PanIndex metadata in your browser:

1. Start the API server:
   ```bash
   python viz/server.py
   ```
2. Open `viz/index.html` in any modern browser.
3. Use **Semantic Zoom** to move between Architectural (HGT/Strains), Functional (Gene clusters), and Atomic (Base pairs) views.

### 5. Running Tests
To verify your installation and implementation:
```bash
python -m pytest tests/
```
The suite covers 69 checkpoints across all biological and mathematical layers.
