
# PanIndex

*Building a searchable index for the tree of life.*

<div align="center">
  <img src="docs/assets/logo.png" width="280" alt="PanIndex Logo" style="border-radius:16px; margin:24px 0; box-shadow: 0 4px 12px rgba(0,0,0,0.1);" />
</div>

---

## 🌍 The Layman's Explanation
Imagine trying to navigate an entire continent using Google Maps, but none of the streets have names, and highways constantly rearrange themselves. 

That is currently what is happening in the science of DNA; Bioinformatics. As we sequence more life on Earth, biologists realize that a single linear string of DNA (like reading a book) isn't enough to capture human diversity. So, they started building massive 3D web-networks called **Pangenome Graphs**. 

But there is a fatal flaw: because it's a massive, tangled web, finding a specific gene requires supercomputers to "walk" the entire web from start to finish. 

**PanIndex solves this.** It is a mathematical engine that gives every single variation, mutation, and gene a permanent, instant address. We took an algorithm originally designed for securing chat protocols (the Fractal Ratchet) and applied it to biology so scientists can pinpoint mutations instantly without supercomputers.

---

## 🔬 The Bioinformatics Bottleneck
For decades, bioinformatics has relied on **FASTA** files—linear strings of text. Because they are 1D lines, finding a mutation is easy: you use a linear coordinate (e.g., `chr1:123,456`). However, this introduces massive **Reference Bias**; if a patient has a genetic sequence that isn't in the reference "string," the software simply ignores it.

To fix this, the field moved to **GFA (Graphical Fragment Assembly)**. A pangenome graph stores overlapping sequences as nodes and paths. But this destroyed our coordinate system. This is known as the **Coordinate Collapse**. If you split a node in a graph, all linear indices shatter. Searching for a Resistance Gene in a 500-isolate bacterial graph suddenly went from simple $O(1)$ indexing to an expensive $O(N)$ graph-alignment problem.

PanIndex is the first **Deterministic, Content-Addressable Meta-Layer** for pangenomic graphs. 
* It completely replaces graph-walking with **O(1) Cryptographic Derivation**.
* It restores stable, FASTA-like pinpointing while keeping 100% of the graph's structural diversity.

---

## ⚙️ The Mathematical Engine
PanIndex does not "assign" IDs to nodes. It mathematically *derives* them.

1. **The Fractal Ratchet ($O(1)$ Lookup):** Computes addresses based on *structural context*. A path like `PangenomeRoot/Chromosome_1/Node_14` mathematically derives a 32-byte address instantly via HMAC-based Key Derivation (HKDF), without walking the graph.
2. **Paninian Meta-Layers:** Inspired by the 2500-year-old linguistics of Pāṇini's Sanskrit grammar, we use invisible meta-tags (Anubandhas) and priority rules (Utsarga/Apavada) to deterministically resolve complex nested variations resulting in a single canonical coordinate for every edge-case.

---

## 🧬 Handling "Impossible" Biological Quirks
DNA isn't straightforward. Standard Directed Acyclic Graphs (DAGs) choke on real-world biology. The PanIndex Engine natively handles life's architectural messiness:

* **Double-Helix Native (Bipartite Expansion):** DNA strands are antiparallel inverses of each other, meaning genomic graphs are inherently *bidirected*. PanIndex expands every physical segment into mathematical states ($N+$ and $N-$), enforcing strict directional flow and tracking strand identity intrinsically within the derived address.
* **Plasmid Loops (Cycle Breakers):** Bacteria have circular DNA (plasmids), creating infinite logic loops that crash algorithms. We apply *Lexicographical Canonical Anchoring*, deterministically "unrolling" circles by electing the alphabetically lowest genetic sequence sequence as the local root.
* **Horizontal Gene Transfer (HGT):** When a bug shares drug-resistance genes with another species, we don't duplicate the data. PanIndex models HGT as a topological symlink—behaving exactly like `git submodule`—so scientists can track the exact evolutionary lineage of a superbug.

---

## ⚡ Quickstart

Get your hands on the code and launch the visualizer using the sample clinical *K. pneumoniae* isolate data.

**1. Add Coordinates to a Raw Graph**
Process a standard `GFA 1.0` graph. PanIndex injects stable 32-byte ratchet addresses right into the sequence tags.

```bash
python src/annotator.py data/kp_pangenome.gfa data/annotated.gfa
```

**2. Query Instantly (No Supercomputing Required)**
Find nodes by categorical path, meta-tag, or raw DNA sequence (using MinHash LSH).

```bash
# Path Search
python src/query.py --gfa data/annotated.gfa --path "PangenomeRoot/Node_42"

# Tag Search (e.g., AMR Resistance Genes)
python src/query.py --gfa data/annotated.gfa --tag "AN:OXA-48"

# Biological Similarity
python src/query.py --gfa data/annotated.gfa --similarity "ATGCGTACT..."
```

**3. Launch the "Warm Tech" Dashboard**
Explore the topological graph in your browser using our custom Cytoscape UI. Toggle between Mathematical Grids, GFA-Style "Hairball" Traces, and Circular Plasmid Rings.

```bash
python viz/server.py
```
*Open `viz/index.html` in any browser to experience the Semantic Zoom interface.*

---

## 📚 Technical Documentation

If you want the deep-dive math, linguistics, and theory, read the canonical documents located in the `docs/` folder:

- 📖 [**The Master Theory**](final.md): The full origin story and Paninian linguistic proofs.
- 🔬 [**Formal Whitepaper**](docs/whitepaper.md): High-level scientific abstract of the implementation.
- ⚙️ [**Technical Specs**](docs/technical_spec.md): Deep dive into the cryptographic hashes and Bipartite rules.

> *Built by Sangeet Sharma.*
