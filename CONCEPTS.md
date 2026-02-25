# Understanding the Vectorized De Bruijn Genome Assembler

This document breaks down the biological terminology, graph theory concepts, and algorithmic steps required to understand how this genome assembler works.

---

## 🧬 Part 1: Biological & Sequencing Terminology

Before building a graph, you must understand the data we are working with.

* **Genome:** The complete set of DNA instructions found in a cell. DNA is made of four chemical bases: Adenine (**A**), Cytosine (**C**), Guanine (**G**), and Thymine (**T**).
* **Double-Stranded / Reverse Complement:** DNA consists of two twisted strands. **A** always pairs with **T**, and **C** always pairs with **G**. Furthermore, the strands run in opposite directions.
* If Strand 1 is `A-T-G-C`, Strand 2 is `G-C-A-T`. We call Strand 2 the **Reverse Complement (RC)**. Sequencers read from *both* strands randomly.
* **Reads:** Modern sequencing machines cannot read an entire genome (millions of bases) from start to finish. Instead, they shatter the genome into millions of tiny fragments (usually 100-250 bases long) and read those. These fragments are called "reads".
* **Coverage / Depth:** Because we randomly shatter the genome, we need to sequence it many times over to ensure every piece is captured. If a genome is 1 million bases long, and we generate 15 million bases worth of reads, we have **15x Coverage**.
* **Contig:** A contiguous (continuous) sequence of DNA created by computationally overlapping our short reads. The goal of an assembler is to merge millions of reads into a handful of massive contigs.
* **Chimeras:** A biological error in assembly where two completely unrelated parts of the genome are accidentally glued together.

---

## 💻 Part 2: Graph Theory & Computing Concepts

To assemble the reads, we use a data structure called a **De Bruijn Graph**.

* **k-mer:** A substring of length $k$. If $k=3$, the k-mers of `ATGC` are `ATG` and `TGC`.
* **Node (Vertex):** A point in the graph. In a De Bruijn graph, a node is a k-mer of length $k-1$.
* **Edge:** A directional line connecting two nodes. In a De Bruijn graph, an edge represents a full k-mer of length $k$, connecting its prefix (first $k-1$ characters) to its suffix (last $k-1$ characters).
* **In-Degree & Out-Degree:** * **In-Degree:** The number of arrows (edges) pointing *into* a node.
* **Out-Degree:** The number of arrows (edges) pointing *out* of a node.


* **Vectorization (NumPy):** Standard Python loops are slow. Vectorization means performing mathematical operations on entire arrays of data simultaneously in C-memory, skipping Python's slow evaluation loop.
* **2-Bit Encoding:** Storing DNA as text (strings) wastes memory. Because there are only 4 bases, they can be perfectly represented by 2 binary bits: `A=00`, `C=01`, `G=10`, `T=11`. `ATGC` becomes `00111001`, drastically shrinking the RAM required.

---

## ⚙️ Part 3: The De Bruijn Algorithm (Step-by-Step)

Here is exactly how `debruijn_assembler.py` transforms raw sequence reads into a finalized genome.

### Step 1: K-mer Extraction and Filtering (The "Shredder")

Instead of trying to find overlapping reads (which is mathematically too slow for millions of reads), we chop every single read into overlapping uniform strings of length $k$ (e.g., $k=4$).

We then count how many times each k-mer appears across the entire dataset. If a k-mer only appears a few times, it is a random machine error and is discarded.

```text
RAW READ: A T G C A T
(Chopping with k=4)

1. [A T G C] (Count: 25) -> TRUSTED (Keep)
2.   [T G C A] (Count: 20) -> TRUSTED (Keep)
3.     [G C A T] (Count: 2)  -> NOISE (Discard & destroy)

```

### Step 2: Building the Graph (Nodes & Edges)

For every surviving (trusted) full k-mer, we split it to define our graph.

* **Node (Prefix):** The first $k-1$ characters.
* **Next Node (Suffix):** The last $k-1$ characters.
* **Edge:** The final character that connects them.

```text
Full K-mer: A T G C
            
   Prefix (Node)      Suffix (Next Node)
     [A T G] ------------> [T G C]
                Edge: C

```

By linking millions of these prefix/suffix pairs together, a massive interconnected web is formed.

### Step 3: Fast Graph Traversal & Bubble Popping

To assemble the genome, the algorithm "walks" along the arrows, recording characters as it goes.

**Standard Walk:**

```text
[ATG] ---> [TGC] ---> [GCA] ---> [CAT]  =  A T G C A T

```

However, biology and sequencing aren't perfect. The graph has tangles that the walker must safely navigate:

**Event A: Bubble Popping (Resolving Minor Errors)**
If a single DNA letter is misread by the machine, the graph temporarily splits into two paths, but immediately merges back together a few steps later (a "bubble"). The algorithm detects this, pops the error path, and bridges the gap.

```text
                   ---> [GCA] (Error path) ---
                  /                           \
[ATG] ---> [TGC]                               ---> [CAT] ---> [ATG]
                  \                           /
                   ---> [GTA] (True path)  ---

Result: Bubble Popped! Traversal safely continues.

```

**Event B: True Divergence (Stop Condition 1)**
If a path splits and *doesn't* come back together, it means the graph hit a repeating DNA sequence that exists in multiple different chromosomes. The walker safely stops to avoid guessing.

```text
                   ---> [GCA] ---> [CAA] (To Chromosome 1)
                  /                           
[ATG] ---> [TGC]                               
                  \                           
                   ---> [GCC] ---> [CCT] (To Chromosome 2)

Result: Stop walking. Output the sequence up to this point.

```

**Event C: Convergence (Stop Condition 2)**
If two completely different paths merge into the same node, walking further could accidentally glue Path X's beginning to Path Y's end, creating a biological **Chimera**. The walker safely stops.

```text
[ATA] ---> \
            ---> [TGC] ---> [GCA]
[GTA] ---> /

Result: Stop walking to prevent creating a Chimera.

```

### Step 4: Tip Removal (Error Pruning)

Even with k-mer filtering, sequencing errors at the very *ends* of reads can create tiny, dead-end branches in the graph known as "tips".

Because a single error corrupts exactly $k$ k-mers, these tips are usually between length $k$ and $2k$. After traversal, the algorithm measures all extracted paths. If a path is a dead end and shorter than $2k$, it is vaporized.

```text
                   ---> [TGA] ---> [GAC] ---> [ACT] (True Genome Path, Length > 2k. KEEP)
                  /
[ATG] ---> [TGC] 
                  \
                   ---> [TGT] (Dead end Tip, Length < 2k. PRUNE & DELETE)

```

### Step 5: Dense K-mer Deduplication (Solving the Double-Strand Problem)

Because the sequencing machine reads from both the forward DNA strand and the reverse complement (RC) strand, our traversal naturally outputs **two complete genomes**: The Forward Genome and the RC Mirror Genome.

To isolate the true haploid genome, the deduplicator compares all contigs against each other. If a contig is found to be a mirror of one we already have, it is discarded.

```text
Extracted Contigs:

Contig 1: [A T G C G T A C G T A G]  --> KEEP (Forward Strand)
Contig 2: [C T A C G T A C G C A T]  --> DISCARD (100% Reverse Complement match to Contig 1)
Contig 3: [G G G C C C A A A T T T]  --> KEEP (New Sequence)

Final Output: Contig 1, Contig 3.

```