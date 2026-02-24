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

**Read:** `ATGCA`
**k-mers (k=4):** 1. `ATGC`
2. `TGCA`

In the codebase, `DNAEncoder` converts these directly into 2-bit integers. We then count how many times each k-mer appears across all millions of reads. If a k-mer only appears 1 or 2 times, it's almost certainly a sequencing machine error. We throw it away (Controlled by `min_coverage`).

### Step 2: Building the Graph (Nodes & Edges)

For every surviving (trusted) full k-mer, we split it to define our graph.

* **Node (Prefix):** The first $k-1$ characters.
* **Next Node (Suffix):** The last $k-1$ characters.

For the k-mer `ATGC` (length 4):

* Prefix Node: `ATG`
* Suffix Node: `TGC`
* Edge: `C` (The character that transitions `ATG` into `TGC`).

```text
       (Edge: C)
[ATG] ----------> [TGC]

```

By doing this for millions of k-mers, a massive interconnected web is formed.

### Step 3: Fast Graph Traversal

To assemble the genome, we just "walk" along the arrows. However, the graph is not a perfect straight line. It has branches and convergences due to repeating sequences in the genome.

**The Traversal Rules (`_extract_contigs`):**

1. Start at a random active node.
2. Follow the outgoing edge to the next node.
3. Keep walking and recording bases until you hit a **Stop Condition**.

**Stop Condition 1: Outward Branch (Out-degree > 1)**
The path splits. We don't know which way the true genome goes, so we safely stop the contig here.

```text
                  ---> [GCA] (Path A)
                 /
[ATG] ---> [TGC] 
                 \
                  ---> [GCC] (Path B)

```

**Stop Condition 2: Convergence (In-degree > 1)**
Two different paths merge into the same node. If we just blindly kept walking, we might accidentally glue Path X's beginning to Path Y's end, creating a biological **Chimera**. We safely stop.

```text
[ATA] ---> \
            ---> [TGC] ---> [GCA]
[GTA] ---> /

```

*Code Note:* In `debruijn_assembler.py`, we pre-compute an `in_degrees` array and transition `pointers`. This allows the walker to move in $O(1)$ constant time, evaluating millions of steps per second.

### Step 4: Tip Removal (Error Correction)

Even with k-mer filtering, sequencing errors at the very ends of reads can create tiny, dead-end branches in the graph known as "tips".

```text
                  ---> [TGA] ---> [GAC] (True Genome Path)
                 /
[ATG] ---> [TGC] 
                 \
                  ---> [TGT] (Dead end / Error Tip)

```

After traversal, the algorithm looks at all the paths it extracted. If a path is shorter than `2 * k`, it is assumed to be an error tip and is permanently deleted.

### Step 5: Dense K-mer Deduplication (Solving the Double-Strand Problem)

Because the sequencing machine reads from both the forward DNA strand and the reverse complement (RC) strand, our traversal naturally outputs **two complete genomes**: The Forward Genome and the RC Genome.

If the true genome is 2.8 million bases, the traversal outputs 5.6 million bases.

To fix this, the assembler runs `_deduplicate_contigs`:

1. We sort the extracted contigs from longest to shortest.
2. We take the longest contig, accept it, and add **every single one of its k-mers** (and their reverse complements) into a `seen_kmers` database.
3. We move to the next contig. We sample roughly 100 test k-mers across its length.
4. If the majority of those test k-mers already exist in our `seen_kmers` database, we know this contig is just the reverse-complement mirror (or an overlapping redundant fragment) of a contig we already accepted.
5. We safely delete the mirror, cutting the total output size precisely in half to reveal the true haploid genome.