# Dual-Algorithm Genome Assembler

A robust, pure-Python bioinformatics pipeline for assembling fragmented DNA sequences into continuous genomic contigs. 

This tool provides a side-by-side implementation of two foundational genome assembly algorithms: the **Greedy Overlap (SCS)** approach and the **De Bruijn Graph** approach. It is designed to be highly configurable, error-tolerant, and observable, featuring real-time progress tracking and extensive debug logging.

## 🧬 Features

* **Multiple File Formats:** Seamlessly reads both FASTA and FASTQ files (powered by Biopython).
* **Two Assembly Engines:** Choose between Overlap-Layout-Consensus (Greedy) or Eulerian Path (De Bruijn) algorithms.
* **Error Tolerance:** Configure allowed mismatches to gracefully handle sequencing errors.
* **Rich Observability:** Live progress bars via `tqdm` and deep debug logging with Python's native `logging` module.
* **Production-Grade Testing:** Backed by a comprehensive `pytest` suite handling biological edge cases like circular genomes and tandem repeats.

---

## ⚙️ Installation & Setup

1. **Clone or download the repository.**
2. **Set up a virtual environment (Recommended):**

```bash
   python -m venv venv
   source venv/bin/activate  # On Windows use: venv\Scripts\activate

```

3. **Install dependencies:**

```bash
pip install -r requirements.txt

```

*(Requires Python 3.6+)*

---

## 🔬 How the Algorithms Work

### 1. The Greedy Overlap Assembler

This algorithm solves the **Shortest Common Superstring (SCS)** problem. It operates like a jigsaw puzzle, looking at entire reads and comparing them against one another to find the longest exact (or near-exact) overlapping regions.

**How it works:**

1. Compares every single fragment to every other fragment in the pool.
2. Finds the pair with the absolute longest overlapping match.
3. Merges the two fragments into a single sequence.
4. Repeats until no fragments share an overlap greater than the minimum threshold.

**Visualizing the Merge:**

```text
Fragment 1:  A T G C G T A
Fragment 2:        C G T A C G G
                   | | | |
Merged:      A T G C G T A C G G

```

* **Pros:** Highly intuitive, tolerates sequencer mismatches well, and preserves the continuity of long reads.
* **Cons:** Computationally expensive. The time complexity is , meaning it becomes exponentially slower on massive datasets.

### 2. The De Bruijn Graph Assembler

Instead of comparing whole reads, this algorithm chops all reads into uniform, tiny sliding windows called **-mers**. It uses these -mers to build a massive directed graph, naturally collapsing redundant data and completely ignoring which original read a sequence came from.

**How it works:**

1. Breaks every read into fragments of length .
2. Creates **Nodes** from the prefixes and suffixes of these -mers (length ).
3. Creates **Edges** using the -mers themselves to connect the nodes.
4. Reconstructs the genome by finding an **Eulerian Path**—a continuous route that visits every single edge in the graph exactly once.

**Visualizing the Graph (Sequence: ATGC, k=3):**

```text
3-mers generated:  ATG, TGC

Nodes (2-mers):    [AT], [TG], [GC]
Edges (3-mers):      (ATG)   (TGC)

Graph Structure:   [AT] -----> [TG] -----> [GC]

```

* **Pros:** Extremely fast and scalable. Time complexity is roughly , making it the industry standard for assembling massive genomes (like the human genome).
* **Cons:** Sensitive to -mer size. If  is too small, repeating genomic regions cause "tangles" (loops) in the graph. Cannot easily tolerate sequencer mismatches without complex "bubble popping" logic.

#### 📏 Choosing the Optimal $k$-mer Size 

Choosing the correct $k$ value is the single most important decision when running a De Bruijn graph assembler. There is no single "magic number"—it depends entirely on your dataset's read length, sequencing error rate, and the organism's genomic complexity (repetitive regions).

##### The "Goldilocks" Trade-off

* **If $k$ is too small (e.g., $k=15$ or $21$):** * **High Sensitivity, Low Specificity.** * The graph easily connects reads, but it cannot span across repetitive sequences in the genome. These repeats collapse into massive "roundabouts" (tangles) in the graph. The algorithm gets stuck looping through them, resulting in a "hairball" graph and an assembled sequence that is wildly bloated (often an order of magnitude larger than the real genome).
* **If $k$ is too large (e.g., $k=99$ on a 100bp read):** * **High Specificity, Low Sensitivity.** * Reads must match almost perfectly to form an edge. A single sequencing error or a tiny gap in coverage will completely shatter the graph. Your output will be hundreds or thousands of tiny, disconnected contigs instead of one continuous genome.

##### Best Practices & Rules of Thumb

1.  **Always use an odd number:** Because DNA is double-stranded, an even-numbered $k$-mer can accidentally be its own reverse-complement (a palindrome). This causes the graph to fold in on itself mathematically. Always pick odd numbers (e.g., 21, 31, 55, 77).
2.  **Aim for 50% to 70% of your read length:** If your sequencing machine outputs reads that are 100 bases long, a good starting point for $k$ is usually between `51` and `71`. 
3.  **To resolve repeats, $k$ must be larger than the repeat:** If the genome has a tandem repeat of length $L$, your $k$ must be at least $L + 1$ to step completely over it and prevent a cycle.

##### The $k$-mer Sweep Strategy

Professional bioinformaticians rarely run an assembly just once. They perform a **$k$-mer sweep**, running the algorithm multiple times to find the mathematical sweet spot.

**How to sweep using this tool:**

Run the CLI with varying odd numbers and watch the `[De Bruijn Graph Assembly Statistics]` output:

```bash
python main.py data.fastq -a debruijn -k 31
python main.py data.fastq -a debruijn -k 45
python main.py data.fastq -a debruijn -k 61
```

---

## 💻 Usage & Command Line Interface

Run the pipeline via `main.py`. The tool features a dynamic CLI built with `argparse`.

### Basic Usage

Run both algorithms with default settings on a FASTQ file:

```bash
python main.py sample.fastq

```

### Advanced Usage

Run only the Greedy assembler, requiring a 20-base overlap and allowing 1 mismatch, outputting a debug log:

```bash
python main.py sample.fasta -f fasta -a greedy -o 20 -m 1 -log

```

### Full Options Reference

| Argument | Flag | Description | Default |
| --- | --- | --- | --- |
| **`filepath`** | (Positional) | The path to your input sequence file. | *Required* |
| **Format** | `-f`, `--format` | Specifies if the file is `fasta` or `fastq`. | `fastq` |
| **Assembler** | `-a`, `--assembler` | Choose which algorithm to run: `greedy`, `debruijn`, or `both`. | `both` |
| **Min Overlap** | `-o`, `--overlap` | (Greedy only) The minimum base overlap required to merge two reads. | `10` |
| **Mismatches** | `-m`, `--mismatches` | (Greedy only) Max sequencing errors allowed in the overlap region. | `0` |
| **K-mer Size** | `-k`, `--kmer` | (De Bruijn only) The length of the -mers used to build the graph. | `5` |
| **Logging** | `-log`, `--log` | Flag to generate a detailed, timestamped `.log` file of the run. | `False` |

---

## 📊 Understanding the Output

When a run finishes, the tool prints a high-level statistical summary to your terminal:

```text
[De Bruijn Graph Assembly Statistics]
  - Settings: k-mer size: 21
  - Initial Reads: 1000
  - Total Edges: 4500
  - Main Genome Size: 4520
  - Data Thrown Out: 2 short reads, 0 unvisited k-mers
----------------------------------------

```

* **Data Thrown Out:** Crucial for QA. In Greedy, this is the number of reads that couldn't be merged. In De Bruijn, this tracks reads shorter than  or edges left stranded by graph tangles.

If you used the `-log` flag, a file (e.g., `assembly_run_2026-02-23_17-32-20.log`) is generated containing granular, base-by-base decisions, edge formations, and pathfinding debug data.

---

## 🧪 Running the Test Suite

The project includes a comprehensive, production-grade test suite covering biological edge cases (tandem repeats, circular plasmids) and strict algorithm bounds.

Run the tests using `pytest`:

```bash
pytest test_assemblers.py -v

```
