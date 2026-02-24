# Vectorized Genome Assembler

A highly memory-optimized, lightning-fast *De Bruijn Graph* genome assembler built in Python.

This assembler is designed to process high-throughput sequencing data (FASTA/FASTQ) using advanced NumPy vectorization, implicit bitmasks, and O(1) graph traversal. It specifically targets memory bottlenecks common in graph-based assembly by using strict 2-bit DNA encoding and post-assembly Dense K-mer Deduplication.

## 🚀 Key Features

* **Vectorized De Bruijn Graph:** Bypasses slow Python dictionaries. Nodes and edges are constructed using raw C-level NumPy arrays and bitwise operations.
* **Extreme Memory Efficiency (2-bit Encoding):** DNA sequences are converted into 64-bit integers (`A=00, C=01, G=10, T=11`). To save 50% of active RAM during graph construction, only the forward strand is built into the graph.
* **O(1) Graph Traversal:** Pre-computes transition pointers for all millions of nodes simultaneously, completely eliminating binary searches during the contig extraction phase.
* **Chimera Prevention (In-Degree Tracking):** Tracks converging paths (in-degree > 1) to stop traversal before merging unrelated genomic regions into fake "chimeric" contigs.
* **Graph Error Correction (Tip Removal):** Automatically detects and drops dead-end error branches (tips) shorter than `2 * k` before they are decoded.
* **Dense K-mer Deduplication:** Seamlessly identifies and removes reverse-complement (RC) and redundant overlapping contigs using a 100-stride dense k-mer anchor check.
* **Automated Logging System:** Generates highly detailed, timestamped debugging logs routed automatically to a generated `logs/` directory.

## 📋 Requirements

The project relies on a few core data science and bioinformatics libraries:

* Python 3.8+
* `numpy` (For vectorized bitwise math)
* `biopython` (For memory-efficient FASTQ/FASTA streaming)
* `tqdm` (For CLI progress bars)

Install the dependencies using:

```bash
pip install numpy biopython tqdm

```

## 🛠️ Usage

Run the assembler via the command line using `main.py`.

### Basic Command

```bash
python main.py ./data/your_reads.fastq -a debruijn -k 31 -log

```

### Full Command-Line Arguments

| Argument | Short | Default | Description |
| --- | --- | --- | --- |
| `filepath` |  | *(Required)* | Path to the input sequence file (`.fasta` or `.fastq`). |
| `--format` | `-f` | `fastq` | Specify the format of the input file (`fasta` or `fastq`). |
| `--assembler` | `-a` | `debruijn` | Algorithm to use (currently optimized for `debruijn`). |
| `--kmer` | `-k` | `5` | K-mer size for the De Bruijn graph. For large genomes, use `31` to `55`. |
| `--save` | `-s` | `assembled_genome.fasta` | Output file path for the final assembled contigs. |
| `--log` | `-log` | `False` | Flag to enable detailed timestamped debugging logs in the `logs/` directory. |

*Note: The legacy Greedy Assembler overlaps (`-o`) and mismatches (`-m`) arguments are currently deprecated in favor of the optimized De Bruijn pipeline.*

## 🧠 How the Pipeline Works

1. **Streaming & Encoding:** Biopython streams reads one at a time to prevent memory overflow. The `DNAEncoder` bit-packs the forward strand of the read into an integer array.
2. **Noise Filtration:** Identical k-mers are grouped. Any k-mer failing to meet the minimum coverage threshold (default: `15`) is discarded as sequencing noise.
3. **Edge Extraction:** Valid k-mers are split into prefixes (nodes) and suffixes (edges) using bit shifts.
4. **O(1) Traversal:** A pointer matrix is built. The graph walks linear paths from active nodes until it hits a dead end, an outward branch (out-degree > 1), or a convergence (in-degree > 1).
5. **Decoding & Error Correction:** Paths shorter than `2*k` are deleted. Surviving paths are decoded back into string DNA.
6. **Dense Deduplication:** Because the graph only processed the forward strand, the raw assembly naturally contains both the full forward genome and the full reverse genome. The deduplicator checks every extracted contig against a database of all previously seen k-mers (and their reverse complements) to safely discard the mirror genome, resulting in the true haploid size.

## 📂 Project Structure

```text
genome-assembler/
├── main.py                  # CLI entry point, argument parsing, and pipeline execution
├── debruijn_assembler.py    # Core vectorized logic, DeBruijnAssembler, and DNAEncoder classes
├── fasta_io.py              # Sequence streaming and saving via Biopython
├── logs/                    # Auto-generated folder containing timestamped run logs
└── assembled_genome.fasta   # Default output destination for assembled contigs

```