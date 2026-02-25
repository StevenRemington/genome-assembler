# Vectorized Genome Assembler

A highly memory-optimized, lightning-fast *De Bruijn Graph* genome assembler built in Python.

This assembler is designed to process high-throughput sequencing data (FASTA/FASTQ) using advanced NumPy vectorization, implicit bitmasks, and O(1) graph traversal. It specifically targets memory bottlenecks common in graph-based assembly by using strict 2-bit DNA encoding, automated graph cleaning, and Dense K-mer Deduplication.

## 🚀 Key Features

* **Vectorized De Bruijn Graph:** Bypasses slow Python dictionaries. Nodes and edges are constructed using raw C-level NumPy arrays and bitwise operations.
* **Extreme Memory Efficiency (2-bit Encoding):** DNA sequences are compressed into 64-bit integers (`A=00, C=01, G=10, T=11`). To save 50% of active RAM during graph construction, only the forward strand is built into the graph.
* **O(1) Graph Traversal:** Pre-computes transition pointers for millions of nodes simultaneously, completely eliminating binary searches during the contig extraction phase.
* **Advanced Graph Cleaning:** * **Bubble Popping:** Automatically detects and resolves internal Single Nucleotide Polymorphisms (SNPs) and minor sequencing errors to prevent the assembly from shattering.
  * **Tip Removal:** Detects and prunes dead-end error branches (tips) shorter than `2 * k` before they are decoded, preventing artificial inflation of the genome size.
* **Dense K-mer Deduplication:** Seamlessly identifies and removes reverse-complement (RC) and redundant overlapping contigs using a 100-stride dense k-mer anchor check.
* **Biological Quality Metrics:** Automatically calculates and logs standard bioinformatics metrics including **N50, L50, and GC-Content**.
* **Automated Logging & I/O:** Generates highly detailed, timestamped debugging logs and dynamically names output FASTA files based on your input data.

## 📋 Requirements

The project relies on a few core data science and bioinformatics libraries:

* Python 3.8+
* `numpy` (For vectorized bitwise math)
* `biopython` (For memory-efficient FASTQ/FASTA streaming)
* `tqdm` (For CLI progress bars)
* `pytest` (For running the test suite)

Install the dependencies using:

```bash
pip install numpy biopython tqdm pytest

```

## 🛠️ Usage

Run the assembler via the command line using `main.py`.

### Basic Command

```bash
python main.py ./data/your_reads.fastq -a debruijn -k 31 -c 20 --log

```

### Full Command-Line Arguments

| Argument | Short | Default | Description |
| --- | --- | --- | --- |
| `filepath` |  | *(Required)* | Path to the input sequence file (`.fasta` or `.fastq`). |
| `--format` | `-f` | `fastq` | Specify the format of the input file (`fasta` or `fastq`). |
| `--assembler` | `-a` | `both` | Algorithm to use (`greedy`, `debruijn`, or `both`). |
| `--kmer` | `-k` | `31` | K-mer size for the De Bruijn graph. **Max: 31** (due to 64-bit limit). |
| `--min-cov` | `-c` | `20` | Minimum coverage threshold to filter out raw sequencing noise. |
| `--save` | `-s` | *Dynamic* | Path to save the assembled FASTA. Defaults to `<input_name>_assembled.fasta`. |
| `--log` |  | `False` | Flag to enable detailed timestamped debugging logs in the `logs/` directory. |

*Note: The legacy Greedy Assembler overlaps (`-o`) and mismatches (`-m`) arguments are still supported but deprecated for large genomes in favor of the optimized De Bruijn pipeline.*

## 🧠 How the Pipeline Works

1. **Streaming & Encoding:** Biopython streams reads one at a time to prevent memory overflow. The `DNAEncoder` bit-packs the forward strand of the read into an integer array.
2. **Noise Filtration:** Identical k-mers are grouped. Any k-mer failing to meet the minimum coverage threshold (`-c`) is discarded as sequencing noise.
3. **Edge Extraction & Pointers:** Valid k-mers are split into prefixes (nodes) and suffixes (edges) using bit shifts. An $O(1)$ transition pointer matrix is built.
4. **Traversal & Bubble Popping:** The graph walks linear paths. If it hits an outward branch, it looks ahead $k+2$ steps; if the branches reconverge, it "pops the bubble" and bridges the gap.
5. **Tip Removal & Decoding:** Paths shorter than `2*k` are deleted as error tips. Surviving paths are decoded back into string DNA.
6. **Dense Deduplication:** The deduplicator checks every extracted contig against a database of all previously seen k-mers (and their reverse complements) to safely discard the mirror genome, resulting in the true haploid assembly.
7. **Metrics Calculation:** Calculates N50, L50, and GC content, saving the final statistics to the terminal and log files.

## 🧪 Testing Suite

This project includes a comprehensive, deterministic unit testing suite that verifies mathematical boundaries and complex graph anomalies.

To ensure the algorithmic engine is working flawlessly, run:

```bash
pytest test_assemblers.py -v

```

The test suite explicitly proves the validity of:

* Symmetric 2-bit DNA encoding and decoding.
* N50, L50, and GC-Content mathematics.
* SNP resolution (Bubble Popping).
* Dead-end error pruning (Tip Removal).
* Double-strand (Reverse Complement) deduplication.

## 📂 Project Structure

```text
genome-assembler/
├── main.py                  # CLI entry point, argument parsing, and pipeline execution
├── debruijn_assembler.py    # Core vectorized logic, DeBruijnAssembler, and DNAEncoder classes
├── assembler.py             # Legacy Greedy Overlap implementation
├── fasta_io.py              # Sequence streaming and saving via Biopython
├── metrics.py               # Biological quality metrics calculator (N50, L50, GC)
├── test_assemblers.py       # Comprehensive pytest suite for graph anomalies
├── logs/                    # Auto-generated folder containing timestamped run logs
└── outputs/                 # (Optional) Clean directory for assembled FASTA files

```
