import os
import pytest
import tempfile
from fasta_io import FastaIO
from assembler import GreedyAssembler
from debruijn_assembler import DeBruijnAssembler

# ==========================================
# FIXTURES (Setup/Teardown for I/O Testing)
# ==========================================

@pytest.fixture
def sample_fasta():
    """Creates a temporary FASTA file and cleans it up after the test."""
    content = ">read1\nATGCGTA\n>read2\nCGTACGG\n"
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as f:
        f.write(content)
        temp_name = f.name
    yield temp_name
    os.remove(temp_name)

@pytest.fixture
def sample_fastq():
    """Creates a temporary FASTQ file and cleans it up after the test."""
    content = "@read1\nATGCGTA\n+\nIIIIIII\n@read2\nCGTACGG\n+\nIIIIIII\n"
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fastq') as f:
        f.write(content)
        temp_name = f.name
    yield temp_name
    os.remove(temp_name)

# ==========================================
# FILE I/O TESTS
# ==========================================

class TestFastaIO:
    
    def test_read_valid_fasta(self, sample_fasta):
        reads = FastaIO.read(sample_fasta, file_format="fasta")
        assert len(reads) == 2
        assert reads[0] == "ATGCGTA"
        assert reads[1] == "CGTACGG"

    def test_read_valid_fastq(self, sample_fastq):
        reads = FastaIO.read(sample_fastq, file_format="fastq")
        assert len(reads) == 2
        assert reads[0] == "ATGCGTA"

    def test_read_missing_file(self):
        """Production code should not crash on missing files, it should handle gracefully."""
        reads = FastaIO.read("nonsense_path.fastq", file_format="fastq")
        assert reads == []

    def test_invalid_format_parameter(self, sample_fasta):
        """Testing how the parser handles a mismatched format expectation."""
        reads = FastaIO.read(sample_fasta, file_format="fastq") 
        # Biopython will fail to parse a FASTA as a FASTQ, should return empty list
        assert reads == []

def test_assembly_metrics():
    from metrics import calculate_assembly_metrics # adjust import as needed
    
    # Total length = 100. Half length = 50.
    contigs = [
        "A" * 40, # 40
        "G" * 30, # 30  (Running sum: 70 -> Passes 50 mark here)
        "C" * 20, # 20
        "T" * 10  # 10
    ]
    
    metrics = calculate_assembly_metrics(contigs)
    
    # The contig that pushes the sum >= 50 is the 30-length one.
    assert metrics["N50"] == "30"
    # It took 2 contigs to reach the 50% mark
    assert metrics["L50"] == "2"
    # GC content should be 50% (30 Gs + 20 Cs out of 100)
    assert metrics["GC Content"] == "50.00%"


# ==========================================
# GREEDY ASSEMBLER TESTS
# ==========================================

class TestGreedyAssembler:

    @pytest.mark.parametrize("seq1, seq2, expected_overlap, expected_mismatches", [
        ("ATGCGTA", "CGTACGG", 4, 0),  # Perfect overlap
        ("ATGCGTA", "CCTACGG", 4, 1),  # 1 mismatch at the start
        ("ATGCGTA", "CGTTCGG", 4, 1),  # 1 mismatch in the middle
        ("ATGCGTA", "ATGCGTA", 7, 0),  # 100% identical reads
        ("AAAA", "TTTT", 0, 0),        # No overlap at all
    ])
    def test_overlap_calculation_boundaries(self, seq1, seq2, expected_overlap, expected_mismatches):
        """Tests the mathematical boundaries of the overlap logic using parametrization."""
        assembler = GreedyAssembler(min_overlap=3, max_mismatches=1)
        overlap, mismatches = assembler._calculate_overlap(seq1, seq2)
        assert overlap == expected_overlap
        assert mismatches == expected_mismatches

    def test_strict_mismatch_rejection(self):
        """Ensures that overlaps exceeding the mismatch threshold are strictly rejected."""
        assembler = GreedyAssembler(min_overlap=4, max_mismatches=0)
        # "CGTA" vs "CATA" (1 mismatch). Since max=0, it should return 0 overlap.
        overlap, mismatches = assembler._calculate_overlap("ATGCGTA", "CATACGG")
        assert overlap == 0

    def test_circular_genome_handling(self):
        reads = ["ATGC", "GCAT"] 
        assembler = GreedyAssembler(min_overlap=2, max_mismatches=0)
        # FIX: Unpack the tuple
        contigs, _ = assembler.assemble(reads)
        assert len(contigs) == 1
        assert contigs[0] == "ATGCAT"

    def test_multiple_disconnected_contigs(self):
        reads = ["ATGCGTA", "CGTACGG", "TTTTTTT", "TTTTAAA"]
        assembler = GreedyAssembler(min_overlap=4, max_mismatches=0)
        # FIX: Unpack the tuple
        contigs, _ = assembler.assemble(reads)
        assert len(contigs) == 2
        assert "ATGCGTACGG" in contigs


# ==========================================
# DE BRUIJN ASSEMBLER TESTS
# ==========================================

class TestDeBruijnAssembler:

    def test_reads_shorter_than_k(self):
        reads = ["AT", "TG", "GC"]
        # FIX: Ensure min_coverage=1 doesn't mask the k-mer length failure
        assembler = DeBruijnAssembler(k=5, min_coverage=1) 
        contigs, _ = assembler.assemble(reads)
        assert contigs == []

    def test_tandem_repeat_cycle(self):
        reads = ["ATATA", "TATAC"]
        assembler = DeBruijnAssembler(k=5, min_coverage=1)
        assembler.min_contig_length = 0
        contigs, _ = assembler.assemble(reads)
        
        assert len(contigs) == 1
        assert "ATATAC" in contigs[0]

    def test_perfect_eulerian_path(self):
        reads = ["ATG", "TGC", "GCA"]
        assembler = DeBruijnAssembler(k=3, min_coverage=1)
        assembler.min_contig_length = 0 
        contigs, _ = assembler.assemble(reads)
        assert contigs == ["ATGCA"]

    def test_disjointed_graph(self):
        reads = ["ATG", "TGC", "CCC", "CCG"]
        assembler = DeBruijnAssembler(k=3, min_coverage=1)
        assembler.min_contig_length = 0 
        contigs, _ = assembler.assemble(reads)
        assert len(contigs) == 1 
        assert contigs[0] in ["ATGC", "CCCG"]

# ==========================================
# 5. INTEGRATION TESTS
# ==========================================

def test_end_to_end_pipeline():
    """
    Verifies that the IO and Assembler modules communicate correctly.
    Uses an isolated 16-base mock genome with no repeating 3-mers to 
    safely pass through the strict De Bruijn Tip-Removal and Convergence filters.
    """
    # Read 1: ATGACCCTAGCAA (13 bases)
    # Read 2:      CCTAGCAATCG (11 bases)
    # True Genome: ATGACCCTAGCAATCG (16 bases)
    content = ">read1\nATGACCCTAGCAA\n>read2\nCCTAGCAATCG\n"
    
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as f:
        f.write(content)
        temp_name = f.name
        
    try:
        reads = FastaIO.read(temp_name, file_format="fasta")
        
        greedy = GreedyAssembler(min_overlap=3, max_mismatches=0)
        greedy_result, _ = greedy.assemble(reads)
        
        # k=4 means paths must be at least 8 bases to survive tip removal.
        # Our 16 base sequence easily passes this safely.
        debruijn = DeBruijnAssembler(k=4, min_coverage=1)
        debruijn_result, _ = debruijn.assemble(reads)
        
        assert greedy_result[0] == "ATGACCCTAGCAATCG"
        assert debruijn_result[0] == "ATGACCCTAGCAATCG"
    finally:
        os.remove(temp_name)