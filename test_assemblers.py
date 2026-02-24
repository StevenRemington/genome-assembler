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
        """Biological Edge Case: Plasmids/circular genomes wrap around on themselves."""
        reads = ["ATGC", "GCAT"] 
        # ATGC and GCAT overlap by 2 ("GC") to make ATGCAT. 
        # But wait, "AT" also overlaps! Production assemblers shouldn't infinite-loop.
        assembler = GreedyAssembler(min_overlap=2, max_mismatches=0)
        contigs = assembler.assemble(reads)
        assert len(contigs) == 1
        assert contigs[0] == "ATGCAT"

    def test_multiple_disconnected_contigs(self):
        """If data represents two different chromosomes, it should yield two contigs."""
        reads = ["ATGCGTA", "CGTACGG", "TTTTTTT", "TTTTAAA"]
        assembler = GreedyAssembler(min_overlap=4, max_mismatches=0)
        contigs = assembler.assemble(reads)
        
        assert len(contigs) == 2
        # Order doesn't matter, so we check inclusion
        assert "ATGCGTACGG" in contigs
        assert "TTTTTTTAAA" in contigs


# ==========================================
# DE BRUIJN ASSEMBLER TESTS
# ==========================================

class TestDeBruijnAssembler:

    def test_reads_shorter_than_k(self):
        """System Edge Case: Reads shorter than the k-mer size cannot form edges."""
        reads = ["AT", "TG", "GC"]
        assembler = DeBruijnAssembler(k=5)
        # Should return empty or handle gracefully without throwing IndexError
        contigs = assembler.assemble(reads)
        assert contigs == []

    def test_tandem_repeat_cycle(self):
        """Biological Edge Case: Tandem repeats create cycles in the graph."""
        
        # "ATATA" has a repeating "AT" sequence.
        reads = ["ATATA", "TATAC"]
        assembler = DeBruijnAssembler(k=4)
        # Graph: AT->TA, TA->AT, AT->TA, TA->AC
        contigs = assembler.assemble(reads)
        assert len(contigs) == 1
        assert "ATATAC" in contigs[0]

    def test_perfect_eulerian_path(self):
        """Tests standard Eulerian path generation."""
        reads = ["ATG", "TGC", "GCA"]
        assembler = DeBruijnAssembler(k=3)
        contigs = assembler.assemble(reads)
        assert contigs == ["ATGCA"]

    def test_disjointed_graph(self):
        """If the graph is split in two, a basic Eulerian path finder will only return one component."""
        reads = ["ATG", "TGC", "CCC", "CCG"]
        assembler = DeBruijnAssembler(k=3)
        contigs = assembler.assemble(reads)
        # Since our current DeBruijn implementation only finds ONE start node and traces it,
        # it will only assemble one of the two disjointed pieces. 
        assert len(contigs) == 1 
        assert contigs[0] in ["ATGC", "CCCG"]


# ==========================================
# INTEGRATION TESTS
# ==========================================

def test_end_to_end_pipeline(sample_fasta):
    """Verifies that the IO and Assembler modules communicate correctly."""
    reads = FastaIO.read(sample_fasta, file_format="fasta")
    
    greedy = GreedyAssembler(min_overlap=3, max_mismatches=0)
    greedy_result = greedy.assemble(reads)
    
    debruijn = DeBruijnAssembler(k=4)
    debruijn_result = debruijn.assemble(reads)
    
    assert greedy_result[0] == "ATGCGTACGG"
    assert debruijn_result[0] == "ATGCGTACGG"