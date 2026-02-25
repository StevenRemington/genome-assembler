try:
    from Bio import SeqIO
except ImportError:
    raise ImportError("Biopython is not installed. Please run 'pip install biopython' first.")

class FastaIO:
    """Handles reading sequence data using Biopython to support FASTA and FASTQ."""
    
    @staticmethod
    def read(filepath, file_format="fastq"):
        """
        Parses a sequence file into a list of sequences using Biopython.
        
        Args:
            filepath (str): The path to the file.
            file_format (str): The format of the file (e.g., 'fasta' or 'fastq').
            
        Returns:
            list[str]: A list of genome sequence fragments as strings.
        """
        reads = []
        try:
            # SeqIO.parse returns an iterator of SeqRecord objects.
            # We iterate through them and extract just the string sequence for now.
            for record in SeqIO.parse(filepath, file_format):
                reads.append(str(record.seq))
                
            return reads
            
        except FileNotFoundError:
            print(f"Error: The file '{filepath}' could not be found.")
            return []
        except ValueError as e:
            print(f"Error parsing file (check if format matches '{file_format}'): {e}")
            return []

    @staticmethod
    def stream(filepath, file_format="fastq"):
        """
        Generates sequence data one read at a time to save memory.
        """
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"Sequence file not found: {filepath}")
            
        try:
            for record in SeqIO.parse(filepath, file_format):
                yield str(record.seq)
        except ValueError as e:
            # Raise a clear diagnostic error if Biopython chokes on the format
            raise ValueError(f"Failed to parse {filepath} as {file_format}. Is the format correct? Details: {e}")
    
    @staticmethod
    def write(filepath, sequences, header_prefix="contig"):
        """Writes a list of sequences to a FASTA file."""
        try:
            with open(filepath, 'w') as f:
                for i, seq in enumerate(sequences):
                    f.write(f">{header_prefix}_{i+1}\n")
                    # Break long sequences into 80-character lines for standard FASTA format
                    for j in range(0, len(seq), 80):
                        f.write(seq[j:j+80] + "\n")
            return True
        except Exception as e:
            print(f"Error writing to file: {e}")
            return False

# --- Quick Test Block ---
if __name__ == "__main__":
    # If you run this file directly, you can test if it reads your FASTQ
    # Make sure to change 'sample.fastq' to your actual file!
    test_reads = FastaIO.read("sample.fastq", file_format="fastq")
    if test_reads:
        print(f"Successfully loaded {len(test_reads)} reads from FASTQ.")
        print(f"First read snippet: {test_reads[0][:50]}...")