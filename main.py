# main.py

import argparse
import logging
from datetime import datetime
from fasta_io import FastaIO
from assembler import GreedyAssembler
from debruijn_assembler import DeBruijnAssembler

def setup_logger(enable_file_logging):
    """Configures the logging system."""
    logger = logging.getLogger("GenomeAssembler")
    logger.setLevel(logging.DEBUG) # Catch everything internally
    
    # Format for the logs
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s - %(name)s - %(message)s', datefmt='%H:%M:%S')
    
    # Console handler (Terminal) - Only show INFO and above
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter('%(message)s')) # Clean output for terminal
    logger.addHandler(console_handler)
    
    # File handler (Log file) - Save all granular DEBUG info
    if enable_file_logging:
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = f"assembly_run_{timestamp}.log"
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(logging.DEBUG) 
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info(f"Detailed debugging logs are being saved to: {filename}")
        
    return logger

def print_stats(stats):
    """Helper to cleanly print the statistics dictionaries."""
    print(f"\n[{stats['Algorithm']} Assembly Statistics]")
    for key, value in stats.items():
        if key != "Algorithm":
            print(f"  - {key}: {value}")
    print("-" * 40)

def run_pipeline(filepath, file_format, min_overlap, max_mismatches, kmer_size, assembler_choice, enable_log, output_path):
    logger = setup_logger(enable_log)
    
    logger.info(f"Loading sequence data from {filepath} (Format: {file_format.upper()})...")
    reads = FastaIO.read(filepath, file_format=file_format)
    
    if not reads:
        logger.error("No reads to process. Exiting.")
        return
        
    logger.info(f"Successfully loaded {len(reads)} sequences.\n")
    
    # --- 1. Run Greedy Assembler ---
    if assembler_choice in ["greedy", "both"]:
        greedy_assembler = GreedyAssembler(min_overlap=min_overlap, max_mismatches=max_mismatches)
        greedy_contigs, greedy_stats = greedy_assembler.assemble(reads)
        print_stats(greedy_stats)

    # --- 2. Run De Bruijn Assembler ---
    if assembler_choice in ["debruijn", "both"]:
        read_stream = FastaIO.stream(filepath, file_format=file_format)
        debruijn_assembler = DeBruijnAssembler(k=kmer_size)
        debruijn_contigs, debruijn_stats = debruijn_assembler.assemble(read_stream)
        print_stats(debruijn_stats)
        
        # Save the results!
        if output_path and debruijn_contigs:
            success = FastaIO.write(output_path, debruijn_contigs, header_prefix="debruijn_contig")
            if success:
                logger.info(f"Assembled sequence saved to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Assemble genome fragments with full logging.")
    parser.add_argument("filepath", type=str, help="Path to the input sequence file (.fasta or .fastq)")
    parser.add_argument("-f", "--format", type=str, choices=["fasta", "fastq"], default="fastq")
    parser.add_argument("-a", "--assembler", type=str, choices=["greedy", "debruijn", "both"], default="both")
    parser.add_argument("-o", "--overlap", type=int, default=10)
    parser.add_argument("-m", "--mismatches", type=int, default=0)
    parser.add_argument("-k", "--kmer", type=int, default=5)
    parser.add_argument("-s", "--save", type=str, help="Path to save the assembled FASTA file", default="assembled_genome.fasta")
    
    # New Log Argument
    parser.add_argument("-log", "--log", action="store_true", help="Generate a detailed timestamped .log file")
    
    args = parser.parse_args()
    
    run_pipeline(
        filepath=args.filepath, file_format=args.format,
        min_overlap=args.overlap, max_mismatches=args.mismatches,
        kmer_size=args.kmer, assembler_choice=args.assembler,
        enable_log=args.log,
        output_path=args.save  # Pass the save path here
    )


if __name__ == "__main__":
    main()