# main.py

import os
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
        log_dir = "logs"
        os.makedirs(log_dir, exist_ok=True) # Creates the 'logs' folder if it doesn't exist
        
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        filename = os.path.join(log_dir, f"assembly_run_{timestamp}.log")
        
        file_handler = logging.FileHandler(filename)
        file_handler.setLevel(logging.DEBUG) 
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        logger.info(f"Detailed debugging logs are being saved to: {filename}")
        
    return logger

def print_stats(stats, logger):
    """Cleanly formats statistics and ensures they are written to the log file."""
    if not stats:
        logger.error("No statistics generated. Assembly may have failed.")
        return

    logger.info(f"\n[{stats.get('Algorithm', 'Unknown')} Assembly Statistics]")
    for key, value in stats.items():
        if key != "Algorithm":
            logger.info(f"  - {key}: {value}")
    logger.info("-" * 40)

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
        try:
            read_stream = FastaIO.stream(filepath, file_format=file_format)
            debruijn_assembler = DeBruijnAssembler(k=kmer_size, min_coverage=args.min_cov)
            debruijn_contigs, debruijn_stats = debruijn_assembler.assemble(read_stream)
            print_stats(debruijn_stats, logger)
        except Exception as e:
            logger.error(f"Assembly pipeline failed: {e}")
            return
        
        # Save the results ONLY if contigs were generated
        if output_path and debruijn_contigs:
            success = FastaIO.write(output_path, debruijn_contigs, header_prefix="debruijn_contig")
            if success:
                logger.info(f"Assembled sequence saved to: {output_path}")
        else:
            logger.warning("No contigs were generated. Output file was not created.")

def main():
    parser = argparse.ArgumentParser(description="Assemble genome fragments with full logging.")
    parser.add_argument("filepath", type=str, help="Path to the input sequence file (.fasta or .fastq)")
    parser.add_argument("-f", "--format", type=str, choices=["fasta", "fastq"], default="fastq")
    parser.add_argument("-a", "--assembler", type=str, choices=["greedy", "debruijn", "both"], default="both")
    
    # Expose coverage threshold
    parser.add_argument("-c", "--min-cov", type=int, default=15, help="Minimum k-mer coverage threshold")
    parser.add_argument("-k", "--kmer", type=int, default=31, help="K-mer size (Max: 31)")
    
    parser.add_argument("-s", "--save", type=str, help="Path to save the assembled FASTA file", default="assembled_genome.fasta")
    # Fixed flag to standard format
    parser.add_argument("--log", action="store_true", help="Generate a detailed timestamped .log file")
    
    args = parser.parse_args()
    
    # Proactive Input Validation
    if args.kmer > 31:
        parser.error("K-mer size cannot exceed 31 due to 64-bit memory optimization limits.")
    if not os.path.exists(args.filepath):
        parser.error(f"Input file not found: {args.filepath}")
    
    run_pipeline(
        filepath=args.filepath, file_format=args.format,
        min_overlap=args.overlap, max_mismatches=args.mismatches,
        kmer_size=args.kmer, assembler_choice=args.assembler,
        enable_log=args.log,
        output_path=args.save  # Pass the save path here
    )

if __name__ == "__main__":
    main()