# assembler.py

import logging
from tqdm import tqdm

# Connect to the main logger
logger = logging.getLogger("GenomeAssembler.Greedy")

class GreedyAssembler:
    def __init__(self, min_overlap=5, max_mismatches=1):
        self.min_overlap = min_overlap
        self.max_mismatches = max_mismatches

    def _calculate_overlap(self, seq1, seq2):
        max_possible_overlap = min(len(seq1), len(seq2))
        for k in range(max_possible_overlap, self.min_overlap - 1, -1):
            suffix = seq1[-k:]
            prefix = seq2[:k]
            mismatches = sum(1 for char1, char2 in zip(suffix, prefix) if char1 != char2)
            if mismatches <= self.max_mismatches:
                return k, mismatches
        return 0, 0

    def assemble(self, fragments):
        pool = fragments[:]
        initial_count = len(pool)
        original_fragments = set(fragments)
        
        logger.info(f"Starting Greedy Assembly with {initial_count} fragments.")
        logger.debug(f"Settings: Min Overlap={self.min_overlap}, Max Mismatches={self.max_mismatches}")
        
        total_merges_possible = initial_count - 1
        
        with tqdm(total=total_merges_possible, desc="Greedy Assembly", unit="merge", leave=False) as pbar:
            while len(pool) > 1:
                best_overlap = -1
                best_pair = (-1, -1)
                merged_seq = ""
                
                for i in range(len(pool)):
                    for j in range(len(pool)):
                        if i == j: continue
                        overlap_len, mismatches = self._calculate_overlap(pool[i], pool[j])
                        if overlap_len > best_overlap:
                            best_overlap = overlap_len
                            best_pair = (i, j)
                            merged_seq = pool[i] + pool[j][overlap_len:]
                
                if best_overlap < self.min_overlap:
                    logger.debug("No valid overlaps remaining meeting the minimum threshold.")
                    break
                    
                i, j = best_pair
                logger.debug(f"Merging indices ({i}, {j}) | Overlap: {best_overlap} bases.")
                
                pool.pop(max(i, j))
                pool.pop(min(i, j))
                pool.append(merged_seq)
                
                pbar.update(1)
        
        # --- Generate Statistics ---
        # Assume the longest sequence is our target genome
        main_genome = max(pool, key=len) if pool else ""
        
        # Thrown out data: Any original fragments that were never merged into a larger contig
        thrown_out = [seq for seq in pool if seq in original_fragments and seq != main_genome]
        
        stats = {
            "Algorithm": "Greedy Overlap",
            "Settings": f"Min Overlap: {self.min_overlap}, Max Mismatches: {self.max_mismatches}",
            "Initial Reads": initial_count,
            "Total Merges": initial_count - len(pool),
            "Final Contigs": len(pool),
            "Main Genome Size": len(main_genome),
            "Data Thrown Out": len(thrown_out) # Number of unmerged reads
        }
        
        logger.info(f"Greedy Assembly completed. Main genome size: {stats['Main Genome Size']} bases.")
        return pool, stats