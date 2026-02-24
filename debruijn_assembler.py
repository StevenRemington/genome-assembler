import logging
from collections import defaultdict
from tqdm import tqdm
import numpy as np

class DNAEncoder:
    """Handles bit-packing DNA sequences into highly memory-efficient integers."""
    
    MAP_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0} # Handle 'N' safely
    MAP_TO_CHAR = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

    @staticmethod
    def encode(sequence):
        """Converts a string of DNA into a 2-bit encoded integer."""
        val = 0
        for char in sequence:
            val = (val << 2) | DNAEncoder.MAP_TO_INT.get(char.upper(), 0)
        return val

    @staticmethod
    def decode(val, length):
        """Converts an integer back into a DNA string of the specified length."""
        chars = []
        for _ in range(length):
            chars.append(DNAEncoder.MAP_TO_CHAR[val & 3])
            val >>= 2
        return "".join(reversed(chars))


class DeBruijnAssembler:
    """Extreme memory-optimized De Bruijn Graph using Implicit Bitmasks."""
    
    def __init__(self, k=31, min_coverage=15, num_anchors=10):
        self.k = k
        self.min_coverage = min_coverage
        self.num_anchors = num_anchors
        self.kmer_mask = (1 << (2 * self.k)) - 1
        self.node_mask = (1 << (2 * (self.k - 1))) - 1
        self.logger = logging.getLogger("GenomeAssembler.DeBruijn")

    def _build_graph(self, read_stream):
        """Builds a vectorized graph with STRICT k-mer edge filtering."""
        self.logger.info(f"Extracting k-mers (Threshold: {self.min_coverage})...")
        
        full_kmers_list = []
        short_reads_discarded = 0

        for read in tqdm(read_stream, desc="1/4 Extracting K-mers", unit="read", leave=False):
            if len(read) < self.k:
                short_reads_discarded += 1
                continue
                
            # Process strictly the forward read to save 50% RAM
            current_kmer_val = DNAEncoder.encode(read[:self.k])
            full_kmers_list.append(current_kmer_val)
            for i in range(self.k, len(read)):
                next_char_val = DNAEncoder.MAP_TO_INT.get(read[i].upper(), 0)
                current_kmer_val = ((current_kmer_val << 2) & self.kmer_mask) | next_char_val
                full_kmers_list.append(current_kmer_val)

        self.logger.info("Sorting raw k-mers...")
        full_kmers_arr = np.array(full_kmers_list, dtype=np.uint64)
        del full_kmers_list # Free up heavy Python list memory immediately
        
        full_kmers_arr.sort()
        
        self.logger.info("Compacting and filtering TRUE noise...")
        changes = np.concatenate(([True], full_kmers_arr[1:] != full_kmers_arr[:-1]))
        unique_indices = np.where(changes)[0]
        counts = np.diff(np.concatenate((unique_indices, [len(full_kmers_arr)])))
        
        keep_mask = counts >= self.min_coverage
        valid_full_kmers = full_kmers_arr[unique_indices][keep_mask]
        del full_kmers_arr
        
        self.logger.info("Extracting graph edges...")
        prefixes = valid_full_kmers >> 2
        chars = (valid_full_kmers & 3).astype(np.uint8)
        edges = (np.array(1, dtype=np.uint8) << chars)
        del valid_full_kmers
        
        prefix_changes = np.concatenate(([True], prefixes[1:] != prefixes[:-1]))
        prefix_unique_indices = np.where(prefix_changes)[0]
        
        final_kmers = prefixes[prefix_unique_indices]
        final_edges = np.bitwise_or.reduceat(edges, prefix_unique_indices)
        
        total_edges = int(np.sum(np.unpackbits(final_edges).reshape(-1, 8).sum(axis=1)))
        
        self.logger.info(f"Graph built: {len(final_kmers):,} nodes, {total_edges:,} edges.")
        return final_kmers, final_edges, total_edges, short_reads_discarded

    def _extract_contigs(self, kmers_arr, edges_arr):
        """Vectorized O(1) Contig Extractor using pre-computed memory pointers."""
        n_nodes = len(kmers_arr)
        self.logger.info(f"Step 2/4: Pre-computing O(1) transition pointers for {n_nodes:,} nodes...")
        
        valid_edges = (edges_arr[:, None] & (1 << np.arange(4, dtype=np.uint8))) > 0
        source_idx, char_idx = np.where(valid_edges)
        
        u_vals = kmers_arr[source_idx]
        v_vals = ((u_vals << 2) & self.node_mask) | char_idx.astype(np.uint64)
        
        dest_idx = np.searchsorted(kmers_arr, v_vals)
        
        valid_dest_mask = dest_idx < n_nodes
        valid_dest_mask[valid_dest_mask] = kmers_arr[dest_idx[valid_dest_mask]] == v_vals[valid_dest_mask]
        
        # In-degree computation to prevent chimeras
        in_degrees = np.zeros(n_nodes, dtype=np.int32)
        np.add.at(in_degrees, dest_idx[valid_dest_mask], 1)

        pointers = np.full((n_nodes, 4), -1, dtype=np.int32)
        pointers[source_idx[valid_dest_mask], char_idx[valid_dest_mask]] = dest_idx[valid_dest_mask]
        
        del valid_edges, source_idx, char_idx, u_vals, v_vals, dest_idx, valid_dest_mask
        
        self.logger.info("Step 3/4: Extracting Linear Contigs (Lightning Speed)...")
        all_paths = []
        working_edges = edges_arr.copy()
        active_nodes = np.where(working_edges > 0)[0]

        with tqdm(total=len(active_nodes), desc="Traversing Graph") as pbar:
            for start_idx in active_nodes:
                if working_edges[start_idx] == 0:
                    continue
                
                curr_idx = start_idx
                path = [kmers_arr[curr_idx]]
                
                while True:
                    out_mask = int(working_edges[curr_idx])
                    
                    if bin(out_mask).count('1') != 1:
                        break
                    
                    char_val = 0
                    for bit in range(4):
                        if out_mask & (1 << bit):
                            char_val = bit
                            break
                            
                    working_edges[curr_idx] &= ~(1 << char_val) & 0xFF
                    pbar.update(1)
                    
                    next_idx = pointers[curr_idx, char_val]
                    
                    if next_idx == -1:
                        break 
                        
                    # Stop if multiple paths merge into this node
                    if in_degrees[next_idx] > 1:
                        break

                    curr_idx = next_idx
                    path.append(kmers_arr[curr_idx])
                    
                if len(path) >= self.k:
                    all_paths.append(path)
                    
        return all_paths

    def _deduplicate_contigs(self, contigs):
        """
        Removes double-stranded duplicates and redundant sub-fragments.
        Stores ALL k-mers of accepted contigs to ensure flawless RC detection.
        """
        self.logger.info("Deduplicating double-stranded and overlapping contigs via Dense K-mer Anchoring...")
        trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
        unique_contigs = []
        seen_kmers = set()
        
        for seq in tqdm(contigs, desc="Deduplicating", unit="contig", leave=False):
            if len(seq) < self.k:
                continue
                
            stride = max(1, (len(seq) - self.k) // 100) 
            test_kmers = [seq[i : i+self.k] for i in range(0, len(seq) - self.k + 1, stride)]
            
            if not test_kmers:
                continue
                
            match_count = sum(1 for kmer in test_kmers if kmer in seen_kmers)
                    
            if match_count > (len(test_kmers) * 0.5):
                continue 
                
            unique_contigs.append(seq)
            
            for i in range(len(seq) - self.k + 1):
                fwd_kmer = seq[i:i+self.k]
                rc_kmer = fwd_kmer.translate(trans)[::-1]
                seen_kmers.add(fwd_kmer)
                seen_kmers.add(rc_kmer)
                
        removed = len(contigs) - len(unique_contigs)
        self.logger.info(f"Deduplication removed {removed} redundant/mirrored contigs.")
        return unique_contigs

    def assemble(self, read_stream):
        """End-to-end execution of the memory-optimized De Bruijn pipeline."""
        self.logger.info(f"Starting Vectorized De Bruijn Assembly (k={self.k}).")

        kmers_arr, edges_arr, total_edges, discarded = self._build_graph(read_stream)
        
        if len(kmers_arr) == 0:
            self.logger.error("Graph is empty. Assembly failed.")
            return [], {}
            
        contig_paths = self._extract_contigs(kmers_arr, edges_arr)
        
        self.logger.info("Step 4/4: Decoding contigs to DNA...")
        final_contigs = []
        for path in contig_paths:
            seq = DNAEncoder.decode(path[0], self.k - 1)
            seq += "".join(DNAEncoder.MAP_TO_CHAR[int(node_int) & 3] for node_int in path[1:])
            final_contigs.append(seq)
        
        final_contigs.sort(key=len, reverse=True)
        
        # PROPERLY SCOPED CALL TO THE INSTANCE METHOD
        final_contigs = self._deduplicate_contigs(final_contigs)
        
        total_bases = sum(len(c) for c in final_contigs)
        longest_contig = len(final_contigs[0]) if final_contigs else 0
            
        stats = {
            "Algorithm": "De Bruijn Graph (Linear Contig Extraction)",
            "Settings": f"k-mer size: {self.k}, min_cov: {self.min_coverage}",
            "TOTAL GENOME SIZE (BASES)": f"{total_bases:,}",
            "Longest Contig": f"{longest_contig:,}",
            "Total Contigs": len(final_contigs),
            "Total Filtered Edges": total_edges,
            "Data Discarded": f"{discarded} short reads"
        }
        
        self.logger.info(f"Assembly completed. Total Bases: {total_bases:,}. Contigs: {len(final_contigs)}.")
        return final_contigs, stats