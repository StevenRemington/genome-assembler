import logging
from collections import defaultdict
from tqdm import tqdm
import numpy as np

logger = logging.getLogger("GenomeAssembler.DeBruijn")

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

def deduplicate_contigs(contigs, k=31):
    """
    Removes reverse-complement duplicates and redundant sub-fragments 
    using K-mer anchoring. Assumes 'contigs' is sorted longest-first.
    """
    logger = logging.getLogger("GenomeAssembler.DeBruijn")
    logger.info("Deduplicating double-stranded and overlapping contigs via K-mer Anchoring...")
    
    trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    unique_contigs = []
    seen_kmers = set()
    
    # We sample up to 10 evenly spaced "anchors" per contig
    num_anchors = 10 
    
    for seq in contigs:
        if len(seq) < k:
            continue
            
        # Extract anchors from the current contig
        step = max(1, (len(seq) - k) // num_anchors)
        anchors = [seq[i : i+k] for i in range(0, len(seq) - k + 1, step)]
        
        # Check if this contig's anchors exist in our 'seen' database
        is_duplicate = False
        for anchor in anchors:
            rc_anchor = anchor.translate(trans)[::-1]
            if anchor in seen_kmers or rc_anchor in seen_kmers:
                is_duplicate = True
                break # We found a match, it's a duplicate!
                
        if not is_duplicate:
            unique_contigs.append(seq)
            # Add all anchors of this NEW unique contig to the database
            for anchor in anchors:
                rc_anchor = anchor.translate(trans)[::-1]
                seen_kmers.add(anchor)
                seen_kmers.add(rc_anchor)
                
    removed = len(contigs) - len(unique_contigs)
    logger.info(f"Deduplication removed {removed} redundant/mirrored contigs.")
    return unique_contigsOka


class DeBruijnAssembler:
    """Extreme memory-optimized De Bruijn Graph using Implicit Bitmasks."""
    
    def __init__(self, k=31, min_coverage=15, num_anchors=10):
        self.k = k
        self.min_coverage = min_coverage
        self.num_anchors = num_anchors
        self.kmer_mask = (1 << (2 * self.k)) - 1
        self.node_mask = (1 << (2 * (self.k - 1))) - 1
        self.logger = logging.getLogger("GenomeAssembler.DeBruijn")

    def _build_graph(self, read_stream, min_coverage=3):
        """Builds a vectorized graph with STRICT k-mer edge filtering."""
        logger.info(f"Extracting k-mers (Threshold: {min_coverage})...")
        
        full_kmers_list = []
        short_reads_discarded = 0

        for read in tqdm(read_stream, desc="1/4 Extracting K-mers", unit="read", leave=False):
            if len(read) < self.k:
                short_reads_discarded += 1
                continue
                
            # Grab the very first full k-mer
            current_kmer_val = DNAEncoder.encode(read[:self.k])
            full_kmers_list.append(current_kmer_val)
            
            # Slide the window for the rest of the read
            for i in range(self.k, len(read)):
                next_char_val = DNAEncoder.MAP_TO_INT.get(read[i].upper(), 0)
                current_kmer_val = ((current_kmer_val << 2) & self.kmer_mask) | next_char_val
                full_kmers_list.append(current_kmer_val)

        logger.info("Sorting raw k-mers...")
        full_kmers_arr = np.array(full_kmers_list, dtype=np.uint64)
        del full_kmers_list
        
        # Sort the FULL k-mers (This is the critical change)
        full_kmers_arr.sort()
        
        logger.info("Compacting and filtering TRUE noise...")
        # Group identical full k-mers
        changes = np.concatenate(([True], full_kmers_arr[1:] != full_kmers_arr[:-1]))
        unique_indices = np.where(changes)[0]
        counts = np.diff(np.concatenate((unique_indices, [len(full_kmers_arr)])))
        
        # Apply coverage threshold to the FULL k-mer
        keep_mask = counts >= min_coverage
        valid_full_kmers = full_kmers_arr[unique_indices][keep_mask]
        del full_kmers_arr
        
        logger.info("Extracting graph edges...")
        # Now we extract the node (prefix) and the edge (last char) from the SURVIVING k-mers
        prefixes = valid_full_kmers >> 2
        chars = (valid_full_kmers & 3).astype(np.uint8)
        edges = (np.array(1, dtype=np.uint8) << chars)
        del valid_full_kmers
        
        # Because we shifted right by 2, the 'prefixes' array is ALREADY SORTED! 
        # We can group them instantly without running another sort.
        prefix_changes = np.concatenate(([True], prefixes[1:] != prefixes[:-1]))
        prefix_unique_indices = np.where(prefix_changes)[0]
        
        final_kmers = prefixes[prefix_unique_indices]
        final_edges = np.bitwise_or.reduceat(edges, prefix_unique_indices)
        
        total_edges = int(np.sum(np.unpackbits(final_edges).reshape(-1, 8).sum(axis=1)))
        
        logger.info(f"Graph built: {len(final_kmers):,} nodes, {total_edges:,} edges.")
        return final_kmers, final_edges, total_edges, short_reads_discarded
    def _extract_contigs(self, kmers_arr, edges_arr):
        """Vectorized O(1) Contig Extractor using pre-computed memory pointers."""
        n_nodes = len(kmers_arr)
        logger.info(f"Step 2/4: Pre-computing O(1) transition pointers for {n_nodes:,} nodes...")
        
        # 1. FIND ALL VALID EDGES IN ONE C-LEVEL OPERATION
        # Create a boolean mask of shape (n_nodes, 4) showing exactly where edges exist
        valid_edges = (edges_arr[:, None] & (1 << np.arange(4, dtype=np.uint8))) > 0
        
        # Get the row (source node index) and column (character index 0-3) for every edge
        source_idx, char_idx = np.where(valid_edges)
        
        # Calculate the integer value of the destination node for every edge simultaneously
        u_vals = kmers_arr[source_idx]
        v_vals = ((u_vals << 2) & self.node_mask) | char_idx.astype(np.uint64)
        
        # Perform ONE massive binary search for all millions of edges at once
        dest_idx = np.searchsorted(kmers_arr, v_vals)
        
        # Filter out dead-ends (where the destination node doesn't exist in our filtered graph)
        valid_dest_mask = dest_idx < n_nodes
        valid_dest_mask[valid_dest_mask] = kmers_arr[dest_idx[valid_dest_mask]] == v_vals[valid_dest_mask]
        
        # 2. BUILD THE O(1) POINTER ARRAY
        # pointers[node_index, char] = next_node_index
        # We use -1 to represent a dead end.
        pointers = np.full((n_nodes, 4), -1, dtype=np.int32)
        pointers[source_idx[valid_dest_mask], char_idx[valid_dest_mask]] = dest_idx[valid_dest_mask]
        
        # Free up heavy memory before traversal
        del valid_edges, source_idx, char_idx, u_vals, v_vals, dest_idx, valid_dest_mask
        
        # 3. FAST TRAVERSAL
        logger.info("Step 3/4: Extracting Linear Contigs (Lightning Speed)...")
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
                    
                    # Stop at a branch (more than 1 edge)
                    if bin(out_mask).count('1') != 1:
                        break
                    
                    # Find the single edge character
                    char_val = 0
                    for bit in range(4):
                        if out_mask & (1 << bit):
                            char_val = bit
                            break
                            
                    # Consume the edge
                    working_edges[curr_idx] &= ~(1 << char_val) & 0xFF
                    pbar.update(1)
                    
                    # O(1) POINTER LOOKUP! No more binary searches in the loop.
                    next_idx = pointers[curr_idx, char_val]
                    
                    if next_idx == -1:
                        break # Dead end
                        
                    curr_idx = next_idx
                    path.append(kmers_arr[curr_idx])
                    
                if len(path) >= self.k:
                    all_paths.append(path)
                    
        return all_paths

    def assemble(self, read_stream):
        """End-to-end execution of the memory-optimized De Bruijn pipeline."""
        logger.info(f"Starting Vectorized De Bruijn Assembly (k={self.k}).")

        kmers_arr, edges_arr, total_edges, discarded = self._build_graph(read_stream, min_coverage=15)
        
        if len(kmers_arr) == 0:
            logger.error("Graph is empty. Assembly failed.")
            return [], {}
            
        # Instead of one Eulerian path, we extract all linear contigs
        contig_paths = self._extract_contigs(kmers_arr, edges_arr)
        
        logger.info("Step 4/4: Decoding contigs to DNA...")
        final_contigs = []
        for path in contig_paths:
            seq = DNAEncoder.decode(path[0], self.k - 1)
            seq += "".join(DNAEncoder.MAP_TO_CHAR[int(node_int) & 3] for node_int in path[1:])
            final_contigs.append(seq)
        
        # Sort contigs by length (longest first)
        final_contigs.sort(key=len, reverse=True)
        
        # Deduplication step
        final_contigs = deduplicate_contigs(final_contigs)
        
        total_bases = sum(len(c) for c in final_contigs)
        longest_contig = len(final_contigs[0]) if final_contigs else 0
            
        stats = {
            "Algorithm": "De Bruijn Graph (Linear Contig Extraction)",
            "Settings": f"k-mer size: {self.k}, min_cov: 15",
            "TOTAL GENOME SIZE (BASES)": f"{total_bases:,}",
            "Longest Contig": f"{longest_contig:,}",
            "Total Contigs": len(final_contigs),
            "Total Filtered Edges": total_edges,
            "Data Discarded": f"{discarded} short reads"
        }
        
        logger.info(f"Assembly completed. Total Bases: {total_bases:,}. Contigs: {len(final_contigs)}.")
        return final_contigs, stats