import logging
from collections import defaultdict
from tqdm import tqdm
import numpy as np

# Import the metrics utility
from metrics import calculate_assembly_metrics

class DNAEncoder:
    """
    Handles bit-packing DNA sequences into highly memory-efficient integers.
    
    Why: Standard Python strings take up roughly 50 bytes plus 1 byte per character. 
    By converting 'A','C','G','T' into 2-bit representations (00, 01, 10, 11), a 31-mer 
    fits perfectly inside a single 64-bit integer, reducing the graph's memory footprint 
    by over 80% and enabling C-level hardware speed during processing.
    """
    
    MAP_TO_INT = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0} 
    MAP_TO_CHAR = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}

    @staticmethod
    def encode(sequence):
        """
        Converts a string of DNA into a 2-bit encoded integer.
        
        Args:
            sequence (str): The raw string sequence of DNA (e.g., "ATGC").
            
        Returns:
            int: The 64-bit integer representation of the sequence.
        """
        val = 0
        for char in sequence:
            # Shift the existing bits left by 2 spaces to make room, 
            # then append the 2-bit value of the current character using bitwise OR.
            val = (val << 2) | DNAEncoder.MAP_TO_INT.get(char.upper(), 0)
        return val

    @staticmethod
    def decode(val, length):
        """
        Converts an integer back into a DNA string of the specified length.
        
        Args:
            val (int): The 2-bit encoded integer.
            length (int): The number of characters to decode (usually k or k-1).
            
        Returns:
            str: The decoded string of DNA bases.
        """
        chars = []
        for _ in range(length):
            # Extract the last 2 bits using bitwise AND (& 3 extracts 00, 01, 10, or 11)
            chars.append(DNAEncoder.MAP_TO_CHAR[val & 3])
            # Shift right by 2 to process the next character in the next loop iteration
            val >>= 2
        # Because we extract from the right (end of the sequence), we must reverse it
        return "".join(reversed(chars))


class DeBruijnAssembler:
    """
    Extreme memory-optimized De Bruijn Graph using Implicit Bitmasks.
    
    Why: Traditional graph construction in Python uses dictionaries (`dict`), which incur 
    massive overhead. This class strictly uses NumPy arrays (`np.uint64`), bitwise shifts, 
    and implicit relationships to bypass Python's slowness entirely.
    """
    
    def __init__(self, k=31, min_coverage=20, num_anchors=10):
        """
        Initializes the Assembler.
        
        Args:
            k (int): The k-mer size. Maximum is 31 because 31*2 = 62 bits, which safely 
                     fits inside a 64-bit unsigned integer (np.uint64).
            min_coverage (int): The minimum number of times a k-mer must appear to be trusted.
            num_anchors (int): Unused here natively, but represents bounds for deduplication checks.
        """
        self.k = k
        self.min_coverage = min_coverage
        self.num_anchors = num_anchors
        
        # Bitmasks used to rapidly drop prefix/suffix characters without string slicing
        self.kmer_mask = (1 << (2 * self.k)) - 1
        self.node_mask = (1 << (2 * (self.k - 1))) - 1
        self.logger = logging.getLogger("GenomeAssembler.DeBruijn")

    def _build_graph(self, read_stream):
        """
        Builds a vectorized graph with STRICT k-mer edge filtering.
        
        Why: Instead of adding nodes one-by-one, we extract every single k-mer from the data, 
        sort the giant array, and use vectorization to instantly find duplicates and build edges.
        
        Args:
            read_stream (iterator): An iterable of string DNA reads from the source file.
            
        Returns:
            tuple: (final_kmers (np.array), final_edges (np.array), total_edges (int), discarded (int))
        """
        self.logger.info(f"Extracting k-mers (Threshold: {self.min_coverage})...")
        
        full_kmers_list = []
        short_reads_discarded = 0

        for read in tqdm(read_stream, desc="1/4 Extracting K-mers", unit="read", leave=False):
            if len(read) < self.k:
                short_reads_discarded += 1
                continue
                
            # We strictly process the forward read here. 
            # Why: Processing both strands during construction doubles RAM usage. 
            # We save 50% RAM by assembling the forward graph and deduplicating at the end.
            current_kmer_val = DNAEncoder.encode(read[:self.k])
            full_kmers_list.append(current_kmer_val)
            for i in range(self.k, len(read)):
                # "Rolling Window" extraction: Shift left to drop the first character, 
                # apply mask to prevent overflow, and append the next character.
                next_char_val = DNAEncoder.MAP_TO_INT.get(read[i].upper(), 0)
                current_kmer_val = ((current_kmer_val << 2) & self.kmer_mask) | next_char_val
                full_kmers_list.append(current_kmer_val)

        self.logger.info("Sorting raw k-mers...")
        full_kmers_arr = np.array(full_kmers_list, dtype=np.uint64)
        del full_kmers_list # Instantly free up standard Python list overhead
        
        # Sorting groups identical k-mers together, enabling fast C-level counting
        full_kmers_arr.sort()
        
        self.logger.info("Compacting and filtering TRUE noise...")
        # Boolean mask indicating where the value changes (i.e., a new unique k-mer starts)
        changes = np.concatenate(([True], full_kmers_arr[1:] != full_kmers_arr[:-1]))
        unique_indices = np.where(changes)[0]
        # Calculate coverage by subtracting the start index of one k-mer from the start of the next
        counts = np.diff(np.concatenate((unique_indices, [len(full_kmers_arr)])))
        
        # Drop anything below min_coverage to completely eliminate random sequencing errors
        keep_mask = counts >= self.min_coverage
        valid_full_kmers = full_kmers_arr[unique_indices][keep_mask]
        del full_kmers_arr
        
        self.logger.info("Extracting graph edges...")
        # The node (prefix) is the first k-1 characters (shift right by 2 to drop the last char)
        prefixes = valid_full_kmers >> 2
        # The edge (suffix character) is the last character
        chars = (valid_full_kmers & 3).astype(np.uint8)
        
        # Convert edge character (0-3) into a bitwise position (1, 2, 4, or 8)
        # Why: This allows us to track multiple outgoing edges from the same node using a single byte (e.g. 0101 means edges C and T exist)
        edges = (np.array(1, dtype=np.uint8) << chars)
        del valid_full_kmers
        
        # Because we shifted right, adjacent elements might now share the same prefix. Group them.
        prefix_changes = np.concatenate(([True], prefixes[1:] != prefixes[:-1]))
        prefix_unique_indices = np.where(prefix_changes)[0]
        
        final_kmers = prefixes[prefix_unique_indices]
        # Bitwise OR collapses multiple edges belonging to the same node into a single byte
        final_edges = np.bitwise_or.reduceat(edges, prefix_unique_indices)
        
        total_edges = int(np.sum(np.unpackbits(final_edges).reshape(-1, 8).sum(axis=1)))
        
        self.logger.info(f"Graph built: {len(final_kmers):,} nodes, {total_edges:,} edges.")
        return final_kmers, final_edges, total_edges, short_reads_discarded

    def _build_transition_pointers(self, kmers_arr, edges_arr):
        """
        Pre-computes O(1) transition pointers and node in-degrees.
        
        Why: Extracts the dense NumPy matrix math out of the traversal logic. 
        This isolated method can now be easily unit-tested for edge-case accuracy.
        """
        n_nodes = len(kmers_arr)
        self.logger.info(f"Step 2/4: Pre-computing O(1) transition pointers for {n_nodes:,} nodes...")
        
        # Identify which specific outgoing edges (0-3) are active for each node
        valid_edges = (edges_arr[:, None] & (1 << np.arange(4, dtype=np.uint8))) > 0
        source_idx, char_idx = np.where(valid_edges)
        
        # Reconstruct the destination node (suffix)
        u_vals = kmers_arr[source_idx]
        v_vals = ((u_vals << 2) & self.node_mask) | char_idx.astype(np.uint64)
        
        # Find exactly where that destination node lives in the kmers_arr array
        dest_idx = np.searchsorted(kmers_arr, v_vals)
        
        # Ensure the destination actually exists
        valid_dest_mask = dest_idx < n_nodes
        valid_dest_mask[valid_dest_mask] = kmers_arr[dest_idx[valid_dest_mask]] == v_vals[valid_dest_mask]
        
        # Track In-Degrees to prevent Chimeras
        in_degrees = np.zeros(n_nodes, dtype=np.int32)
        np.add.at(in_degrees, dest_idx[valid_dest_mask], 1)

        # Build the pointer matrix [Node Index, Edge Character] -> [Destination Index]
        pointers = np.full((n_nodes, 4), -1, dtype=np.int32)
        pointers[source_idx[valid_dest_mask], char_idx[valid_dest_mask]] = dest_idx[valid_dest_mask]
        
        return pointers, in_degrees

    def _try_pop_bubble(self, start_idx, out_mask, pointers, working_edges):
        """
        Graph Cleaning Part 1: Bubble Popping.
        
        Attempts to find a converging path within k steps (a sequencing error bubble).
        Extracted to a class method to allow for isolated unit testing.
        """
        branches = [bit for bit in range(4) if out_mask & (1 << bit)]
        
        def trace_path(start_char):
            p = []
            curr = pointers[start_idx, start_char]
            if curr == -1: return p
            p.append(curr)
            
            # Search up to k+2 steps ahead to account for complex but tiny variants
            for _ in range(self.k + 2): 
                o_mask = int(working_edges[curr])
                if bin(o_mask).count('1') != 1:
                    break
                next_char = 0
                for b in range(4):
                    if o_mask & (1 << b):
                        next_char = b
                        break
                curr = pointers[curr, next_char]
                if curr == -1: break
                p.append(curr)
            return p
            
        p1 = trace_path(branches[0])
        p2 = trace_path(branches[1])
        
        if not p1 or not p2:
            return None, None, None, None, None
            
        # Check if the paths intersect (converge)
        set_p1 = set(p1)
        for i, node in enumerate(p2):
            if node in set_p1:
                idx1 = p1.index(node)
                idx2 = i
                return node, p1[:idx1], p2[:idx2], branches[0], branches[1]
                
        return None, None, None, None, None

    def _extract_contigs(self, kmers_arr, edges_arr):
        """
        Vectorized O(1) Contig Extractor.
        
        Why: Now strictly responsible for orchestrating the traversal loop. 
        It delegates pointer math and bubble logic to their respective helpers.
        """
        # 1. Delegate pointer building
        pointers, in_degrees = self._build_transition_pointers(kmers_arr, edges_arr)
        
        self.logger.info("Step 3/4: Extracting Linear Contigs (Lightning Speed)...")
        all_paths = []
        tips_removed = 0 
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
                    out_count = bin(out_mask).count('1')
                    
                    if out_count == 0:
                        break # End of the line
                    
                    if out_count == 1:
                        # STANDARD TRAVERSAL
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
                            
                        # Stop to prevent Chimeras
                        if in_degrees[next_idx] > 1:
                            break

                        curr_idx = next_idx
                        path.append(kmers_arr[curr_idx])
                        
                    elif out_count == 2:
                        # OUTWARD BRANCH: Delegate bubble popping
                        conv_node, p1_nodes, p2_nodes, b1, b2 = self._try_pop_bubble(
                            curr_idx, out_mask, pointers, working_edges
                        )
                        
                        if conv_node is not None:
                            # It was an error bubble! Safely resolve it.
                            working_edges[curr_idx] &= ~(1 << b1) & 0xFF
                            working_edges[curr_idx] &= ~(1 << b2) & 0xFF
                            pbar.update(1)
                            
                            for n in p1_nodes + p2_nodes:
                                working_edges[n] = 0
                                
                            for n in p1_nodes:
                                path.append(kmers_arr[n])
                                
                            curr_idx = conv_node
                            path.append(kmers_arr[curr_idx])
                        else:
                            # True biological repeat; stop safely.
                            break
                    else:
                        break # Out-degree > 2 (Complex repetitive region, stop safely)
                    
                # Graph Cleaning Part 2: Tip Removal
                if len(path) >= 2 * self.k:
                    all_paths.append(path)
                elif len(path) >= self.k:
                    tips_removed += 1
                    
        self.logger.info(f"Graph Cleaning: Pruned {tips_removed:,} dead-end tips.")
        return all_paths, tips_removed
    def _deduplicate_contigs(self, contigs):
        """
        Removes double-stranded duplicates and redundant sub-fragments.
        
        Why: We saved RAM by only processing the forward strand during construction. This naturally 
        creates a "mirror" assembly where the output contains both the forward genome AND its 
        reverse-complement twin. We use dense K-mer anchoring to detect and discard the twin, 
        yielding the true haploid genome.
        
        Args:
            contigs (list[str]): The fully extracted and decoded string DNA contigs.
            
        Returns:
            list[str]: The strictly unique, deduplicated contigs.
        """
        self.logger.info("Deduplicating double-stranded and overlapping contigs via Dense K-mer Anchoring...")
        trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
        unique_contigs = []
        seen_kmers = set()
        
        for seq in tqdm(contigs, desc="Deduplicating", unit="contig", leave=False):
            if len(seq) < self.k:
                continue
                
            # Performance Optimization: Instead of checking every single k-mer in a 50,000 base contig,
            # we sample test k-mers across its length (striding up to 100 times).
            stride = max(1, (len(seq) - self.k) // 100) 
            test_kmers = [seq[i : i+self.k] for i in range(0, len(seq) - self.k + 1, stride)]
            
            if not test_kmers:
                continue
                
            match_count = sum(1 for kmer in test_kmers if kmer in seen_kmers)
                    
            # If the majority of the sampled anchors already exist in the database, 
            # this contig is redundant or a mirror. Toss it.
            if match_count > (len(test_kmers) * 0.5):
                continue 
                
            unique_contigs.append(seq)
            
            # If accepted, we must add ALL of its k-mers (and their Reverse Complements)
            # to the database so future fragments know they have already been covered.
            for i in range(len(seq) - self.k + 1):
                fwd_kmer = seq[i:i+self.k]
                rc_kmer = fwd_kmer.translate(trans)[::-1]
                seen_kmers.add(fwd_kmer)
                seen_kmers.add(rc_kmer)
                
        removed = len(contigs) - len(unique_contigs)
        self.logger.info(f"Deduplication removed {removed} redundant/mirrored contigs.")
        return unique_contigs

    def assemble(self, read_stream):
        """
        End-to-end execution of the memory-optimized De Bruijn pipeline.
        
        Orchestrates the four core phases: Extract -> Traversal -> Decode -> Deduplicate.
        
        Args:
            read_stream (iterator): Sequence generator from the source FASTQ/FASTA.
            
        Returns:
            tuple: (final_contigs (list[str]), stats (dict))
        """
        self.logger.info(f"Starting Vectorized De Bruijn Assembly (k={self.k}).")

        kmers_arr, edges_arr, total_edges, discarded = self._build_graph(read_stream)
        
        if len(kmers_arr) == 0:
            self.logger.error("Graph is empty. Assembly failed.")
            return [], {}
            
        contig_paths, tips_removed = self._extract_contigs(kmers_arr, edges_arr)
        
        self.logger.info("Step 4/4: Decoding contigs to DNA...")
        final_contigs = []
        for path in contig_paths:
            # Decode the first node fully (length k-1)
            seq = DNAEncoder.decode(path[0], self.k - 1)
            # For every subsequent node in the path, only decode and append the final character 
            # (since k-1 characters already overlap by definition)
            seq += "".join(DNAEncoder.MAP_TO_CHAR[int(node_int) & 3] for node_int in path[1:])
            final_contigs.append(seq)
        
        # Sort largest to smallest so Deduplication evaluates massive anchor contigs first
        final_contigs.sort(key=len, reverse=True)
        final_contigs = self._deduplicate_contigs(final_contigs)
        
        total_bases = sum(len(c) for c in final_contigs)
        longest_contig = len(final_contigs[0]) if final_contigs else 0
        
        # Generate biological quality metrics
        quality_metrics = calculate_assembly_metrics(final_contigs)
            
        stats = {
            "Algorithm": "De Bruijn Graph (Linear Contig Extraction)",
            "Settings": f"k-mer size: {self.k}, min_cov: {self.min_coverage}",
            "TOTAL GENOME SIZE (BASES)": f"{total_bases:,}",
            "Longest Contig": f"{longest_contig:,}",
            "Total Contigs": len(final_contigs),
            "Tips Pruned": tips_removed, 
            "Total Filtered Edges": total_edges,
            "Data Discarded": f"{discarded} short reads",
            "N50": quality_metrics["N50"],
            "L50": quality_metrics["L50"],
            "GC Content": quality_metrics["GC Content"]
        }
        
        self.logger.info(f"Assembly completed. Total Bases: {total_bases:,}. Contigs: {len(final_contigs)}.")
        return final_contigs, stats