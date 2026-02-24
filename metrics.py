def calculate_assembly_metrics(contigs):
    """
    Calculates standard bioinformatics assembly metrics: N50, L50, and GC-Content.
    
    Args:
        contigs (list[str]): A list of assembled DNA sequences.
        
    Returns:
        dict: A dictionary of computed metrics.
    """
    if not contigs:
        return {"N50": 0, "L50": 0, "GC Content (%)": 0.0}

    # 1. Sort lengths in descending order
    lengths = sorted([len(c) for c in contigs], reverse=True)
    total_bases = sum(lengths)
    
    # 2. Calculate N50 and L50
    n50 = 0
    l50 = 0
    running_sum = 0
    half_length = total_bases / 2.0

    for i, length in enumerate(lengths):
        running_sum += length
        if running_sum >= half_length:
            n50 = length
            l50 = i + 1
            break

    # 3. Calculate GC Content
    # We use a generator expression to avoid concatenating all strings into RAM
    g_count = sum(c.upper().count('G') for c in contigs)
    c_count = sum(c.upper().count('C') for c in contigs)
    gc_content = ((g_count + c_count) / total_bases) * 100 if total_bases > 0 else 0.0

    return {
        "N50": f"{n50:,}",
        "L50": f"{l50:,}",
        "GC Content": f"{gc_content:.2f}%"
    }