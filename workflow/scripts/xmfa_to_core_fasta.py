input_file = snakemake.input[0]
output_file = snakemake.output[0]

def parse_xmfa(path):
    """Parse XMFA file into blocks with proper multi-line sequence handling"""
    with open(path) as f:
        blocks = []
        current_block = []
        current_header = None
        current_sequence = ""
        
        for line in f:
            line = line.strip()
            
            # Block separator
            if line.startswith('='):
                # Save current sequence if exists
                if current_header and current_sequence:
                    current_block.append(current_header)
                    current_block.append(current_sequence)
                
                # Save current block if it has content
                if current_block:
                    blocks.append(current_block)
                    current_block = []
                
                current_header = None
                current_sequence = ""
                
            # Header line
            elif line.startswith('>'):
                # Save previous sequence if exists
                if current_header and current_sequence:
                    current_block.append(current_header)
                    current_block.append(current_sequence)
                
                current_header = line
                current_sequence = ""
                
            # Sequence line (or continuation)
            elif line:
                current_sequence += line
        
        # Don't forget the last sequence and block
        if current_header and current_sequence:
            current_block.append(current_header)
            current_block.append(current_sequence)
        
        if current_block:
            blocks.append(current_block)
    
    return blocks

def extract_sequences_from_block(block):
    """Extract header-sequence pairs from a block"""
    sequences = []
    for i in range(0, len(block), 2):
        if i + 1 < len(block):
            header = block[i]
            sequence = block[i + 1]
            sequences.append((header, sequence))
    return sequences

def all_sequences_identical(sequences):
    """Check if all sequences are 100% identical (case-insensitive)"""
    if not sequences:
        return False
    
    first_seq = sequences[0][1].upper()
    return all(seq[1].upper() == first_seq for seq in sequences)

def find_conserved_subregions(sequences, min_length=18):
    """Find conserved subregions of at least min_length within aligned sequences"""
    if not sequences:
        return []
    
    # Get the length of aligned sequences (should be same for all)
    seq_length = len(sequences[0][1])
    
    # Check that all sequences have the same length (aligned)
    if not all(len(seq[1]) == seq_length for seq in sequences):
        print("Warning: Sequences in block have different lengths, skipping subregion analysis")
        return []
    
    conserved_regions = []
    i = 0
    
    while i < seq_length:
        # Check if position i is conserved across all sequences
        ref_char = sequences[0][1][i].upper()
        is_conserved = all(seq[1][i].upper() == ref_char for seq in sequences)
        
        if is_conserved:
            # Start of a potentially conserved region
            start = i
            
            # Extend as long as positions remain conserved
            while i < seq_length:
                ref_char = sequences[0][1][i].upper()
                if all(seq[1][i].upper() == ref_char for seq in sequences):
                    i += 1
                else:
                    break
            
            # Check if the conserved region is long enough
            region_length = i - start
            if region_length >= min_length:
                # Extract the conserved subsequence
                conserved_subseq = sequences[0][1][start:i]
                conserved_regions.append((start, i, conserved_subseq, region_length))
        else:
            i += 1
    
    return conserved_regions

def main():
    blocks = parse_xmfa(input_file)
    print(f"Parsed {len(blocks)} blocks from XMFA")
    
    # Determine total number of genomes from the first complete block
    total_genomes = None
    for block in blocks:
        sequences = extract_sequences_from_block(block)
        if sequences:
            total_genomes = len(sequences)
            break
    
    if total_genomes is None:
        print("Error: Could not determine number of genomes")
        return
    
    print(f"Expected number of genomes: {total_genomes}")
    
    fully_conserved_blocks = 0
    subregion_blocks = 0
    total_conserved_length = 0
    output_count = 0
    
    with open(output_file, "w") as out:
        for block_idx, block in enumerate(blocks):
            # Extract sequences from this block
            sequences = extract_sequences_from_block(block)
            
            # Skip if not all genomes are present
            if len(sequences) != total_genomes:
                continue
            
            # Check if all sequences are identical (fully conserved block)
            if all_sequences_identical(sequences):
                fully_conserved_blocks += 1
                sequence = sequences[0][1]
                output_count += 1
                
                out.write(f">fully_conserved_block_{block_idx+1:04d}_length_{len(sequence)}\n")
                
                # Write sequence in 80-character lines
                for i in range(0, len(sequence), 80):
                    out.write(sequence[i:i+80] + "\n")
                
                total_conserved_length += len(sequence)
                print(f"Block {block_idx + 1}: Fully conserved, {len(sequence)} bp")
            
            else:
                # Look for conserved subregions within this partially conserved block
                conserved_subregions = find_conserved_subregions(sequences, min_length=18)
                
                if conserved_subregions:
                    subregion_blocks += 1
                    print(f"Block {block_idx + 1}: Found {len(conserved_subregions)} conserved subregions")
                    
                    for subregion_idx, (start, end, subseq, length) in enumerate(conserved_subregions):
                        output_count += 1
                        
                        out.write(f">subregion_block_{block_idx+1:04d}_{subregion_idx+1:02d}_pos_{start+1}-{end}_length_{length}\n")
                        
                        # Write sequence in 80-character lines
                        for i in range(0, len(subseq), 80):
                            out.write(subseq[i:i+80] + "\n")
                        
                        total_conserved_length += length
                        print(f"  Subregion {subregion_idx + 1}: positions {start+1}-{end}, {length} bp")
    
    print(f"\nSummary:")
    print(f"Fully conserved blocks: {fully_conserved_blocks}")
    print(f"Blocks with conserved subregions: {subregion_blocks}")
    print(f"Total conserved sequences extracted: {output_count}")
    print(f"Total conserved sequence length: {total_conserved_length} bp")
    print(f"Output written to: {output_file}")

if __name__ == "__main__":
    main()

# Run the main function
main()
