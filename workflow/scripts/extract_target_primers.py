"""
Extract 100% conserved 18nt target sequences with 60nt primer regions
Allows up to 6 mismatches in primer regions for subsequent PrimedRPA validation
"""

import sys
import re
import pandas as pd
from Bio.SeqUtils import gc_fraction

def parse_xmfa(xmfa_file):
    """Parse XMFA file and yield alignment blocks"""
    current_block = {}
    current_seq_name = None
    
    with open(xmfa_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                # Extract genome name: >1:start-end strand genome_name
                parts = line[1:].split()
                if len(parts) >= 2:
                    genome_name = ' '.join(parts[1:])
                    current_block[genome_name] = []
                    current_seq_name = genome_name
            
            elif line.startswith('='):
                # End of block
                if current_block:
                    processed_block = {
                        name: ''.join(seq_parts).upper()
                        for name, seq_parts in current_block.items()
                    }
                    yield processed_block
                    current_block = {}
                    current_seq_name = None
            
            elif line and not line.startswith('#') and current_seq_name:
                # Sequence data - keep gaps, remove spaces
                current_block[current_seq_name].append(line.replace(' ', ''))
        
        # Yield last block if exists
        if current_block:
            processed_block = {
                name: ''.join(seq_parts).upper()
                for name, seq_parts in current_block.items()
            }
            yield processed_block

def get_consensus(sequences):
    """Get consensus sequence from multiple sequences"""
    if not sequences:
        return ""
    
    consensus = []
    for i in range(len(sequences[0])):
        bases = [seq[i] for seq in sequences if i < len(seq)]
        if bases:
            # Most frequent base
            most_common = max(set(bases), key=bases.count)
            consensus.append(most_common)
    
    return ''.join(consensus)

def count_mismatches(seq1, seq2):
    """Count mismatches between two sequences"""
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def is_primer_region_suitable(sequences, max_mismatches=6):
    """
    Check if primer region is suitable for RPA
    Allows up to 6 mismatches total across all sequences
    """
    if not sequences or len(sequences) < 2:
        return False
    
    # Check for gaps or Ns - these make regions unsuitable
    for seq in sequences:
        if '-' in seq or 'N' in seq:
            return False
    
    # Get consensus sequence
    consensus = get_consensus(sequences)
    if not consensus:
        return False
    
    # Count total mismatches across all sequences
    total_mismatches = 0
    for seq in sequences:
        mismatches = count_mismatches(seq, consensus)
        total_mismatches += mismatches
        # Early exit if too many mismatches
        if total_mismatches > max_mismatches:
            return False
    
    # Basic quality checks for consensus
    gc_content = gc_fraction(consensus) * 100
    if gc_content < 20 or gc_content > 80:  # Relaxed for initial screening
        return False
    
    # Check for extreme homopolymers (6+ bases)
    if re.search(r'([ATGC])\1{5,}', consensus):
        return False
    
    return True

def extract_conserved_targets(xmfa_file, target_size=18, flank_min=30, flank_max=70):
    """
    Extract 100% conserved targets with suitable primer regions
    """
    candidates = []
    block_id = 0
    
    print(f"Extracting {target_size}nt conserved targets with flanks from {flank_min}–{flank_max} bp...")
    
    for block in parse_xmfa(xmfa_file):
        block_id += 1
        if block_id % 10 == 0:
            print(f"Processing block {block_id}...")
        
        sequences = list(block.values())
        genome_names = list(block.keys())
        
        if len(sequences) < 2:
            continue
        
        min_length = min(len(seq) for seq in sequences)
        sequences = [seq[:min_length] for seq in sequences]  # Align lengths

        min_required = flank_min + target_size + flank_min
        if min_length < min_required:
            continue

        
        # Scan for conserved targets
        for pos in range(flank_max, min_length - target_size - flank_max):
            target_seqs = [seq[pos:pos + target_size] for seq in sequences]
            if any('-' in s or 'N' in s for s in target_seqs):
                continue
            if len(set(target_seqs)) != 1:
                continue
            target_seq = target_seqs[0]

            found_valid_flanks = False
            for left_len in range(flank_min, flank_max + 1):
                for right_len in range(flank_min, flank_max + 1):
                    start = pos - left_len
                    end = pos + target_size + right_len
                    if start < 0 or end > min_length:
                        continue
                    left_flank = [seq[start:pos] for seq in sequences]
                    right_flank = [seq[pos + target_size:end] for seq in sequences]

                    if is_primer_region_suitable(left_flank) and is_primer_region_suitable(right_flank):
                        left_consensus = get_consensus(left_flank)
                        right_consensus = get_consensus(right_flank)
                        left_mismatches = sum(count_mismatches(seq, left_consensus) for seq in left_flank)
                        right_mismatches = sum(count_mismatches(seq, right_consensus) for seq in right_flank)

                        candidates.append({
                            'block_id': block_id,
                            'position': pos,
                            'target_sequence': target_seq,
                            'target_gc': round(gc_fraction(target_seq) * 100, 1),
                            'left_region_consensus': left_consensus,
                            'right_region_consensus': right_consensus,
                            'left_region_gc': round(gc_fraction(left_consensus) * 100, 1),
                            'right_region_gc': round(gc_fraction(right_consensus) * 100, 1),
                            'left_region_mismatches': left_mismatches,
                            'right_region_mismatches': right_mismatches,
                            'total_mismatches': left_mismatches + right_mismatches,
                            'flank_left': left_len,
                            'flank_right': right_len,
                            'genome_count': len(sequences),
                            'genome_names': '|'.join(genome_names[:3]) + ('|...' if len(genome_names) > 3 else '')
                        })
                        found_valid_flanks = True
                        break
                if found_valid_flanks:
                    break
        
        if block_id % 50 == 0:
            print(f"  Processed {block_id} blocks, found {len(candidates)} candidates")
    
    print(f"Found {len(candidates)} candidate targets for RPA validation")
    return candidates

def main():
    if len(sys.argv) != 3:
        print("Usage: python extract_conserved_targets.py input.xmfa output.csv")
        print("Extracts 100% conserved 18nt targets with 60nt primer regions (≤6 mismatches)")
        sys.exit(1)
    
    xmfa_file = sys.argv[1]
    output_file = sys.argv[2]
    
    print("Starting conserved target extraction...")
    
    # Extract candidates
    candidates = extract_conserved_targets(xmfa_file)
    
    if not candidates:
        print("No suitable targets found!")
        # Create empty file for Snakemake
        pd.DataFrame().to_csv(output_file, index=False)
        sys.exit(0)
    
    # Sort by quality metrics
    df = pd.DataFrame(candidates)
    df = df.sort_values(['total_mismatches', 'target_gc', 'genome_count'], 
                       ascending=[True, False, False])
    
    # Export results
    df.to_csv(output_file, index=False)
    
    print(f"\nExtracted {len(candidates)} targets to {output_file}")
    
    # Summary statistics
    unique_targets = df['target_sequence'].nunique()
    avg_mismatches = df['total_mismatches'].mean()
    
    print(f"Summary:")
    print(f"  Unique target sequences: {unique_targets}")
    print(f"  Average mismatches in primer regions: {avg_mismatches:.1f}")
    print(f"  Target GC range: {df['target_gc'].min():.1f}% - {df['target_gc'].max():.1f}%")
    print(f"  Genome count range: {df['genome_count'].min()} - {df['genome_count'].max()}")
    
    # Show top candidates
    print(f"\nTop 5 candidates (lowest mismatch count):")
    print("=" * 100)
    for _, row in df.head().iterrows():
        print(f"Target: {row['target_sequence']} (GC: {row['target_gc']}%, Genomes: {row['genome_count']})")
        print(f"  Left region:  {row['left_region_consensus'][:30]}... (mismatches: {row['left_region_mismatches']})")
        print(f"  Right region: {row['right_region_consensus'][:30]}... (mismatches: {row['right_region_mismatches']})")
        print(f"  Total mismatches: {row['total_mismatches']}")
        print()

if __name__ == "__main__":
    main()