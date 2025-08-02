#!/usr/bin/env python3
"""
RPA Primer Designer for Variable Target Sequences
Designed for Rickettsia species-specific amplification

Requirements:
- Primer length: 30-35 bp
- Primer GC: 30-70%
- Amplicon size: 200-300 bp
- Amplicon GC: 35-60%
- Avoid repeats, secondary structures, 3' complementarity
- Prefer G/C at 3' end
"""

import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re
from collections import defaultdict
from itertools import combinations
import numpy as np

def extract_genome_name(seq_id):
    """Extract genome name from sequence ID"""
    return seq_id.split(';')[0] if ';' in seq_id else seq_id

def parse_alignment_blocks_dynamic(fasta_file):
    """
    Parse alignment file into blocks dynamically
    Returns generator of (block_id, block_sequences)
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        print("ERROR: Empty alignment file!")
        return
    
    # Find first genome to use as block marker
    first_genome = extract_genome_name(records[0].id)
    
    blocks = []
    current_block = []
    
    for rec in records:
        genome_name = extract_genome_name(rec.id)
        
        if genome_name == first_genome and current_block:
            blocks.append(current_block)
            current_block = [rec]
        else:
            current_block.append(rec)
    
    if current_block:
        blocks.append(current_block)
    
    print(f"Parsed {len(blocks)} alignment blocks")
    
    for i, block in enumerate(blocks, 1):
        yield i, block

def find_target_in_block(block, target_sequence, tolerance=2):
    """
    Find target sequence position in alignment block
    Returns (start_pos, end_pos) or None if not found
    """
    if not block:
        return None
    
    # Get consensus sequence from first record (remove gaps)
    consensus = str(block[0].seq).replace('-', '').upper()
    target = target_sequence.upper()
    
    # Try exact match first
    pos = consensus.find(target)
    if pos != -1:
        return pos, pos + len(target)
    
    # Try with some mismatches allowed
    target_len = len(target)
    for start in range(len(consensus) - target_len + 1):
        segment = consensus[start:start + target_len]
        mismatches = sum(1 for a, b in zip(target, segment) if a != b)
        if mismatches <= tolerance:
            return start, start + target_len
    
    return None

def get_block_consensus_regions(block, target_start, target_end, amplicon_size):
    """
    Extract consensus regions around target for primer design
    Returns (left_region, target_region, right_region)
    """
    if not block:
        return None, None, None
    
    # Get gapped alignment sequence
    alignment_seq = str(block[0].seq).upper()
    
    # Convert ungapped positions to gapped positions
    ungapped_to_gapped = {}
    ungapped_pos = 0
    for gapped_pos, char in enumerate(alignment_seq):
        if char != '-':
            ungapped_to_gapped[ungapped_pos] = gapped_pos
            ungapped_pos += 1
    
    # Find gapped coordinates
    try:
        gapped_target_start = ungapped_to_gapped[target_start]
        gapped_target_end = ungapped_to_gapped[target_end - 1] + 1
    except KeyError:
        return None, None, None
    
    # Calculate region boundaries for amplicon size
    max_primer_len = 35
    needed_flanking = (amplicon_size - (target_end - target_start)) // 2
    
    # Extract regions with some buffer
    buffer = max(needed_flanking + 50, 150)  # Extra buffer for primer search
    
    left_start = max(0, gapped_target_start - buffer)
    right_end = min(len(alignment_seq), gapped_target_end + buffer)
    
    left_region = alignment_seq[left_start:gapped_target_start]
    target_region = alignment_seq[gapped_target_start:gapped_target_end]
    right_region = alignment_seq[gapped_target_end:right_end]
    
    return left_region, target_region, right_region

class RPAPrimerDesigner:
    """RPA Primer Designer with specific requirements"""
    
    def __init__(self):
        self.primer_requirements = {
            'min_length': 30,
            'max_length': 35,
            'min_gc': 30,
            'max_gc': 70,
            'preferred_3prime': ['G', 'C'],
            'avoid_3prime': ['T']
        }
        
        self.amplicon_requirements = {
            'min_size': 200,
            'max_size': 300,
            'min_gc': 35,
            'max_gc': 60
        }
        
        self.validation_thresholds = {
            'max_homopolymer': 4,
            'max_dinuc_repeat': 6,
            'max_hairpin_score': 8,
            'max_self_comp': 4,
            'max_3prime_comp': 2,
            'max_pair_3prime_comp': 3
        }
    
    def validate_primer_sequence(self, primer_seq):
        """Validate individual primer against requirements"""
        issues = []
        penalty = 0
        
        # Length check
        length = len(primer_seq)
        if not (self.primer_requirements['min_length'] <= length <= self.primer_requirements['max_length']):
            issues.append(f"Invalid length ({length}bp, need {self.primer_requirements['min_length']}-{self.primer_requirements['max_length']})")
            penalty += 50
        
        # GC content
        gc_content = gc_fraction(primer_seq) * 100
        if not (self.primer_requirements['min_gc'] <= gc_content <= self.primer_requirements['max_gc']):
            issues.append(f"Invalid GC content ({gc_content:.1f}%, need {self.primer_requirements['min_gc']}-{self.primer_requirements['max_gc']}%)")
            penalty += 30
        
        # 3' end preference
        three_prime = primer_seq[-1]
        if three_prime in self.primer_requirements['preferred_3prime']:
            penalty -= 5  # Bonus
        elif three_prime in self.primer_requirements['avoid_3prime']:
            issues.append(f"Poor 3' end ({three_prime})")
            penalty += 15
        
        # Homopolymer check
        for base in 'ATGC':
            max_homo = self.find_max_homopolymer(primer_seq, base)
            if max_homo > self.validation_thresholds['max_homopolymer']:
                issues.append(f"Homopolymer: {base}×{max_homo}")
                penalty += 20
        
        # Dinucleotide repeats
        max_dinuc = self.find_max_dinucleotide_repeat(primer_seq)
        if max_dinuc > self.validation_thresholds['max_dinuc_repeat']:
            issues.append(f"Dinucleotide repeat: {max_dinuc}bp")
            penalty += 15
        
        # Secondary structures
        hairpin_score = self.calculate_hairpin_potential(primer_seq)
        if hairpin_score > self.validation_thresholds['max_hairpin_score']:
            issues.append(f"Hairpin potential: {hairpin_score}")
            penalty += 25
        
        # Self-complementarity
        self_comp = self.calculate_self_complementarity(primer_seq)
        if self_comp > self.validation_thresholds['max_self_comp']:
            issues.append(f"Self-complementarity: {self_comp}bp")
            penalty += 20
        
        # 3' self-complementarity
        three_prime_comp = self.calculate_3prime_self_complementarity(primer_seq)
        if three_prime_comp > self.validation_thresholds['max_3prime_comp']:
            issues.append(f"3' self-complementarity: {three_prime_comp}bp")
            penalty += 30
        
        return {
            'is_valid': len([i for i in issues if 'Invalid' in i]) == 0,
            'issues': issues,
            'penalty': penalty,
            'gc_content': gc_content,
            'three_prime': three_prime,
            'hairpin_score': hairpin_score,
            'self_comp': self_comp,
            'three_prime_comp': three_prime_comp
        }
    
    def validate_amplicon(self, amplicon_seq):
        """Validate amplicon against requirements"""
        issues = []
        penalty = 0
        
        # Remove gaps
        clean_amplicon = amplicon_seq.replace('-', '')
        
        # Size check
        size = len(clean_amplicon)
        if not (self.amplicon_requirements['min_size'] <= size <= self.amplicon_requirements['max_size']):
            issues.append(f"Invalid amplicon size ({size}bp, need {self.amplicon_requirements['min_size']}-{self.amplicon_requirements['max_size']})")
            penalty += 50
        
        # GC content
        gc_content = gc_fraction(clean_amplicon) * 100
        if not (self.amplicon_requirements['min_gc'] <= gc_content <= self.amplicon_requirements['max_gc']):
            issues.append(f"Invalid amplicon GC ({gc_content:.1f}%, need {self.amplicon_requirements['min_gc']}-{self.amplicon_requirements['max_gc']}%)")
            penalty += 30
        
        # Check for direct repeats
        repeat_score = self.find_direct_repeats(clean_amplicon)
        if repeat_score > 15:
            issues.append(f"Direct repeats detected: {repeat_score}")
            penalty += 20
        
        # Check for palindromes
        palindrome_score = self.find_palindromes(clean_amplicon)
        if palindrome_score > 12:
            issues.append(f"Palindromes detected: {palindrome_score}")
            penalty += 15
        
        return {
            'is_valid': len([i for i in issues if 'Invalid' in i]) == 0,
            'issues': issues,
            'penalty': penalty,
            'size': size,
            'gc_content': gc_content,
            'repeat_score': repeat_score,
            'palindrome_score': palindrome_score
        }
    
    def validate_primer_pair(self, forward_primer, reverse_primer):
        """Validate primer pair interactions"""
        issues = []
        penalty = 0
        
        # 3' end complementarity between primers
        pair_3prime_comp = self.calculate_primer_pair_3prime_complementarity(forward_primer, reverse_primer)
        if pair_3prime_comp > self.validation_thresholds['max_pair_3prime_comp']:
            issues.append(f"Primer pair 3' complementarity: {pair_3prime_comp}bp")
            penalty += 40
        
        # Cross-complementarity
        cross_comp = self.calculate_cross_complementarity(forward_primer, reverse_primer)
        if cross_comp > 6:
            issues.append(f"Cross-complementarity: {cross_comp}bp")
            penalty += 25
        
        return {
            'issues': issues,
            'penalty': penalty,
            'pair_3prime_comp': pair_3prime_comp,
            'cross_comp': cross_comp
        }
    
    def design_primers_for_target(self, left_region, target_region, right_region, target_info):
        """Design primer pairs for a specific target"""
        
        # Remove gaps for primer design
        left_clean = left_region.replace('-', '')
        target_clean = target_region.replace('-', '')
        right_clean = right_region.replace('-', '')
        
        primer_pairs = []
        
        # Try different primer lengths
        for f_len in range(self.primer_requirements['min_length'], self.primer_requirements['max_length'] + 1):
            for r_len in range(self.primer_requirements['min_length'], self.primer_requirements['max_length'] + 1):
                
                # Calculate positions for target amplicon sizes
                for target_amp_size in range(self.amplicon_requirements['min_size'], 
                                           self.amplicon_requirements['max_size'] + 1, 10):
                    
                    # Calculate how much flanking sequence we need
                    needed_flanking = target_amp_size - len(target_clean) - f_len - r_len
                    if needed_flanking < 0:
                        continue
                    
                    # Try different distributions of flanking sequence
                    for left_flank in range(max(0, needed_flanking - 50), min(needed_flanking + 1, 51)):
                        right_flank = needed_flanking - left_flank
                        
                        # Check if we have enough sequence
                        if left_flank > len(left_clean) - f_len or right_flank > len(right_clean) - r_len:
                            continue
                        
                        # Extract primer sequences
                        f_start = len(left_clean) - f_len - left_flank
                        f_end = len(left_clean) - left_flank
                        forward_primer = left_clean[f_start:f_end]
                        
                        r_start = right_flank
                        r_end = right_flank + r_len
                        reverse_region = right_clean[r_start:r_end]
                        reverse_primer = str(Seq(reverse_region).reverse_complement())
                        
                        # Skip if primers contain invalid characters
                        if any(char not in 'ATGC' for char in forward_primer + reverse_primer):
                            continue
                        
                        # Extract amplicon components
                        left_flanking = left_clean[f_start:len(left_clean) - left_flank]
                        right_flanking = right_clean[:r_end]
                        
                        # Construct full amplicon
                        amplicon = left_flanking + target_clean + right_flanking
                        
                        # Calculate positions within amplicon
                        forward_start_in_amplicon = 0
                        forward_end_in_amplicon = len(forward_primer)
                        target_start_in_amplicon = len(left_flanking)
                        target_end_in_amplicon = len(left_flanking) + len(target_clean)
                        reverse_start_in_amplicon = len(amplicon) - len(reverse_primer)
                        reverse_end_in_amplicon = len(amplicon)
                        
                        # Create reverse complement of the actual reverse region for clarity
                        actual_reverse_region = right_clean[r_start:r_end]
                        
                        # Validate primers
                        f_validation = self.validate_primer_sequence(forward_primer)
                        r_validation = self.validate_primer_sequence(reverse_primer)
                        amp_validation = self.validate_amplicon(amplicon)
                        pair_validation = self.validate_primer_pair(forward_primer, reverse_primer)
                        
                        # Calculate total penalty
                        total_penalty = (f_validation['penalty'] + r_validation['penalty'] + 
                                       amp_validation['penalty'] + pair_validation['penalty'])
                        
                        # Check if design is valid
                        is_valid = (f_validation['is_valid'] and r_validation['is_valid'] and 
                                  amp_validation['is_valid'] and len(pair_validation['issues']) == 0)
                        
                        primer_data = {
                            # === MAIN SEQUENCES ===
                            'forward_primer': forward_primer,
                            'reverse_primer': reverse_primer,
                            'target_sequence_in_amplicon': target_clean,
                            'full_amplicon_sequence': amplicon,
                            
                            # === AMPLICON BREAKDOWN ===
                            'left_flanking_sequence': left_flanking,
                            'right_flanking_sequence': right_flanking,
                            'reverse_region_on_template': actual_reverse_region,  # Actual sequence on template (not rev-comp)
                            
                            # === POSITIONS IN AMPLICON ===
                            'forward_start_pos': forward_start_in_amplicon,
                            'forward_end_pos': forward_end_in_amplicon,
                            'target_start_pos': target_start_in_amplicon,
                            'target_end_pos': target_end_in_amplicon,
                            'reverse_start_pos': reverse_start_in_amplicon,
                            'reverse_end_pos': reverse_end_in_amplicon,
                            
                            # === SEQUENCE ANNOTATIONS ===
                            'amplicon_annotation': f"Forward({forward_start_in_amplicon}-{forward_end_in_amplicon})|Target({target_start_in_amplicon}-{target_end_in_amplicon})|Reverse({reverse_start_in_amplicon}-{reverse_end_in_amplicon})",
                            
                            # === LENGTH AND COMPOSITION ===
                            'forward_length': f_len,
                            'reverse_length': r_len,
                            'target_length_in_amplicon': len(target_clean),
                            'left_flanking_length': len(left_flanking),
                            'right_flanking_length': len(right_flanking),
                            'amplicon_size': len(amplicon),
                            
                            # === GC CONTENT ===
                            'forward_gc': f_validation['gc_content'],
                            'reverse_gc': r_validation['gc_content'],
                            'target_gc': round(gc_fraction(target_clean) * 100, 1),
                            'amplicon_gc': amp_validation['gc_content'],
                            
                            # === PRIMER CHARACTERISTICS ===
                            'forward_3prime': f_validation['three_prime'],
                            'reverse_3prime': r_validation['three_prime'],
                            
                            # === VALIDATION RESULTS ===
                            'total_penalty': total_penalty,
                            'is_valid': is_valid,
                            'forward_issues': f_validation['issues'],
                            'reverse_issues': r_validation['issues'],
                            'amplicon_issues': amp_validation['issues'],
                            'pair_issues': pair_validation['issues'],
                            'all_issues': (f_validation['issues'] + r_validation['issues'] + 
                                         amp_validation['issues'] + pair_validation['issues']),
                            **target_info
                        }
                        
                        primer_pairs.append(primer_data)
        
        # Sort by validity and penalty
        primer_pairs.sort(key=lambda x: (not x['is_valid'], x['total_penalty']))
        
        return primer_pairs[:10]  # Return top 10 designs
    
    def find_max_homopolymer(self, sequence, base):
        """Find maximum homopolymer length for a specific base"""
        max_len = current_len = 0
        for nucleotide in sequence:
            if nucleotide == base:
                current_len += 1
                max_len = max(max_len, current_len)
            else:
                current_len = 0
        return max_len
    
    def find_max_dinucleotide_repeat(self, sequence):
        """Find maximum dinucleotide repeat length"""
        max_repeat = 0
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            repeat_len = 2
            j = i + 2
            while j + 1 < len(sequence) and sequence[j:j+2] == dinuc:
                repeat_len += 2
                j += 2
            max_repeat = max(max_repeat, repeat_len)
        return max_repeat
    
    def calculate_hairpin_potential(self, sequence):
        """Calculate hairpin formation potential"""
        max_score = 0
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        
        for stem_len in range(4, min(8, len(sequence) // 3)):
            for loop_len in range(3, 8):
                for start in range(len(sequence) - 2 * stem_len - loop_len):
                    if start + 2 * stem_len + loop_len > len(sequence):
                        continue
                    
                    stem1 = sequence[start:start + stem_len]
                    stem2 = sequence[start + stem_len + loop_len:start + 2 * stem_len + loop_len]
                    
                    matches = sum(1 for a, b in zip(stem1, stem2[::-1]) 
                                if complement.get(a) == b)
                    
                    if matches >= stem_len - 1:
                        score = matches + (stem_len - loop_len) * 0.5
                        max_score = max(max_score, score)
        
        return max_score
    
    def calculate_self_complementarity(self, sequence):
        """Calculate self-complementarity excluding hairpins"""
        max_comp = 0
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        
        for i in range(1, len(sequence)):
            current_comp = 0
            for j in range(min(i, len(sequence) - i)):
                if j >= len(sequence) or (len(sequence) - i + j) >= len(sequence):
                    break
                if sequence[j] == complement.get(sequence[len(sequence) - i + j], 'N'):
                    current_comp += 1
                else:
                    break
            max_comp = max(max_comp, current_comp)
        
        return max_comp
    
    def calculate_3prime_self_complementarity(self, sequence):
        """Calculate 3' end self-complementarity"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        three_prime = sequence[-5:]
        max_comp = 0
        
        for i in range(len(sequence) - 5):
            comp_count = 0
            for j, base in enumerate(three_prime):
                if i + j < len(sequence) and sequence[i + j] == complement.get(base):
                    comp_count += 1
                else:
                    break
            max_comp = max(max_comp, comp_count)
        
        return max_comp
    
    def calculate_primer_pair_3prime_complementarity(self, primer1, primer2):
        """Calculate 3' complementarity between primer pair"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        
        p1_3prime = primer1[-5:]
        p2_3prime = primer2[-5:]
        
        max_comp = 0
        
        # Check forward primer 3' against reverse primer
        for i in range(len(primer2) - 4):
            comp_count = 0
            for j, base in enumerate(p1_3prime):
                if i + j < len(primer2) and primer2[i + j] == complement.get(base):
                    comp_count += 1
                else:
                    break
            max_comp = max(max_comp, comp_count)
        
        # Check reverse primer 3' against forward primer
        for i in range(len(primer1) - 4):
            comp_count = 0
            for j, base in enumerate(p2_3prime):
                if i + j < len(primer1) and primer1[i + j] == complement.get(base):
                    comp_count += 1
                else:
                    break
            max_comp = max(max_comp, comp_count)
        
        return max_comp
    
    def calculate_cross_complementarity(self, primer1, primer2):
        """Calculate cross-complementarity between primers"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        max_comp = 0
        
        # Check primer1 against primer2
        for i in range(len(primer2)):
            current_comp = 0
            for j in range(min(len(primer1), len(primer2) - i)):
                if primer1[j] == complement.get(primer2[i + j], 'N'):
                    current_comp += 1
                else:
                    if current_comp > 0:
                        break
                    current_comp = 0
            max_comp = max(max_comp, current_comp)
        
        return max_comp
    
    def find_direct_repeats(self, sequence):
        """Find direct repeats in sequence"""
        max_repeat_score = 0
        
        for length in range(6, min(len(sequence) // 2, 20)):
            for i in range(len(sequence) - length):
                pattern = sequence[i:i + length]
                for j in range(i + length, len(sequence) - length + 1):
                    if sequence[j:j + length] == pattern:
                        max_repeat_score = max(max_repeat_score, length)
        
        return max_repeat_score
    
    def find_palindromes(self, sequence):
        """Find palindromic sequences"""
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        max_palindrome = 0
        
        for length in range(6, min(len(sequence) + 1, 16)):
            for start in range(len(sequence) - length + 1):
                subseq = sequence[start:start + length]
                rev_comp = ''.join(complement.get(base, 'N') for base in subseq[::-1])
                if subseq == rev_comp:
                    max_palindrome = max(max_palindrome, length)
        
        return max_palindrome

def main():
    parser = argparse.ArgumentParser(description="Design RPA primers for variable target sequences")
    parser.add_argument("targets_csv", help="CSV file with variable targets from find_variable_targets.py")
    parser.add_argument("alignment_fasta", help="Partitioned alignment FASTA file")
    parser.add_argument("output_csv", help="Output CSV file for primer designs")
    args = parser.parse_args()
    
    # Load targets
    try:
        targets_df = pd.read_csv(args.targets_csv)
        print(f"Loaded {len(targets_df)} target candidates")
    except Exception as e:
        print(f"Error loading targets CSV: {e}")
        return
    
    # Initialize primer designer
    designer = RPAPrimerDesigner()
    
    # Process each unique target position
    unique_targets = targets_df.groupby(['block_id', 'start_pos', 'end_pos']).first().reset_index()
    print(f"Processing {len(unique_targets)} unique target positions")
    
    all_primer_designs = []
    
    # Parse alignment blocks
    alignment_blocks = {}
    for block_id, block in parse_alignment_blocks_dynamic(args.alignment_fasta):
        alignment_blocks[block_id] = block
    
    for idx, target_row in unique_targets.iterrows():
        if (idx + 1) % 10 == 0:
            print(f"Processing target {idx + 1}/{len(unique_targets)}")
        
        block_id = target_row['block_id']
        target_seq = target_row['sequence_variant']
        
        # Get alignment block
        if block_id not in alignment_blocks:
            print(f"Warning: Block {block_id} not found in alignment")
            continue
        
        block = alignment_blocks[block_id]
        
        # Find target position in block
        target_pos = find_target_in_block(block, target_seq)
        if target_pos is None:
            print(f"Warning: Target sequence not found in block {block_id}")
            continue
        
        target_start, target_end = target_pos
        
        # Get consensus regions
        regions = get_block_consensus_regions(block, target_start, target_end, 250)  # Middle of 200-300 range
        if regions[0] is None:
            print(f"Warning: Could not extract regions for block {block_id}")
            continue
        
        left_region, target_region, right_region = regions
        
        # Prepare target info
        target_info = {
            'candidate_id': idx + 1,
            'block_id': block_id,
            'target_sequence': target_seq,
            'target_start': target_start,
            'target_end': target_end,
            'target_length': len(target_seq),
            'uniqueness_score': target_row.get('uniqueness_score', 'unknown'),
            'num_variants': target_row.get('num_variants', 'unknown'),
            'genome_count': target_row.get('genome_count', 'unknown')
        }
        
        # Design primers
        try:
            primer_designs = designer.design_primers_for_target(
                left_region, target_region, right_region, target_info
            )
            
            # Add rank to each design
            for rank, design in enumerate(primer_designs, 1):
                design['primer_rank'] = rank
                all_primer_designs.append(design)
                
        except Exception as e:
            print(f"Error designing primers for target {idx + 1}: {e}")
            continue
    
    # Save results
    if all_primer_designs:
        results_df = pd.DataFrame(all_primer_designs)
        results_df.to_csv(args.output_csv, index=False)
        
        # Print summary
        valid_designs = results_df[results_df['is_valid'] == True]
        print(f"\n=== RPA PRIMER DESIGN SUMMARY ===")
        print(f"Total designs generated: {len(results_df)}")
        print(f"Valid designs: {len(valid_designs)} ({len(valid_designs)/len(results_df)*100:.1f}%)")
        
        if len(valid_designs) > 0:
            print(f"\nValid designs statistics:")
            print(f"  Average amplicon size: {valid_designs['amplicon_size'].mean():.1f} bp")
            print(f"  Average forward primer length: {valid_designs['forward_length'].mean():.1f} bp")
            print(f"  Average reverse primer length: {valid_designs['reverse_length'].mean():.1f} bp")
            print(f"  Average forward GC: {valid_designs['forward_gc'].mean():.1f}%")
            print(f"  Average reverse GC: {valid_designs['reverse_gc'].mean():.1f}%")
            print(f"  Average amplicon GC: {valid_designs['amplicon_gc'].mean():.1f}%")
            
            # Show top designs
            top_designs = valid_designs.head(5)
            print(f"\nTop 5 valid designs:")
            for i, (_, row) in enumerate(top_designs.iterrows(), 1):
                print(f"{i}. Target {row['candidate_id']}, Rank {row['primer_rank']}")
                print(f"   Amplicon: {row['amplicon_size']}bp, GC: {row['amplicon_gc']:.1f}%")
                print(f"   Forward:  {row['forward_primer']} ({row['forward_gc']:.1f}% GC, 3'={row['forward_3prime']})")
                print(f"   Target:   {row['target_sequence_in_amplicon']} ({row['target_gc']:.1f}% GC)")
                print(f"   Reverse:  {row['reverse_primer']} ({row['reverse_gc']:.1f}% GC, 3'={row['reverse_3prime']})")
                print(f"   Annotation: {row['amplicon_annotation']}")
                print()
            
            print(f"\nAmplicon structure example (first valid design):")
            if len(valid_designs) > 0:
                first_design = valid_designs.iloc[0]
                print(f"Full amplicon: {first_design['full_amplicon_sequence']}")
                print(f"Forward:       {first_design['forward_primer']}")
                print(f"Left flank:    {first_design['left_flanking_sequence']}")
                print(f"Target:        {first_design['target_sequence_in_amplicon']}")
                print(f"Right flank:   {first_design['right_flanking_sequence']}")
                print(f"Reverse (RC):  {first_design['reverse_primer']}")
                print(f"Reverse (Tmpl):{first_design['reverse_region_on_template']}")
                print()
        
        print(f"Results saved to: {args.output_csv}")
    else:
        print("No primer designs generated!")
        # Save empty file
        pd.DataFrame().to_csv(args.output_csv, index=False)

if __name__ == "__main__":
    main()