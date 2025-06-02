"""
Validate RPA primer candidates using PrimedRPA
Takes CSV from extract_conserved_targets.py and validates with PrimedRPA
"""

import sys
import os
import subprocess
import tempfile
import shutil
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

def check_primedrpa():
    """Check if PrimedRPA is installed"""
    try:
        result = subprocess.run(['PrimedRPA', '--help'], 
                              capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except:
        return False

def create_fasta_for_validation(target_seq, left_region, right_region, candidate_id):
    """
    Create FASTA sequence for PrimedRPA validation
    Format: left_region + target + right_region (gap-free)
    """
    # Remove gaps and create full sequence
    full_sequence = (left_region + target_seq + right_region).replace('-', '')
    
    # Create sequence record
    record = SeqRecord(
        Seq(full_sequence),
        id=f"candidate_{candidate_id}",
        description=f"RPA_target_{candidate_id} target:{target_seq}"
    )
    
    return record

def run_primedrpa(fasta_file, output_dir, candidate_id):
    """
    Run PrimedRPA on a single candidate
    Returns list of primer pairs or empty list if failed
    """
    candidate_dir = os.path.join(output_dir, f"candidate_{candidate_id}")
    os.makedirs(candidate_dir, exist_ok=True)
    
    cmd = [
        'PrimedRPA',
        '--input', fasta_file,
        '--output', candidate_dir,
        '--minlength', '30',    # Minimum primer length
        '--maxlength', '35',    # Maximum primer length  
        '--maxsize', '120',     # Maximum amplicon size
        '--minsize', '80',      # Minimum amplicon size
        '--gc_min', '30',       # Minimum GC content
        '--gc_max', '70',       # Maximum GC content
        '--max_hairpin', '-8',  # Maximum hairpin stability
        '--max_homodimer', '-8', # Maximum homodimer stability
        '--max_heterodimer', '-8' # Maximum heterodimer stability
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        
        if result.returncode == 0:
            return parse_primedrpa_results(candidate_dir)
        else:
            print(f"    PrimedRPA failed for candidate {candidate_id}: {result.stderr}")
            return []
            
    except subprocess.TimeoutExpired:
        print(f"    PrimedRPA timeout for candidate {candidate_id}")
        return []
    except Exception as e:
        print(f"    PrimedRPA error for candidate {candidate_id}: {e}")
        return []

def parse_primedrpa_results(output_dir):
    """
    Parse PrimedRPA output files and extract primer information
    """
    primer_pairs = []
    
    # Look for CSV output files
    for filename in os.listdir(output_dir):
        if filename.endswith('.csv'):
            filepath = os.path.join(output_dir, filename)
            
            try:
                df = pd.read_csv(filepath)
                
                # Find relevant columns (PrimedRPA output format may vary)
                forward_col = None
                reverse_col = None
                length_col = None
                
                for col in df.columns:
                    col_lower = col.lower()
                    if 'forward' in col_lower or 'primer' in col_lower:
                        if forward_col is None:
                            forward_col = col
                    elif 'reverse' in col_lower:
                        reverse_col = col
                    elif 'length' in col_lower or 'size' in col_lower or 'amplicon' in col_lower:
                        length_col = col
                
                if forward_col and reverse_col:
                    for _, row in df.iterrows():
                        primer_pair = {
                            'forward_primer': str(row[forward_col]),
                            'reverse_primer': str(row[reverse_col]),
                            'amplicon_length': row[length_col] if length_col else None,
                            'validation_source': filename
                        }
                        primer_pairs.append(primer_pair)
                        
            except Exception as e:
                print(f"    Error parsing {filepath}: {e}")
                continue
    
    return primer_pairs

def validate_targets(input_csv, temp_dir):
    """
    Validate all targets from CSV using PrimedRPA
    """
    if not check_primedrpa():
        print("ERROR: PrimedRPA not found!")
        print("Install with: conda install -c bioconda primedrpa")
        print("Or: pip install PrimedRPA")
        return []
    
    # Read candidates
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        print(f"Error reading input CSV: {e}")
        return []
    
    if df.empty:
        print("No candidates to validate (empty input file)")
        return []
    
    print(f"Validating {len(df)} candidates with PrimedRPA...")
    
    validated_results = []
    os.makedirs(temp_dir, exist_ok=True)
    
    for idx, row in df.iterrows():
        candidate_id = idx + 1
        
        if candidate_id % 10 == 0:
            print(f"  Processing candidate {candidate_id}/{len(df)}")
        
        # Create FASTA file for this candidate
        fasta_record = create_fasta_for_validation(
            row['target_sequence'],
            row['left_region_consensus'], 
            row['right_region_consensus'],
            candidate_id
        )
        
        fasta_file = os.path.join(temp_dir, f"candidate_{candidate_id}.fasta")
        SeqIO.write([fasta_record], fasta_file, "fasta")
        
        # Run PrimedRPA
        primer_pairs = run_primedrpa(fasta_file, temp_dir, candidate_id)
        
        # Process results
        for primer_pair in primer_pairs:
            result = {
                'candidate_id': candidate_id,
                'block_id': row['block_id'],
                'position': row['position'],
                'target_sequence': row['target_sequence'],
                'target_gc': row['target_gc'],
                'forward_primer': primer_pair['forward_primer'],
                'reverse_primer': primer_pair['reverse_primer'],
                'forward_primer_gc': round(gc_fraction(primer_pair['forward_primer']), 1),
                'reverse_primer_gc': round(gc_fraction(primer_pair['reverse_primer']), 1),
                'amplicon_length': primer_pair['amplicon_length'],
                'genome_count': row['genome_count'],
                'original_total_mismatches': row['total_mismatches'],
                'genome_names': row['genome_names'],
                'validation_source': primer_pair['validation_source']
            }
            validated_results.append(result)
    
    print(f"PrimedRPA validation complete: {len(validated_results)} valid primer pairs found")
    return validated_results

def export_results(results, output_file):
    """
    Export validated results with summary statistics
    """
    if not results:
        print("No validated primers to export!")
        # Create empty CSV for Snakemake
        pd.DataFrame().to_csv(output_file, index=False)
        return
    
    # Create DataFrame and sort by quality
    df = pd.DataFrame(results)
    df = df.sort_values(['target_sequence', 'amplicon_length', 'original_total_mismatches'])
    
    # Export main results
    df.to_csv(output_file, index=False)
    
    # Create additional output files
    base_name = output_file.replace('.csv', '')
    
    # Export unique targets FASTA
    targets_fasta = f"{base_name}_targets.fasta"
    export_targets_fasta(results, targets_fasta)
    
    # Export primers FASTA  
    primers_fasta = f"{base_name}_primers.fasta"
    export_primers_fasta(results, primers_fasta)
    
    print(f"\nValidated RPA results exported:")
    print(f"  Main results: {output_file}")
    print(f"  Target sequences: {targets_fasta}")
    print(f"  Primer sequences: {primers_fasta}")
    
    # Summary statistics
    unique_targets = df['target_sequence'].nunique()
    avg_amplicon = df['amplicon_length'].mean() if 'amplicon_length' in df.columns else 0
    
    print(f"\nSummary:")
    print(f"  Total validated primer pairs: {len(results)}")
    print(f"  Unique target sequences: {unique_targets}")
    print(f"  Average amplicon length: {avg_amplicon:.1f} bp")
    
    # Show top results
    print(f"\nTop 5 validated RPA targets:")
    print("=" * 100)
    for _, row in df.head().iterrows():
        print(f"Target: {row['target_sequence']} (GC: {row['target_gc']}%)")
        print(f"  Forward:  {row['forward_primer']} (GC: {row['forward_primer_gc']}%)")
        print(f"  Reverse:  {row['reverse_primer']} (GC: {row['reverse_primer_gc']}%)")
        print(f"  Amplicon: {row['amplicon_length']} bp, Genomes: {row['genome_count']}")
        print()

def export_targets_fasta(results, fasta_file):
    """Export unique target sequences"""
    seen_targets = set()
    records = []
    
    for i, result in enumerate(results):
        target_seq = result['target_sequence']
        if target_seq not in seen_targets:
            seen_targets.add(target_seq)
            record = SeqRecord(
                Seq(target_seq),
                id=f"RPA_target_{len(records)+1}",
                description=f"GC:{result['target_gc']}% genomes:{result['genome_count']}"
            )
            records.append(record)
    
    SeqIO.write(records, fasta_file, "fasta")
    print(f"  Exported {len(records)} unique targets")

def export_primers_fasta(results, fasta_file):
    """Export all primer sequences"""
    records = []
    
    for i, result in enumerate(results):
        # Forward primer
        forward_record = SeqRecord(
            Seq(result['forward_primer']),
            id=f"RPA_primer_{i+1}_F",
            description=f"Forward for {result['target_sequence']}"
        )
        records.append(forward_record)
        
        # Reverse primer
        reverse_record = SeqRecord(
            Seq(result['reverse_primer']),
            id=f"RPA_primer_{i+1}_R",
            description=f"Reverse for {result['target_sequence']}"
        )
        records.append(reverse_record)
    
    SeqIO.write(records, fasta_file, "fasta")
    print(f"  Exported {len(records)} primer sequences")

def main():
    if len(sys.argv) != 4:
        print("Usage: python validate_rpa_primers.py input.csv output.csv temp_dir")
        print("Validates RPA candidates using PrimedRPA")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_csv = sys.argv[2] 
    temp_dir = sys.argv[3]
    
    print("Starting RPA primer validation with PrimedRPA...")
    
    # Validate targets
    validated_results = validate_targets(input_csv, temp_dir)
    
    # Export results
    export_results(validated_results, output_csv)
    
    # Cleanup temp directory if requested
    if os.path.exists(temp_dir) and len(validated_results) > 0:
        try:
            shutil.rmtree(temp_dir)
            print(f"Cleaned up temporary directory: {temp_dir}")
        except:
            print(f"Could not clean up {temp_dir} - manual cleanup may be needed")
    
    print("RPA primer validation complete!")

if __name__ == "__main__":
    main()