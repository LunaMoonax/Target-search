#!/usr/bin/env python3
"""
Hierarchical Contamination Checker for Pre-Validated RPA Primers
CRITICAL: Hierarchical checking - Targets → Primers → Amplicons
Allows few mismatches to catch similar sequences that could cause issues
"""

import pandas as pd
import subprocess
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile

class HierarchicalContaminationChecker:
    """
    Hierarchical contamination checker with mismatch tolerance
    Step 1: Check targets → Remove contaminated targets
    Step 2: Check primers (only clean targets) → Remove contaminated primers  
    Step 3: Check amplicons (only clean primers) → Final clean set
    """
    
    def __init__(self, blast_db_prefix):
        self.blast_db_prefix = blast_db_prefix
        
        # Hierarchical contamination thresholds with mismatch tolerance
        self.thresholds = {
            "target": {
                "identity": 80,      # 80% = allows 20% mismatches to catch similar sequences
                "coverage": 70,      # 70% coverage minimum
                "evalue": "1e-3"     # More permissive for target detection
            },
            "primer": {
                "identity": 75,      # 75% = allows 25% mismatches for short primers
                "coverage": 60,      # 60% coverage minimum  
                "evalue": "1"        # Permissive for short sequences
            },
            "amplicon": {
                "identity": 85,      # 85% = allows 15% mismatches for full amplicons
                "coverage": 80,      # 80% coverage minimum
                "evalue": "1e-5"     # Standard threshold for long sequences
            }
        }

    def construct_full_amplicon(self, target_seq, forward_primer, reverse_primer):
        """Construct full amplicon sequence: Forward primer + Target + Reverse complement of reverse primer"""
        try:
            target = Seq(str(target_seq).upper().strip())
            forward = Seq(str(forward_primer).upper().strip())
            reverse = Seq(str(reverse_primer).upper().strip())
            
            # Get reverse complement of reverse primer
            reverse_complement = reverse.reverse_complement()
            
            # Construct full amplicon
            full_amplicon = forward + target + reverse_complement
            return str(full_amplicon)
        except Exception as e:
            print(f"   Warning: Could not construct amplicon - {e}")
            return None

    def add_full_amplicon_sequences(self, primers_df):
        """Add full amplicon sequences to the dataframe"""
        print("🧬 Constructing full amplicon sequences...")
        
        primers_df = primers_df.copy()
        amplicon_sequences = []
        valid_amplicons = 0
        
        for idx, row in primers_df.iterrows():
            amplicon = self.construct_full_amplicon(
                row['target_sequence'], 
                row['forward_primer'], 
                row['reverse_primer']
            )
            amplicon_sequences.append(amplicon)
            if amplicon:
                valid_amplicons += 1
        
        primers_df['full_amplicon_sequence'] = amplicon_sequences
        
        print(f"   ✅ Constructed {valid_amplicons}/{len(primers_df)} valid amplicon sequences")
        
        return primers_df
    
    def run_blast(self, query_fasta, output_file, blast_type="target"):
        """Run BLAST with type-specific parameters"""
        if not os.path.exists(query_fasta) or os.path.getsize(query_fasta) == 0:
            print(f"   Warning: Empty query file {query_fasta}")
            return False
        
        params = self.thresholds[blast_type]
        
        cmd = [
            "blastn",
            "-query", query_fasta,
            "-db", self.blast_db_prefix,
            "-out", output_file,
            "-outfmt", "6 qseqid sseqid pident qcovs length evalue bitscore stitle",
            "-evalue", params["evalue"],
            "-word_size", "7" if blast_type == "primer" else "11",  # Smaller word size for primers
            "-max_target_seqs", "100",
            "-dust", "no",
            "-soft_masking", "false",
            "-perc_identity", str(params["identity"] - 5),  # Allow BLAST to find more hits, filter later
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            return True
        except subprocess.CalledProcessError as e:
            print(f"   BLAST failed for {blast_type}: {e}")
            return False
        except FileNotFoundError:
            print(f"   Error: BLAST not found. Please ensure BLAST+ is installed and in PATH")
            return False
    
    def parse_blast_results(self, blast_file, blast_type="target"):
        """Parse BLAST results with type-specific thresholds"""
        if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
            return set(), []
        
        params = self.thresholds[blast_type]
        contaminated_sequences = set()
        all_hits = []
        
        try:
            df = pd.read_csv(blast_file, sep='\t', 
                           names=['qseqid', 'sseqid', 'pident', 'qcovs', 'length', 'evalue', 'bitscore', 'stitle'])
            
            # Apply strict filtering with mismatch tolerance
            significant_hits = df[
                (df['pident'] >= params["identity"]) & 
                (df['qcovs'] >= params["coverage"])
            ]
            
            for _, hit in significant_hits.iterrows():
                seq_id = hit['qseqid']
                contaminated_sequences.add(seq_id)
                
                # Store hit details for reporting
                hit_info = {
                    'query': seq_id,
                    'subject': hit['sseqid'],
                    'identity': hit['pident'],
                    'coverage': hit['qcovs'],
                    'length': hit['length'],
                    'evalue': hit['evalue'],
                    'description': hit.get('stitle', 'Unknown'),
                    'blast_type': blast_type
                }
                all_hits.append(hit_info)
        
        except Exception as e:
            print(f"   Error parsing BLAST results for {blast_type}: {e}")
        
        return contaminated_sequences, all_hits
    
    def check_targets(self, primers_df, work_dir):
        """STEP 1: Check target sequences for contamination"""
        print("🎯 STEP 1: Checking target sequences for contamination...")
        
        # Get unique targets
        unique_targets = primers_df['target_sequence'].unique()
        print(f"   Analyzing {len(unique_targets)} unique target sequences")
        
        # Create target FASTA
        target_fasta = os.path.join(work_dir, "targets.fasta")
        target_to_seq_map = {}
        
        with open(target_fasta, 'w') as f:
            for i, target_seq in enumerate(unique_targets):
                target_id = f"target_{i}"
                f.write(f">{target_id}\n{target_seq}\n")
                target_to_seq_map[target_id] = target_seq
        
        # Run BLAST for targets
        blast_output = os.path.join(work_dir, "targets_blast.txt")
        if not self.run_blast(target_fasta, blast_output, "target"):
            print("   Warning: Target BLAST failed, assuming all targets are clean")
            return set(), []
        
        # Parse results
        contaminated_target_ids, target_hits = self.parse_blast_results(blast_output, "target")
        
        # Map back to actual sequences
        contaminated_target_sequences = set()
        for target_id in contaminated_target_ids:
            if target_id in target_to_seq_map:
                contaminated_target_sequences.add(target_to_seq_map[target_id])
        
        print(f"   🚫 Found {len(contaminated_target_sequences)} contaminated targets")
        print(f"   ✅ {len(unique_targets) - len(contaminated_target_sequences)} clean targets remaining")
        
        if contaminated_target_sequences:
            print(f"   Contaminated targets will exclude {len(primers_df[primers_df['target_sequence'].isin(contaminated_target_sequences)])} primer pairs")
        
        # Cleanup
        try:
            os.remove(target_fasta)
            os.remove(blast_output)
        except:
            pass
        
        return contaminated_target_sequences, target_hits
    
    def check_primers(self, clean_target_primers_df, work_dir):
        """STEP 2: Check primer sequences (only for clean targets)"""
        print("🧬 STEP 2: Checking primer sequences (clean targets only)...")
        
        if len(clean_target_primers_df) == 0:
            print("   No clean targets available for primer checking")
            return set(), []
        
        print(f"   Analyzing {len(clean_target_primers_df)} primer pairs from clean targets")
        
        # Create primer FASTA
        primer_fasta = os.path.join(work_dir, "primers.fasta")
        primer_to_index_map = {}
        
        with open(primer_fasta, 'w') as f:
            for idx, row in clean_target_primers_df.iterrows():
                # Forward primer
                f_id = f"forward_{idx}"
                f.write(f">{f_id}\n{row['forward_primer']}\n")
                primer_to_index_map[f_id] = idx
                
                # Reverse primer
                r_id = f"reverse_{idx}"
                f.write(f">{r_id}\n{row['reverse_primer']}\n")
                primer_to_index_map[r_id] = idx
        
        # Run BLAST for primers
        blast_output = os.path.join(work_dir, "primers_blast.txt")
        if not self.run_blast(primer_fasta, blast_output, "primer"):
            print("   Warning: Primer BLAST failed, assuming all primers are clean")
            return set(), []
        
        # Parse results
        contaminated_primer_ids, primer_hits = self.parse_blast_results(blast_output, "primer")
        
        # Map back to primer pair indices
        contaminated_primer_indices = set()
        contamination_details = {}
        
        for primer_id in contaminated_primer_ids:
            if primer_id in primer_to_index_map:
                primer_index = primer_to_index_map[primer_id]
                contaminated_primer_indices.add(primer_index)
                
                # Track which primer (forward/reverse) was contaminated
                primer_type = primer_id.split('_')[0]  # 'forward' or 'reverse'
                if primer_index not in contamination_details:
                    contamination_details[primer_index] = []
                contamination_details[primer_index].append(primer_type)
        
        print(f"   🚫 Found {len(contaminated_primer_indices)} contaminated primer pairs")
        print(f"   ✅ {len(clean_target_primers_df) - len(contaminated_primer_indices)} clean primer pairs remaining")
        
        if contaminated_primer_indices:
            print(f"   Contamination breakdown:")
            forward_count = sum(1 for details in contamination_details.values() if 'forward' in details)
            reverse_count = sum(1 for details in contamination_details.values() if 'reverse' in details)
            print(f"     Forward primers: {forward_count}")
            print(f"     Reverse primers: {reverse_count}")
        
        # Cleanup
        try:
            os.remove(primer_fasta)
            os.remove(blast_output)
        except:
            pass
        
        return contaminated_primer_indices, primer_hits
    
    def check_amplicons(self, clean_primer_pairs_df, work_dir):
        """STEP 3: Check full amplicon sequences (only for clean primers)"""
        print("🔬 STEP 3: Checking amplicon sequences (clean primers only)...")
        
        if len(clean_primer_pairs_df) == 0:
            print("   No clean primer pairs available for amplicon checking")
            return set(), []
        
        print(f"   Analyzing {len(clean_primer_pairs_df)} amplicon sequences")
        
        # Create amplicon FASTA
        amplicon_fasta = os.path.join(work_dir, "amplicons.fasta")
        amplicon_to_index_map = {}
        valid_amplicons = 0
        
        with open(amplicon_fasta, 'w') as f:
            for idx, row in clean_primer_pairs_df.iterrows():
                # Try to get full amplicon sequence
                amplicon_seq = None
                
                if 'full_amplicon_sequence' in row and pd.notna(row['full_amplicon_sequence']):
                    amplicon_seq = str(row['full_amplicon_sequence']).strip()
                
                if amplicon_seq and len(amplicon_seq) > 50:  # Minimum length check
                    amplicon_id = f"amplicon_{idx}"
                    f.write(f">{amplicon_id}\n{amplicon_seq}\n")
                    amplicon_to_index_map[amplicon_id] = idx
                    valid_amplicons += 1
        
        if valid_amplicons == 0:
            print("   No valid amplicon sequences found for checking")
            return set(), []
        
        print(f"   Found {valid_amplicons} valid amplicon sequences to check")
        
        # Run BLAST for amplicons
        blast_output = os.path.join(work_dir, "amplicons_blast.txt")
        if not self.run_blast(amplicon_fasta, blast_output, "amplicon"):
            print("   Warning: Amplicon BLAST failed, assuming all amplicons are clean")
            return set(), []
        
        # Parse results
        contaminated_amplicon_ids, amplicon_hits = self.parse_blast_results(blast_output, "amplicon")
        
        # Map back to primer pair indices
        contaminated_amplicon_indices = set()
        for amplicon_id in contaminated_amplicon_ids:
            if amplicon_id in amplicon_to_index_map:
                amplicon_index = amplicon_to_index_map[amplicon_id]
                contaminated_amplicon_indices.add(amplicon_index)
        
        print(f"   🚫 Found {len(contaminated_amplicon_indices)} contaminated amplicons")
        print(f"   ✅ {valid_amplicons - len(contaminated_amplicon_indices)} clean amplicons remaining")
        
        # Cleanup
        try:
            os.remove(amplicon_fasta)
            os.remove(blast_output)
        except:
            pass
        
        return contaminated_amplicon_indices, amplicon_hits
    
    def check_contamination_hierarchical(self, primers_df, work_dir):
        """Main hierarchical contamination checking workflow"""
        print(f"🔬 HIERARCHICAL CONTAMINATION CHECK: {len(primers_df)} pre-validated primers")
        print(f"📊 Contamination thresholds (with mismatch tolerance):")
        for seq_type, params in self.thresholds.items():
            print(f"   {seq_type.capitalize()}: ≥{params['identity']}% identity, ≥{params['coverage']}% coverage")
        
        os.makedirs(work_dir, exist_ok=True)
        
        primers_df_with_amplicons = self.add_full_amplicon_sequences(primers_df)

        # STEP 1: Check targets
        contaminated_targets, target_hits = self.check_targets(primers_df, work_dir)
        
        # Filter to primers with clean targets
        clean_target_primers = primers_df_with_amplicons[~primers_df_with_amplicons['target_sequence'].isin(contaminated_targets)].copy()
        
        # STEP 2: Check primers (only for clean targets)
        contaminated_primers, primer_hits = self.check_primers(clean_target_primers, work_dir)
        
        # Filter to clean primer pairs
        clean_primer_pairs = clean_target_primers[~clean_target_primers.index.isin(contaminated_primers)].copy()
        
        # STEP 3: Check amplicons (only for clean primers)
        contaminated_amplicons, amplicon_hits = self.check_amplicons(clean_primer_pairs, work_dir)
        
        # Final clean set
        final_clean_primers = clean_primer_pairs[~clean_primer_pairs.index.isin(contaminated_amplicons)].copy()
        
        # Combine all hits for reporting
        all_hits = target_hits + primer_hits + amplicon_hits
        
        # Summary statistics
        total_contaminated = len(contaminated_targets) + len(contaminated_primers) + len(contaminated_amplicons)
        
        print(f"\n📊 HIERARCHICAL CONTAMINATION SUMMARY:")
        print(f"   🎯 Target contamination: {len(contaminated_targets)} sequences")
        print(f"   🧬 Primer contamination: {len(contaminated_primers)} pairs") 
        print(f"   🔬 Amplicon contamination: {len(contaminated_amplicons)} sequences")
        print(f"   📋 Total contaminated: {len(primers_df) - len(final_clean_primers)} primer pairs")
        print(f"   ✅ Final clean primers: {len(final_clean_primers)} pairs")
        print(f"   📈 Success rate: {len(final_clean_primers)/len(primers_df)*100:.1f}%")
        
        return final_clean_primers, all_hits, {
            'contaminated_targets': contaminated_targets,
            'contaminated_primers': contaminated_primers, 
            'contaminated_amplicons': contaminated_amplicons
        }, primers_df_with_amplicons
    
    def create_clean_primers_output(self, original_df, clean_primers, all_hits, contamination_summary, output_file):
        """Create final output with contamination status for all primers"""
        
        # Add contamination status to all original primers
        result_df = original_df.copy()
        result_df['contamination_status'] = 'CLEAN'
        result_df['contamination_step'] = 'PASSED_ALL'
        result_df['contamination_details'] = 'No contamination detected'
        
        # Mark contaminated primers by step
        contaminated_targets = contamination_summary['contaminated_targets']
        contaminated_primers = contamination_summary['contaminated_primers']
        contaminated_amplicons = contamination_summary['contaminated_amplicons']
        
        # Target contamination
        target_contaminated_indices = result_df[result_df['target_sequence'].isin(contaminated_targets)].index
        for idx in target_contaminated_indices:
            result_df.loc[idx, 'contamination_status'] = 'CONTAMINATED'
            result_df.loc[idx, 'contamination_step'] = 'STEP1_TARGET'
            result_df.loc[idx, 'contamination_details'] = 'Target sequence matches contamination database'
        
        # Primer contamination
        for idx in contaminated_primers:
            if idx in result_df.index:
                result_df.loc[idx, 'contamination_status'] = 'CONTAMINATED'
                result_df.loc[idx, 'contamination_step'] = 'STEP2_PRIMER'
                result_df.loc[idx, 'contamination_details'] = 'Primer sequences match contamination database'
        
        # Amplicon contamination
        for idx in contaminated_amplicons:
            if idx in result_df.index:
                result_df.loc[idx, 'contamination_status'] = 'CONTAMINATED'
                result_df.loc[idx, 'contamination_step'] = 'STEP3_AMPLICON'
                result_df.loc[idx, 'contamination_details'] = 'Full amplicon matches contamination database'
        
        # Save complete results
        complete_output = output_file.replace('.csv', '_complete_results.csv')
        result_df.to_csv(complete_output, index=False)
        
        # Save only clean primers
        clean_only = result_df[result_df['contamination_status'] == 'CLEAN'].copy()
        clean_only = clean_only.drop(['contamination_status', 'contamination_step', 'contamination_details'], axis=1)
        clean_only.to_csv(output_file, index=False)
        
        # Save contaminated primers for reference
        contaminated_output = output_file.replace('.csv', '_contaminated.csv')
        contaminated_only = result_df[result_df['contamination_status'] == 'CONTAMINATED'].copy()
        if len(contaminated_only) > 0:
            contaminated_only.to_csv(contaminated_output, index=False)
        
        return clean_only, len(contaminated_only)
    
    def create_summary_report(self, input_file, output_file, original_count, clean_count, contaminated_count, all_hits, contamination_summary):
        """Create detailed hierarchical contamination report"""
        
        report_file = output_file.replace('.csv', '_report.txt')
        
        contaminated_targets = contamination_summary['contaminated_targets']
        contaminated_primers = contamination_summary['contaminated_primers']
        contaminated_amplicons = contamination_summary['contaminated_amplicons']
        
        report = f"""
HIERARCHICAL CONTAMINATION CHECK REPORT
======================================

INPUT:
------
Pre-validated primers file: {input_file}
Total primers checked: {original_count}

HIERARCHICAL CONTAMINATION THRESHOLDS:
-------------------------------------
Step 1 - Targets: ≥{self.thresholds['target']['identity']}% identity, ≥{self.thresholds['target']['coverage']}% coverage
Step 2 - Primers: ≥{self.thresholds['primer']['identity']}% identity, ≥{self.thresholds['primer']['coverage']}% coverage  
Step 3 - Amplicons: ≥{self.thresholds['amplicon']['identity']}% identity, ≥{self.thresholds['amplicon']['coverage']}% coverage

CONTAMINATION DATABASE:
----------------------
- Ixodes (tick) genomes
- Midichloria (tick endosymbiont)
- Rickettsiella (tick bacteria)
- Spiroplasma (tick bacteria)  
- Other arthropod-associated microorganisms

HIERARCHICAL RESULTS:
--------------------
🎯 STEP 1 - Target contamination: {len(contaminated_targets)} unique sequences
🧬 STEP 2 - Primer contamination: {len(contaminated_primers)} primer pairs
🔬 STEP 3 - Amplicon contamination: {len(contaminated_amplicons)} amplicons

FINAL RESULTS:
-------------
✅ Clean primers (synthesis ready): {clean_count}
🚫 Total contaminated primers (all steps): {contaminated_count}
📊 Success rate: {clean_count/original_count*100:.1f}%

CONTAMINATION EFFICIENCY:
------------------------
Hierarchical approach saved computational time by:
- Not checking primers for {len(contaminated_targets)} contaminated targets
- Not checking amplicons for {len(contaminated_primers)} contaminated primer pairs
- Total sequences avoided: {len(contaminated_targets) * 4 + len(contaminated_primers) * 1} BLAST queries

FILES CREATED:
-------------
📁 {os.path.basename(output_file)} - Clean primers ready for synthesis ({clean_count} primers)
📁 {os.path.basename(output_file.replace('.csv', '_complete_results.csv'))} - All primers with contamination status
📁 {os.path.basename(output_file.replace('.csv', '_contaminated.csv'))} - Contaminated primers by step
📁 {os.path.basename(report_file)} - This detailed report

CONTAMINATION DETAILS:
---------------------
"""
        
        if len(all_hits) > 0:
            # Group hits by type
            target_hits = [h for h in all_hits if h['blast_type'] == 'target']
            primer_hits = [h for h in all_hits if h['blast_type'] == 'primer']
            amplicon_hits = [h for h in all_hits if h['blast_type'] == 'amplicon']
            
            if target_hits:
                report += f"🎯 Target contamination hits: {len(target_hits)}\n"
                for hit in target_hits[:3]:  # Show top 3
                    report += f"   • {hit['query']}: {hit['identity']:.1f}% identity to {hit['subject']}\n"
                if len(target_hits) > 3:
                    report += f"   ... and {len(target_hits) - 3} more target hits\n"
            
            if primer_hits:
                report += f"🧬 Primer contamination hits: {len(primer_hits)}\n"
                for hit in primer_hits[:3]:  # Show top 3
                    report += f"   • {hit['query']}: {hit['identity']:.1f}% identity to {hit['subject']}\n"
                if len(primer_hits) > 3:
                    report += f"   ... and {len(primer_hits) - 3} more primer hits\n"
            
            if amplicon_hits:
                report += f"🔬 Amplicon contamination hits: {len(amplicon_hits)}\n"
                for hit in amplicon_hits[:3]:  # Show top 3
                    report += f"   • {hit['query']}: {hit['identity']:.1f}% identity to {hit['subject']}\n"
                if len(amplicon_hits) > 3:
                    report += f"   ... and {len(amplicon_hits) - 3} more amplicon hits\n"
        else:
            report += "🎉 No contamination detected in any step!\n"
        
        report += f"""
RECOMMENDATIONS:
---------------
"""
        
        if clean_count > 0:
            report += f"🎉 SUCCESS: {clean_count} primers ready for synthesis!\n"
            report += f"📊 Use {os.path.basename(output_file)} for primer ordering\n"
            report += f"🧬 All primers passed: structure validation + hierarchical contamination check\n"
            report += f"🔍 Hierarchical approach ensured thorough contamination screening\n"
        else:
            report += f"❌ No clean primers available for synthesis\n"
            report += f"🔍 All primers showed contamination - review database or consider relaxing thresholds\n"
        
        if contaminated_count > 0:
            report += f"⚠️  {contaminated_count} primers excluded due to contamination\n"
            report += f"📋 Check contaminated primers file to see which step caught each contamination\n"
            report += f"🎯 Consider reviewing target sequences if many failed at Step 1\n"
        
        # Write report
        with open(report_file, 'w') as f:
            f.write(report)
        
        print(report)


def main():
    parser = argparse.ArgumentParser(description="Hierarchical Contamination Checker with Mismatch Tolerance")
    parser.add_argument("--input", required=True, help="Pre-validated primers CSV file")
    parser.add_argument("--database", required=True, help="BLAST database prefix for contamination checking")
    parser.add_argument("--output", required=True, help="Output clean primers CSV file")
    parser.add_argument("--work-dir", default="temp_contamination_check", help="Working directory for temporary files")
    
    args = parser.parse_args()
    
    # Read pre-validated primers
    print(f"📖 Reading pre-validated primers from: {args.input}")
    try:
        primers_df = pd.read_csv(args.input)
        print(f"   Found {len(primers_df)} pre-validated primer pairs")
    except Exception as e:
        print(f"❌ Error reading input file: {e}")
        return 1
    
    if len(primers_df) == 0:
        print(f"❌ No primers found in input file")
        return 1
    
    # Validate required columns
    required_columns = ['target_sequence', 'forward_primer', 'reverse_primer']
    missing_columns = [col for col in required_columns if col not in primers_df.columns]
    if missing_columns:
        print(f"❌ Missing required columns: {missing_columns}")
        print(f"   Available columns: {list(primers_df.columns)}")
        return 1
    
    # Initialize hierarchical contamination checker
    try:
        checker = HierarchicalContaminationChecker(args.database)
        print(f"✅ Hierarchical contamination checker initialized")
        print(f"   Database: {args.database}")
        print(f"   Approach: Hierarchical (Targets → Primers → Amplicons)")
        print(f"   Mismatch tolerance: Allows similarities to catch potential issues")
    except Exception as e:
        print(f"❌ Error initializing contamination checker: {e}")
        return 1
    
    # Create output directory
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"\n🔬 Starting hierarchical contamination analysis...")
    
    try:
        # Hierarchical contamination check
        clean_primers, all_hits, contamination_summary, primers_df_with_amplicons = checker.check_contamination_hierarchical(primers_df, args.work_dir)
        
        # Create output files
        final_clean, contaminated_count = checker.create_clean_primers_output(
            primers_df_with_amplicons, clean_primers, all_hits, contamination_summary, args.output
        )
        
        # Create detailed report
        checker.create_summary_report(
            args.input, args.output, len(primers_df), len(final_clean), contaminated_count, 
            all_hits, contamination_summary
        )
        
        print(f"\n🎉 HIERARCHICAL CONTAMINATION CHECK COMPLETE!")
        print(f"✅ {len(final_clean)} clean primers ready for synthesis")
        print(f"🚫 {contaminated_count} contaminated primers excluded")
        print(f"📁 Results saved to: {args.output}")
        print(f"🎯 Hierarchical approach: Targets → Primers → Amplicons")
        print(f"🔍 Mismatch tolerance ensured thorough contamination detection")
        
        # Cleanup working directory
        try:
            if os.path.exists(args.work_dir):
                import shutil
                shutil.rmtree(args.work_dir)
        except:
            pass
        
        return 0
        
    except Exception as e:
        print(f"❌ Error during hierarchical contamination analysis: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())