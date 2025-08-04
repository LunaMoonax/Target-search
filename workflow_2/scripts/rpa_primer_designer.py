#!/usr/bin/env python3
"""
RPA Primer Designer - EXHAUSTIVE SEARCH VERSION
"""
import sys
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re

# --- PAGALBINĖS FUNKCIJOS ---
def extract_genome_name(seq_id):
    return seq_id.split(';')[0] if ';' in seq_id else seq_id

def parse_alignment_blocks(fasta_file):
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        print("KLAIDA: Įvesties FASTA failas yra tuščias.")
        return
    
    first_genome_name = extract_genome_name(records[0].id)
    current_block = []
    for rec in records:
        genome_name = extract_genome_name(rec.id)
        if genome_name == first_genome_name and current_block:
            yield current_block
            current_block = []
        current_block.append(rec)
    if current_block:
        yield current_block

def map_ungapped_to_gapped(sequence):
    mapping = {}
    ungapped_idx = 0
    for gapped_idx, char in enumerate(sequence):
        if char != '-':
            mapping[ungapped_idx] = gapped_idx
            ungapped_idx += 1
    return mapping

# --- PRADMENŲ KŪRIMO KLASĖ SU PILNU VALIDAVIMU ---
class RPAPrimerDesigner:
    def __init__(self):
        # --- VALDYMO PULTAS: Čia galite keisti reikalavimus ---
        self.P_MIN_LEN, self.P_MAX_LEN = 30, 35
        self.P_MIN_GC, self.P_MAX_GC = 30, 70
        self.A_MIN_SIZE, self.A_MAX_SIZE = 200, 300
        self.A_MIN_GC, self.A_MAX_GC = 35, 60
        self.CONSERVATION_THRESHOLD = 0.90

    @staticmethod
    def calculate_conservation(sequences):
        if len(sequences) < 2: return 1.0
        total_identity, num_pairs = 0.0, 0
        for i in range(len(sequences)):
            for j in range(i + 1, len(sequences)):
                seq1 = sequences[i].replace('-', '')
                seq2 = sequences[j].replace('-', '')
                min_len = min(len(seq1), len(seq2))
                if min_len == 0: continue
                matches = sum(1 for k in range(min_len) if seq1[k] == seq2[k])
                total_identity += matches / min_len
                num_pairs += 1
        return total_identity / num_pairs if num_pairs > 0 else 0.0

    def _validate_primer(self, primer_seq):
        if not primer_seq: return False
        if not (self.P_MIN_LEN <= len(primer_seq) <= self.P_MAX_LEN): return False
        gc = gc_fraction(primer_seq) * 100
        if not (self.P_MIN_GC <= gc <= self.P_MAX_GC): return False
        if re.search(r'([ATGC])\1{4,}', primer_seq): return False
        return True

    def _validate_amplicon(self, amplicon_seq):
        if not amplicon_seq: return False
        size = len(amplicon_seq)
        if not (self.A_MIN_SIZE <= size <= self.A_MAX_SIZE): return False
        gc = gc_fraction(amplicon_seq) * 100
        if not (self.A_MIN_GC <= gc <= self.A_MAX_GC): return False
        return True

    def design(self, block, target_start, target_end, target_info, stats):
        primer_pairs = []
        reference_seq_gapped = str(block[0].seq)
        reference_seq_ungapped = reference_seq_gapped.replace('-', '')
        target_len_ungapped = target_end - target_start
        
        ungapped_map = map_ungapped_to_gapped(reference_seq_gapped)

        # Ciklas per visus galimus amplikono dydžius
        for amp_len in range(self.A_MIN_SIZE, self.A_MAX_SIZE + 1):
            # Ciklas per visus galimus Forward pradmens ilgius
            for f_len in range(self.P_MIN_LEN, self.P_MAX_LEN + 1):
                # Ciklas per visus galimus Reverse pradmens ilgius
                for r_len in range(self.P_MIN_LEN, self.P_MAX_LEN + 1):
                    
                    flanking_needed = amp_len - target_len_ungapped - f_len - r_len
                    if flanking_needed < 0: continue
                    
                    # Ciklas per visus galimus kairės srities ilgius
                    for left_flank_len in range(flanking_needed + 1):
                        right_flank_len = flanking_needed - left_flank_len
                        stats["Total Pairs Tested"] += 1
                        
                        f_end_pos = target_start - left_flank_len
                        f_start_pos = f_end_pos - f_len
                        
                        r_start_pos = target_end + right_flank_len
                        r_end_pos = r_start_pos + r_len

                        if f_start_pos < 0 or r_end_pos > len(reference_seq_ungapped): continue
                        
                        try:
                            f_start_gapped = ungapped_map[f_start_pos]
                            f_end_gapped = ungapped_map[f_end_pos - 1] + 1
                            r_start_gapped = ungapped_map[r_start_pos]
                            r_end_gapped = ungapped_map[r_end_pos - 1] + 1
                        except KeyError:
                            continue

                        forward_aln = [str(rec.seq[f_start_gapped:f_end_gapped]) for rec in block]
                        if self.calculate_conservation(forward_aln) < self.CONSERVATION_THRESHOLD:
                            stats["Rejected: Primer Conservation"] += 1; continue

                        reverse_aln = [str(rec.seq[r_start_gapped:r_end_gapped]) for rec in block]
                        if self.calculate_conservation(reverse_aln) < self.CONSERVATION_THRESHOLD:
                            stats["Rejected: Primer Conservation"] += 1; continue
                            
                        forward_primer = forward_aln[0].replace('-', '')
                        reverse_template = reverse_aln[0].replace('-', '')
                        
                        if not self._validate_primer(forward_primer) or not self._validate_primer(str(Seq(reverse_template).reverse_complement())):
                            stats["Rejected: Primer Properties (GC/Length)"] += 1; continue
                        
                        amplicon = reference_seq_ungapped[f_start_pos:r_end_pos]
                        if not self._validate_amplicon(amplicon):
                            stats["Rejected: Amplicon Properties (Size/GC)"] += 1; continue

                        stats["Valid Designs Found"] += 1
                        design = {"forward_primer": forward_primer, "reverse_primer": str(Seq(reverse_template).reverse_complement()),
                                  "amplicon_sequence": amplicon, "amplicon_size": len(amplicon)}
                        design.update(target_info)
                        primer_pairs.append(design)
                        #if len(primer_pairs) >= 1: return primer_pairs, stats # Grįžtame po pirmo sėkmingo
        return primer_pairs, stats

def main():
    parser = argparse.ArgumentParser(description="Kuriami RPA pradmenys su pilnu validavimu.")
    parser.add_argument("targets_csv")
    parser.add_argument("alignment_fasta")
    parser.add_argument("output_csv")
    args = parser.parse_args()

    try:
        targets_df = pd.read_csv(args.targets_csv)
    except FileNotFoundError:
        print(f"KLAIDA: Nerastas taikinių failas: {args.targets_csv}"); sys.exit(1)
        
    designer = RPAPrimerDesigner()
    all_designs = []

    alignment_blocks = {i: block for i, block in enumerate(parse_alignment_blocks(args.alignment_fasta), 1)}
    unique_target_positions = targets_df.drop_duplicates(subset=['block_id', 'start_pos'])
    print(f"Bus apdorojama {len(unique_target_positions)} unikalių taikinių pozicijų iš {len(alignment_blocks)} blokų.")

    total_stats = { "Total Targets Processed": len(unique_target_positions), "Rejected: Target Not Mapped": 0, "Total Pairs Tested": 0,
                    "Rejected: Primer Conservation": 0, "Rejected: Primer Properties (GC/Length)": 0, 
                    "Rejected: Amplicon Properties (Size/GC)": 0, "Valid Designs Found": 0 }

    for i, target_row in unique_target_positions.iterrows():
        block_id, start_pos = target_row['block_id'], target_row['start_pos']
        block = alignment_blocks.get(block_id)
        if not block: continue
        end_pos = start_pos + len(target_row['sequence_variant'])
        
        designs, stats = designer.design(block, start_pos, end_pos, target_row.to_dict(), total_stats)
        if designs: all_designs.extend(designs)

    print("\n--- DIAGNOSTIKOS SUVESTINĖ ---")
    for reason, count in total_stats.items():
        print(f"  {reason}: {count}")
    
    if not all_designs:
        print("\nTinkamų pradmenų porų nerasta."); pd.DataFrame().to_csv(args.output_csv, index=False)
        return
        
    results_df = pd.DataFrame(all_designs).drop_duplicates(subset=["forward_primer", "reverse_primer"])
    results_df.to_csv(args.output_csv, index=False)
    print(f"\nSukurta {len(results_df)} unikalių pradmenų porų. Rezultatai išsaugoti: {args.output_csv}")

if __name__ == "__main__":
    main()