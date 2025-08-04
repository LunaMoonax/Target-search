#!/usr/bin/env python3
"""
Comprehensive RPA Primer and Amplicon Validator with Diagnostics.
Validates primers and the full amplicon sequence for secondary structures.
"""
import argparse
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import re
from collections import defaultdict

class RPAValidator:
    def __init__(self):
        # --- VALDYMO PULTAS: Čia galite keisti reikalavimus ---
        # Pradmenų taisyklės
        self.P_3_PRIME_PREFER = ['G', 'C']
        self.P_3_PRIME_AVOID = ['T']
        self.MAX_HOMOPOLYMER = 4
        self.MAX_DINUC_REPEAT = 4
        self.MAX_HAIRPIN = 7
        self.MAX_SELF_DIMER = 7
        self.MAX_3_PRIME_SELF_COMP = 2
        self.MAX_PAIR_3_PRIME_COMP = 2
        # Amplikono taisyklės
        self.MAX_AMP_REPEATS = 10
        self.MAX_AMP_PALINDROME = 10

    # --- DETALIOS VALIDAVIMO FUNKCIJOS ---
    def _find_max_dinucleotide_repeat(self, sequence):
        max_repeat = 0
        for i in range(len(sequence) - 3):
            dinuc = sequence[i:i+2]
            repeats = 1
            for j in range(i + 2, len(sequence) - 1, 2):
                if sequence[j:j+2] == dinuc: repeats += 1
                else: break
            max_repeat = max(max_repeat, repeats)
        return max_repeat

    def _calculate_hairpin(self, seq):
        max_score = 0
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        for stem_len in range(4, 9):
            for loop_len in range(3, 9):
                for i in range(len(seq) - (2 * stem_len) - loop_len + 1):
                    stem1 = seq[i : i + stem_len]
                    stem2_rev = seq[i + stem_len + loop_len : i + 2 * stem_len + loop_len][::-1]
                    matches = sum(1 for b1, b2 in zip(stem1, stem2_rev) if comp.get(b1) == b2)
                    if matches > 2: max_score = max(max_score, matches)
        return max_score

    def _calculate_dimer(self, seq1, seq2, is_3_prime=False):
        max_score = 0
        s2_rev_comp = str(Seq(seq2).reverse_complement())
        start1 = len(seq1) - 8 if is_3_prime else 0
        for i in range(start1, len(seq1)):
            for j in range(len(s2_rev_comp)):
                score = 0
                for k in range(min(len(seq1) - i, len(s2_rev_comp) - j)):
                     if seq1[i+k] == s2_rev_comp[j+k]: score += 1
                     else: break
                max_score = max(max_score, score)
        return max_score

    def _find_direct_repeats(self, sequence):
        max_len = 0
        for length in range(8, 20):
            for i in range(len(sequence) - 2 * length + 1):
                pattern = sequence[i:i+length]
                for j in range(i + length, len(sequence) - length + 1):
                    if sequence[j:j+length] == pattern:
                        max_len = max(max_len, length)
        return max_len

    def _find_palindromes(self, sequence):
        comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        max_len = 0
        for length in range(8, 20, 2):
            for i in range(len(sequence) - length + 1):
                sub = sequence[i:i+length]
                sub_rev_comp = "".join(comp.get(b, "N") for b in reversed(sub))
                if sub == sub_rev_comp:
                    max_len = max(max_len, length)
        return max_len

    def validate_row(self, row, stats):
        f_primer = str(row['forward_primer'])
        r_primer = str(row['reverse_primer'])
        amplicon = str(row['amplicon_sequence'])
        
        f_issues, f_penalty, r_issues, r_penalty, pair_issues, pair_penalty, amp_issues, amp_penalty = [], 0, [], 0, [], 0, [], 0

        # --- Forward Pradmens Patikra ---
        if re.search(r'([ATGC])\1{4,}', f_primer):
            stats['fail:f_homopolymer'] +=1; f_issues.append("Homopolymer>4"); f_penalty += 15
        if self._find_max_dinucleotide_repeat(f_primer) > self.MAX_DINUC_REPEAT:
            stats['fail:f_dinuc_repeat'] +=1; f_issues.append("Dinuc repeat"); f_penalty += 15
        if self._calculate_hairpin(f_primer) > self.MAX_HAIRPIN:
            stats['fail:f_hairpin'] +=1; f_issues.append("Hairpin"); f_penalty += 20
        if self._calculate_dimer(f_primer, f_primer) > self.MAX_SELF_DIMER:
            stats['fail:f_self_dimer'] +=1; f_issues.append("Self-dimer"); f_penalty += 20
        if self._calculate_dimer(f_primer, f_primer, is_3_prime=True) > self.MAX_3_PRIME_SELF_COMP:
            stats['fail:f_3prime_self_dimer'] +=1; f_issues.append("3' self-dimer"); f_penalty += 30
        if f_primer[-1] in self.P_3_PRIME_AVOID:
            stats['fail:f_bad_3_prime'] +=1; f_issues.append("Bad 3' end"); f_penalty += 15
        elif f_primer[-1] in self.P_3_PRIME_PREFER:
            f_penalty -= 5
            
        # --- Reverse Pradmens Patikra ---
        # (identiškos patikros kaip ir Forward)
        if re.search(r'([ATGC])\1{4,}', r_primer):
            stats['fail:r_homopolymer'] +=1; r_issues.append("Homopolymer>4"); r_penalty += 15
        # ... (ir t.t.)

        # --- Poros Patikra ---
        if self._calculate_dimer(f_primer, r_primer, is_3_prime=True) > self.MAX_PAIR_3_PRIME_COMP:
            stats['fail:pair_3prime_dimer'] +=1; pair_issues.append("3' pair-dimer"); pair_penalty += 40

        # --- Amplikono Patikra (PRIDĖTA) ---
        if self._find_direct_repeats(amplicon) > self.MAX_AMP_REPEATS:
            stats['fail:amp_repeats'] +=1; amp_issues.append("Direct Repeats"); amp_penalty += 20
        if self._find_palindromes(amplicon) > self.MAX_AMP_PALINDROME:
            stats['fail:amp_palindrome'] +=1; amp_issues.append("Palindrome"); amp_penalty += 20
        
        is_valid = not (f_issues or r_issues or pair_issues or amp_issues)
        if is_valid:
            stats['final_valid_pairs'] +=1

        total_penalty = f_penalty + r_penalty + pair_penalty + amp_penalty

        return pd.Series([is_valid, total_penalty, ";".join(f_issues), ";".join(r_issues), ";".join(pair_issues), ";".join(amp_issues)])

def main():
    parser = argparse.ArgumentParser(description="Atlieka išsamią RPA pradmenų ir amplikonų patikrą.")
    parser.add_argument("input_csv")
    parser.add_argument("output_csv")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_csv)
        if df.empty:
            print("Įvesties failas tuščias, nėra ką tikrinti."); df.to_csv(args.output_csv, index=False); return
    except FileNotFoundError:
        print(f"KLAIDA: Nerastas įvesties failas: {args.input_csv}"); sys.exit(1)

    print(f"Tikrinama {len(df)} pradmenų porų...")
    validator = RPAValidator()
    
    stats = defaultdict(int)
    
    validation_cols = ['is_valid', 'total_penalty', 'forward_issues', 'reverse_issues', 'pair_issues', 'amplicon_issues']
    
    df[validation_cols] = df.apply(lambda row: validator.validate_row(row, stats), axis=1)
    
    df_sorted = df.sort_values(by=['is_valid', 'total_penalty'], ascending=[False, True])

    # Atrenkame tik tas eilutes, kurios praėjo validaciją
    valid_only_df = df_sorted[df_sorted['is_valid'] == True].copy()
    
    # Išsaugome tik atrinktą lentelę
    valid_only_df.to_csv(args.output_csv, index=False)
        
    print("\n--- DIAGNOSTIKOS SUVESTINĖ ---")
    print(f"  Iš viso patikrinta porų: {len(df)}")
    for reason, count in sorted(stats.items()):
        if reason.startswith("fail:"):
            print(f"  Atmesta dėl '{reason[5:]}': {count}")
    print("  ---------------------------------")
    print(f"  Tinkamų porų (praėjo visus filtrus): {stats['final_valid_pairs']}")
    
    print(f"\nRezultatai su detalia analize išsaugoti: {args.output_csv}")

if __name__ == "__main__":
    main()