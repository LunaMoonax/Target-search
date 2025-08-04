#!/usr/bin/env python3
"""
Generates a final report for top primer candidates by finding the primer
positions in the alignment and extracting all corresponding amplicon variants.
"""
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

def extract_genome_name(seq_id):
    """Extracts the genome name from a FASTA header."""
    return seq_id.split(';')[0] if ';' in seq_id else seq_id

def parse_alignment_blocks(fasta_file):
    """Parses a concatenated multi-FASTA alignment file into blocks."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records: return
    
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
    """Creates a map from ungapped coordinates to gapped coordinates."""
    mapping = {}
    ungapped_idx = 0
    for gapped_idx, char in enumerate(sequence):
        if char != '-':
            mapping[ungapped_idx] = gapped_idx
            ungapped_idx += 1
    return mapping

def main():
    parser = argparse.ArgumentParser(description="Generates a final report with all amplicon variants.")
    parser.add_argument("--primers", required=True, help="Input CSV file with top primer pairs (e.g., top_primers.csv).")
    parser.add_argument("--alignment", required=True, help="The partitioned alignment FASTA file.")
    parser.add_argument("--output", required=True, help="Output CSV file for the final report.")
    args = parser.parse_args()

    try:
        primers_df = pd.read_csv(args.primers)
    except FileNotFoundError:
        print(f"KLAIDA: Nerastas pradmenų failas: {args.primers}"); return

    print("Analizuojami sulyginimo blokai...")
    alignment_blocks = {i: block for i, block in enumerate(parse_alignment_blocks(args.alignment), 1)}
    
    new_column_data = []

    for _, row in primers_df.iterrows():
        block_id = int(row['block_id'])
        f_primer = str(row['forward_primer'])
        r_primer = str(row['reverse_primer'])
        
        block = alignment_blocks.get(block_id)
        if not block:
            new_column_data.append("Error: Block not found")
            continue

        # --- Logika koordinačių radimui ---
        ref_gapped = str(block[0].seq)
        ref_ungapped = ref_gapped.replace('-', '')
        ungapped_map = map_ungapped_to_gapped(ref_gapped)

        f_start_ungapped = ref_ungapped.find(f_primer)
        r_template = str(Seq(r_primer).reverse_complement())
        r_start_ungapped = ref_ungapped.find(r_template, f_start_ungapped)

        if f_start_ungapped == -1 or r_start_ungapped == -1:
            new_column_data.append("Error: Primers not found in reference sequence")
            continue
            
        r_end_ungapped = r_start_ungapped + len(r_template)

        try:
            start_gapped = ungapped_map[f_start_ungapped]
            end_gapped = ungapped_map[r_end_ungapped - 1] + 1
        except KeyError:
            new_column_data.append("Error: Could not map coordinates")
            continue
        # --- Koordinačių radimo pabaiga ---

        genome_specific_amplicons = []
        for record in block:
            genome_name = extract_genome_name(record.id)
            gapped_amplicon = str(record.seq[start_gapped:end_gapped])
            clean_amplicon = gapped_amplicon.replace('-', '')
            genome_specific_amplicons.append(f"{genome_name}:{clean_amplicon}")

        new_column_data.append(";".join(genome_specific_amplicons))

    primers_df['genome_specific_amplicons'] = new_column_data
    primers_df.to_csv(args.output, index=False)
    
    print(f"\nAtaskaita sėkmingai sukurta.")
    print(f"Rezultatai su visais variantais išsaugoti: {args.output}")

if __name__ == "__main__":
    main()