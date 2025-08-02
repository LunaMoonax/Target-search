import sys
import argparse
import pandas as pd
from Bio import SeqIO
from itertools import combinations
from collections import defaultdict

def extract_genome_name(seq_id):
    """
    Ištraukia genomo pavadinimą iš sekos ID.
    Pvz.: "Rickettsia_rhipicephali;0_0_236" → "Rickettsia_rhipicephali"
    """
    return seq_id.split(';')[0] if ';' in seq_id else seq_id

def parse_alignment_blocks_smart(fasta_file):
    """
    Skaito multi-FASTA failą ir suskirsto jį į blokus dinamiškai.
    Naujas blokas prasideda, kai pasikartoja pirmas genomas (Rickettsia_rhipicephali).
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        print("KLAIDA: Failas tuščias!")
        return
    
    print(f"Nuskaityta {len(records)} sekų iš {fasta_file}")
    
    # Randame pirmą genomą
    first_genome = extract_genome_name(records[0].id)
    print(f"Pirmas genomas (blokų žymeklis): {first_genome}")
    
    # Dinamiškai dalijame į blokus pagal pirmo genomo pasikartojimą
    blocks = []
    current_block = []
    
    for i, rec in enumerate(records):
        genome_name = extract_genome_name(rec.id)
        
        # Jei tai pirmas genomas ir jau turime sekų bloke, pradedame naują bloką
        if genome_name == first_genome and current_block:
            blocks.append(current_block)
            current_block = [rec]
        else:
            current_block.append(rec)
    
    # Pridedame paskutinį bloką
    if current_block:
        blocks.append(current_block)
    
    print(f"Rasta {len(blocks)} blokų")
    
    # Analizuojame blokų struktūrą
    block_sizes = [len(block) for block in blocks]
    unique_sizes = set(block_sizes)
    
    print(f"Blokų dydžiai: {dict(zip(range(1, len(blocks)+1), block_sizes))}")
    print(f"Unikalūs blokų dydžiai: {sorted(unique_sizes)}")
    
    if len(unique_sizes) > 1:
        print(f"ĮSPĖJIMAS: Skirtingi blokų dydžiai! Tai normalu, jei skirtingose dalyse trūksta kai kurių genomų.")
    
    blocks_created = 0
    for block in blocks:
        if not block:
            continue
            
        blocks_created += 1
        
        # Parodome bloko sudėtį kas 20 blokų arba jei blogai atrodo
        if blocks_created % 20 == 0 or len(block) != block_sizes[0]:
            block_genomes = [extract_genome_name(rec.id) for rec in block]
            print(f"  Blokas #{blocks_created}: {len(block)} genomų - {block_genomes[:3]}...")
        
        yield block
    
    print(f"Sukurta {blocks_created} tinkamų blokų")

def calculate_avg_pairwise_identity(sequences):
    """
    Apskaičiuoja vidutinį porinį identiškumą duotam sekų sąrašui.
    Grąžina reikšmę nuo 0.0 (visiškai skirtingi) iki 1.0 (identiški).
    """
    if len(sequences) < 2:
        return 1.0

    pairs = combinations(sequences, 2)
    total_pairs = 0
    total_identity = 0

    for seq1, seq2 in pairs:
        # Ignoruojame pozicijas su tarpais ar N
        valid_positions = 0
        matches = 0
        
        for a, b in zip(seq1, seq2):
            if a not in '-N' and b not in '-N':
                valid_positions += 1
                if a == b:
                    matches += 1
        
        if valid_positions > 0:
            identity = matches / valid_positions
            total_identity += identity
            total_pairs += 1

    return total_identity / total_pairs if total_pairs > 0 else 0.0

def has_valid_nucleotides(sequence):
    """
    Tikrina ar seka turi tik validžius nukleotidus (A, T, G, C).
    """
    valid_chars = set('ATGC')
    return all(c in valid_chars for c in sequence.upper())

def calculate_uniqueness_score(variants):
    """
    Apskaičiuoja, kiek unikalių variantų yra (idealiai - kiekvienas genomas turi unikalų).
    Grąžina santykį: unique_variants / total_genomes
    """
    total_genomes = sum(len(genomes) for genomes in variants.values())
    unique_variants = len(variants)
    return unique_variants / total_genomes if total_genomes > 0 else 0

def find_species_specific_targets(fasta_file, target_size=18, min_block_len=200, 
                                min_similarity=0.65, max_similarity=0.75, min_uniqueness=0.8):
    """
    Skenuoja sulyginimo blokus ir ieško 18 nt taikinių, kurie:
    1. Turi vidutinį panašumą ~70% (pakankamai panašūs, bet ne identiški)
    2. Turi kuo daugiau unikalių variantų (idealiai - kiekvienas genomas skirtingas)
    """
    candidates = []
    block_id = 0
    total_windows_checked = 0
    
    print(f"Ieškoma rūšims specifinių {target_size} nt taikinių...")
    print(f"Panašumas: {min_similarity*100:.1f}% - {max_similarity*100:.1f}%")
    print(f"Minimalus unikalumas: {min_uniqueness*100:.1f}%")
    print(f"Minimalus bloko ilgis: {min_block_len} nt")

    for block in parse_alignment_blocks_smart(fasta_file):
        block_id += 1
        
        # 1. Patikriname bloko ilgį
        if not block:  # Tuščias blokas
            continue
            
        block_len = len(block[0].seq)
        if block_len < min_block_len:
            print(f"Blokas #{block_id} per trumpas ({block_len} < {min_block_len}), praleidžiama")
            continue

        if block_id % 20 == 0:  # Dažnesni pranešimai
            print(f"  Apdorojamas blokas #{block_id}, ilgis: {block_len} nt")

        # 2. Slenkame per bloką target_size nt langu
        block_candidates = 0
        for pos in range(block_len - target_size + 1):
            total_windows_checked += 1
            
            target_seqs_map = {}
            skip_window = False
            
            for rec in block:
                seq_segment = str(rec.seq[pos : pos + target_size]).upper()
                
                # Praleidžiame langus su tarpais ('-') arba 'N'
                if not has_valid_nucleotides(seq_segment):
                    skip_window = True
                    break
                
                # Išsaugome su genomo pavadinimu (be koordinačių)
                genome_name = extract_genome_name(rec.id)
                target_seqs_map[genome_name] = seq_segment
            
            if skip_window:
                continue

            # 3. Apskaičiuojame panašumą
            sequences = list(target_seqs_map.values())
            avg_identity = calculate_avg_pairwise_identity(sequences)

            # 4. Tikriname, ar panašumas atitinka kriterijų
            if min_similarity <= avg_identity <= max_similarity:
                # Grupuojame genomus pagal rastą sekos variantą
                variants = defaultdict(list)
                for genome, seq in target_seqs_map.items():
                    variants[seq].append(genome)

                # 5. Apskaičiuojame unikalumo balą
                uniqueness_score = calculate_uniqueness_score(variants)
                
                # 6. Tikriname unikalumo kriterijų
                if uniqueness_score >= min_uniqueness:
                    candidate_info = {
                        "block_id": block_id,
                        "start_pos": pos,
                        "end_pos": pos + target_size,
                        "avg_pairwise_identity": round(avg_identity, 3),
                        "num_variants": len(variants),
                        "uniqueness_score": round(uniqueness_score, 3),
                        "variants": variants,
                        "block_coordinates": f"{block[0].id.split(';')[1] if ';' in block[0].id else 'unknown'}"
                    }
                    candidates.append(candidate_info)
                    block_candidates += 1

        if block_candidates > 0:
            print(f"    Bloke #{block_id} rasta {block_candidates} kandidatų")

    print(f"Analizė baigta. Patikrinta {total_windows_checked:,} langų.")
    print(f"Rasta {len(candidates)} kandidatų iš {block_id} blokų.")
    return candidates

def main():
    # --- PRADŽIA: Rankinis parametrų nustatymas ---
    
    TARGET_SIMILARITY = 0.70  # Vidutinis panašumas ~70%
    SIMILARITY_WINDOW = 0.05  # ±5%
    MIN_UNIQUENESS = 0.7      # Bent 70% variantų turi būti unikalūs
    MIN_BLOCK_LENGTH = 100    # Minimalus bloko ilgis

    # --- PABAIGA: Rankinis parametrų nustatymas ---

    parser = argparse.ArgumentParser(description="Randa rūšims specifinius 18 nt taikinius blokiniame alignment faile.")
    parser.add_argument("input_fasta", help="Blokinis alignment FASTA failas.")
    parser.add_argument("output_csv", help="CSV failas rezultatams išsaugoti.")
    args = parser.parse_args()

    min_sim = TARGET_SIMILARITY - SIMILARITY_WINDOW
    max_sim = TARGET_SIMILARITY + SIMILARITY_WINDOW

    print("=== RŪŠIMS SPECIFINIŲ TAIKINIŲ PAIEŠKA ===")
    print(f"Įvesties failas: {args.input_fasta}")
    print(f"Panašumo riba: {min_sim*100:.1f}% - {max_sim*100:.1f}%")
    print(f"Min. unikalumas: {MIN_UNIQUENESS*100:.1f}%")
    print(f"Min. bloko ilgis: {MIN_BLOCK_LENGTH} nt")
    print("==========================================\n")

    try:
        candidates = find_species_specific_targets(
            args.input_fasta, 
            min_block_len=MIN_BLOCK_LENGTH,
            min_similarity=min_sim, 
            max_similarity=max_sim,
            min_uniqueness=MIN_UNIQUENESS
        )
    except Exception as e:
        print(f"KLAIDA: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

    if not candidates:
        print(f"Tinkamų taikinių nerasta su nustatytais kriterijais.")
        print("Pabandykite:")
        print("  - Sumažinti MIN_UNIQUENESS (pvz., 0.5)")
        print("  - Pakeisti TARGET_SIMILARITY reikšmę")
        print("  - Padidinti SIMILARITY_WINDOW")
        print("  - Sumažinti MIN_BLOCK_LENGTH")
        
        # Sukuriame tuščią failą
        pd.DataFrame(columns=["block_id", "start_pos", "end_pos", "avg_pairwise_identity", 
                             "num_variants", "uniqueness_score", "sequence_variant", 
                             "genome_count", "genomes", "block_coordinates"]).to_csv(args.output_csv, index=False)
        return

    # Formuojame galutinį CSV failą
    output_data = []
    for cand in candidates:
        for seq_variant, genomes in cand['variants'].items():
            output_data.append({
                "block_id": cand['block_id'],
                "start_pos": cand['start_pos'],
                "end_pos": cand['end_pos'],
                "avg_pairwise_identity": cand['avg_pairwise_identity'],
                "num_variants": cand['num_variants'],
                "uniqueness_score": cand['uniqueness_score'],
                "sequence_variant": seq_variant,
                "genome_count": len(genomes),
                "genomes": ";".join(sorted(genomes)),
                "block_coordinates": cand['block_coordinates']
            })

    df = pd.DataFrame(output_data)
    
    # Rūšiuojame pagal unikalumo balą (didėjimo tvarka) ir panašumą
    df = df.sort_values(['uniqueness_score', 'avg_pairwise_identity', 'block_id', 'start_pos'], ascending=[False, True, True, True])
    
    df.to_csv(args.output_csv, index=False)
    
    print(f"\n=== REZULTATAI ===")
    print(f"Rasta {len(candidates)} unikalių pozicijų")
    print(f"Iš viso variantų: {len(df)}")
    print(f"Blokų spektras: {df['block_id'].min()} - {df['block_id'].max()}")
    print(f"Panašumo spektras: {df['avg_pairwise_identity'].min():.3f} - {df['avg_pairwise_identity'].max():.3f}")
    print(f"Unikalumo spektras: {df['uniqueness_score'].min():.3f} - {df['uniqueness_score'].max():.3f}")
    print(f"Variantų per taikinį: {df['num_variants'].value_counts().to_dict()}")
    
    # Papildoma statistika
    best_candidates = df[df['uniqueness_score'] == df['uniqueness_score'].max()]
    print(f"Geriausi kandidatai (unikalumas = {df['uniqueness_score'].max():.3f}): {len(best_candidates)}")
    
    print(f"Rezultatai išsaugoti: {args.output_csv}")

if __name__ == "__main__":
    main()