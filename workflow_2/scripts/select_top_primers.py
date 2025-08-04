#!/usr/bin/env python3
"""
Selects the best primer pair for each unique target location based on the lowest penalty score.
"""
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Selects the best primer pair for each unique target.")
    parser.add_argument("input_csv", help="CSV file with validated primer candidates.")
    parser.add_argument("output_csv", help="Output CSV for the top primer for each target.")
    args = parser.parse_args()

    try:
        df = pd.read_csv(args.input_csv)
        if df.empty:
            print("Įvesties failas tuščias, nėra ką atrinkti.")
            df.to_csv(args.output_csv, index=False)
            return
    except FileNotFoundError:
        print(f"KLAIDA: Nerastas įvesties failas: {args.input_csv}")
        return

    print(f"Atrankomi geriausi kandidatai iš {len(df)} eilučių...")

    # Paliekame tik tas eilutes, kurios praėjo validaciją
    valid_df = df[df['is_valid'] == True].copy()

    if valid_df.empty:
        print("Nerasta nei vienos validžios pradmenų poros.")
        valid_df.to_csv(args.output_csv, index=False)
        return

    # Rūšiuojame pagal baudos balą (nuo mažiausio)
    valid_df.sort_values(by='total_penalty', inplace=True)

    # Kiekvienai unikaliai taikinio pozicijai ('block_id', 'start_pos') paliekame tik pirmą (t.y. geriausią) eilutę
    top_primers_df = valid_df.drop_duplicates(subset=['block_id', 'start_pos'], keep='first')

    top_primers_df.to_csv(args.output_csv, index=False)

    print(f"\nAtrinkta {len(top_primers_df)} geriausių pradmenų porų (po vieną kiekvienam taikiniui).")
    print(f"Galutinis rezultatas išsaugotas: {args.output_csv}")

if __name__ == "__main__":
    main()