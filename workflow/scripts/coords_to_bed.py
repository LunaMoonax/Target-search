input_file = snakemake.input[0]
output_file = snakemake.output[0]

with open(input_file) as fin, open(output_file, "w") as fout:
    for i, line in enumerate(fin):
        if i < 5:
            continue  # Skip header
        parts = line.strip().split()
        if len(parts) < 13:
            continue  # Skip malformed lines
        try:
            s1 = int(parts[0])
            e1 = int(parts[1])
            chrom = parts[12]  # Reference sequence ID (Borrelia block)
            start = min(s1, e1) - 1  # BED is 0-based
            end = max(s1, e1)
            fout.write(f"{chrom}\t{start}\t{end}\n")
        except ValueError:
            continue  # Skip lines with non-integer coordinates
