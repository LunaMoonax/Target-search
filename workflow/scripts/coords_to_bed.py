input_file = snakemake.input[0]
output_file = snakemake.output[0]

regions = []

with open(input_file) as fin:
    for i, line in enumerate(fin):
        if i < 5:
            continue  # Skip header
        parts = line.strip().split()
        if len(parts) < 13:
            continue  # Skip malformed lines
        try:
            s1 = int(parts[0])
            e1 = int(parts[1])
            chrom = parts[11]  # Reference sequence ID
            start = min(s1, e1) - 1  # BED is 0-based
            end = max(s1, e1)
            regions.append((chrom, start, end))
        except ValueError:
            continue  # Skip bad lines

# Sort by chromosome and then by start coordinate
regions.sort(key=lambda x: (x[0], x[1]))

# Write to output
with open(output_file, "w") as fout:
    for chrom, start, end in regions:
        fout.write(f"{chrom}\t{start}\t{end}\n")

