import sys

def xmfa_to_fasta(input_file, output_file):
    """Convert XMFA to FASTA format - just format conversion, no filtering"""
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        current_header = None
        current_sequence = ""
        sequence_count = 0
        
        for line in infile:
            line = line.strip()
            
            # Skip empty lines and block separators
            if not line or line.startswith('='):
                continue
            
            # Header line
            if line.startswith('>'):
                # Write previous sequence if it exists
                if current_header and current_sequence:
                    outfile.write(f"{current_header}\n")
                    # Write sequence in 80-character lines
                    for i in range(0, len(current_sequence), 80):
                        outfile.write(current_sequence[i:i+80] + "\n")
                    sequence_count += 1
                
                # Start new sequence
                current_header = line
                current_sequence = ""
            
            # Sequence line
            else:
                current_sequence += line
        
        # Don't forget the last sequence
        if current_header and current_sequence:
            outfile.write(f"{current_header}\n")
            for i in range(0, len(current_sequence), 80):
                outfile.write(current_sequence[i:i+80] + "\n")
            sequence_count += 1
    
    print(f"Converted {sequence_count} sequences from {input_file} to {output_file}")

def main():
    if len(sys.argv) != 3:
        print("Usage: python xmfa_to_fasta.py input.xmfa output.fasta")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        xmfa_to_fasta(input_file, output_file)
        print("Conversion completed successfully!")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
