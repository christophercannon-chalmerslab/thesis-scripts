import os

def concatenate_fasta_files(directory, output_file):
    fasta_files = [f for f in os.listdir(directory) if f.endswith('.fasta')]
    
    if not fasta_files:
        print("No FASTA files found in the directory.")
        return

    with open(output_file, 'w') as outfile:
        for fasta_file in fasta_files:
            with open(os.path.join(directory, fasta_file), 'r') as infile:
                outfile.write(infile.read())

    print(f"Concatenation complete! Output saved as: {output_file}")

def main():
    directory = input("Please enter the directory path containing the fasta files: ").strip()
    
    if not os.path.isdir(directory):
        print("Invalid directory path.")
        return
    
    output_file = input("Please enter the name for the output fasta file (e.g., output.fasta): ").strip()
    
    if not output_file.endswith('.fasta'):
        output_file += '.fasta'

    concatenate_fasta_files(directory, output_file)

if __name__ == "__main__":
    main()
