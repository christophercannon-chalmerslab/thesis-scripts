from Bio import SeqIO
from Bio import pairwise2

input_fasta = input("Enter the input FASTA file name (including path if necessary): ").strip()
output_fasta = input("Enter the output FASTA file name (including path if necessary): ").strip()

reference_seq = "gctttcgcagacgcgcggcgatacgctcacgca"
threshold_identity = 80  # Percent identity threshold

def calculate_percent_identity(seq1, seq2):
    """Calculate the percent identity between two sequences."""
    alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True, score_only=True)
    percent_identity = (alignment / max(len(seq1), len(seq2))) * 100
    return percent_identity

filtered_sequences = []
try:
    for record in SeqIO.parse(input_fasta, "fasta"):
        identity = calculate_percent_identity(str(record.seq), reference_seq)
        if identity < threshold_identity:
            filtered_sequences.append(record)

    SeqIO.write(filtered_sequences, output_fasta, "fasta")
    print(f"Filtered sequences have been saved to {output_fasta}")
except FileNotFoundError:
    print(f"Error: File {input_fasta} not found. Please check the file name and path.")
except Exception as e:
    print(f"An error occurred: {e}")
