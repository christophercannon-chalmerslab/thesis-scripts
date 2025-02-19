from Bio import SeqIO

def remove_duplicates_with_count_sorted(input_fasta, output_fasta):
    sequence_counts = {} 
    sequence_records = {} 

    for record in SeqIO.parse(input_fasta, "fasta"):
        seq_str = str(record.seq)  
        if seq_str in sequence_counts:
            sequence_counts[seq_str] += 1
        else:
            sequence_counts[seq_str] = 1
            sequence_records[seq_str] = record 

    sorted_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)

    records_to_write = []
    for seq_str, count in sorted_sequences:
        record = sequence_records[seq_str] 
        record.id += f"_count{count}" 
        record.description = "" 
        records_to_write.append(record)

    SeqIO.write(records_to_write, output_fasta, "fasta")
    print(f"Processed sequences. Unique sequences sorted by count written to {output_fasta}")

input_fasta = input("Enter the input FASTA file name: ")
output_fasta = input("Enter the output FASTA file name: ")

remove_duplicates_with_count_sorted(input_fasta, output_fasta)
