import subprocess
import os
import csv

def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.
    
    Parameters:
    - seq: The DNA sequence to reverse complement.
    
    Returns:
    - The reverse complement of the sequence.
    """
    complement = str.maketrans('ACGTacgt', 'TGCAtgca')
    return seq.translate(complement)[::-1]

def extract_pam(sequence, sstart, send, sstrand):
    """
    Extract the PAM for a sequence alignment.
    
    For the plus strand, extract the three bases before sstart.
    For the minus strand, extract the three bases after sstart and return the reverse complement.
    
    Parameters:
    - sequence: The full sequence from the database (subject).
    - sstart: Start position of the alignment (1-based).
    - send: End position of the alignment (1-based).
    - sstrand: Strand direction of the alignment ('plus' or 'minus').
    
    Returns:
    - The PAM sequence (three bases).
    """
    if sstrand == 'plus':
        pam_position = sstart - 4 
        if pam_position >= 0:
            return sequence[pam_position:sstart-1]
    elif sstrand == 'minus':
        pam_position = sstart
        if pam_position + 2 < len(sequence):
            pam_sequence = sequence[pam_position:pam_position + 3]
            return reverse_complement(pam_sequence)
    return "N/A"

def read_sequences(file_path):
    """
    Read sequences from a FASTA file and return a dictionary mapping sequence IDs to sequences.
    
    Parameters:
    - file_path: Path to the FASTA file.
    
    Returns:
    - A dictionary mapping sequence IDs to sequences.
    """
    sequences = {}
    print(f"Reading sequences from {file_path}...")
    
    try:
        with open(file_path, 'r') as f:
            seq_id = None
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    seq_id = line[1:]
                    sequences[seq_id] = ''
                elif seq_id:
                    sequences[seq_id] += line
        print(f"Loaded {len(sequences)} sequences.")
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
    return sequences

def perform_blast(query_fasta, db_name):
    """
    Perform a BLAST search using a query FASTA file and a BLAST database.
    
    Parameters:
    - query_fasta: Path to the query FASTA file.
    - db_name: The base name of the BLAST database.
    
    Returns:
    - A temporary file with the BLAST results.
    """
    temp_output_file = "temp_results.txt"
    
    blastn_command = [
        "blastn",
        "-query", query_fasta,
        "-db", db_name,
        "-out", temp_output_file, 
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand"
    ]
    
    print(f"Running BLAST command: {' '.join(blastn_command)}")
    
    try:
        subprocess.run(blastn_command, check=True)
        print("BLAST search completed.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred during BLAST search: {e}")
        return None
    
    return temp_output_file

def process_blast_results(temp_file, query_sequences, db_sequences):
    """
    Process the BLAST results, extracting relevant information and generating the PAM sequences.
    
    Parameters:
    - temp_file: Path to the temporary file containing BLAST results.
    - query_sequences: A dictionary of query sequences.
    - db_sequences: A dictionary of database sequences.
    
    Returns:
    - Two lists: one for p1656 results and one for JB028 results.
    - A set of sequences that have already aligned to p1656.
    """
    jb028_results = []
    p1656_results = []
    aligned_to_p1656 = set()

    with open(temp_file, 'r') as temp_file:
        lines = temp_file.readlines()
    
    for line in lines:
        fields = line.strip().split('\t')
        if len(fields) >= 13:
            try:
                pident = float(fields[2])
                length = int(fields[3])
                if pident == 100.0 and length == 33:
                    qseqid = fields[0]
                    sseqid = fields[1]
                    sstart = int(fields[8])  
                    send = int(fields[9]) 
                    sstrand = fields[12]  
                    
                    query_sequence = query_sequences.get(qseqid, 'Unknown')
                    
                    db_sequence = db_sequences.get(sseqid, '')
                    pam_sequence = extract_pam(db_sequence, sstart, send, sstrand)
                    
                    if query_sequence:
                        pam_sequence += query_sequence[0]
                    pam_sequence = pam_sequence.upper() 

                    if "JB028" in sseqid:
                        jb028_results.append([qseqid, query_sequence, sstart] + fields[1:] + [pam_sequence])
                    elif "p1656" in sseqid:
                        p1656_results.append([qseqid, query_sequence, sstart] + fields[1:] + [pam_sequence])
                        aligned_to_p1656.add(qseqid)
            except ValueError:
                pass
    
    return jb028_results, p1656_results, aligned_to_p1656

def write_output(jb028_results, p1656_results, aligned_to_p1656, output_file):
    """
    Write the BLAST results to a CSV file in the required format.
    
    Parameters:
    - jb028_results: Results from the JB028 database.
    - p1656_results: Results from the p1656 database.
    - aligned_to_p1656: Set of sequences that have already aligned to p1656.
    - output_file: Path to the output CSV file.
    """
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        column_headings = ["Count", "Spacer ID", "Sequence", "Database", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sstrand", "PAM"]
        csvwriter.writerow(column_headings)

        count = 1
        for row in p1656_results:
            csvwriter.writerow([count, row[0], row[1]] + row[3:])
            count += 1
        
        for row in jb028_results:
            if row[0] not in aligned_to_p1656:
                csvwriter.writerow([count, row[0], row[1]] + row[3:])
                count += 1
    
    print(f"Results saved to {output_file}.")

def main():
    print("Script started.")
    
    query_fasta = input("Enter the filename of the fasta file with small sequences: ")
    db_name = input("Enter the base name of the BLAST database (e.g., JB028-p1656): ")

    if not os.path.isfile(query_fasta):
        print(f"Error: The file {query_fasta} does not exist.")
        return
    db_fasta_file = f"{db_name}.fasta"
    if not os.path.isfile(db_fasta_file):
        print(f"Error: The file {db_fasta_file} does not exist.")
        return
    
    query_sequences = read_sequences(query_fasta)
    db_sequences = read_sequences(db_fasta_file)
    
    temp_results = perform_blast(query_fasta, db_name)
    if temp_results is None:
        return
    
    jb028_results, p1656_results, aligned_to_p1656 = process_blast_results(temp_results, query_sequences, db_sequences)
    
    output_file = f"{os.path.splitext(query_fasta)[0]}-blast-{db_name}.csv"
    
    write_output(jb028_results, p1656_results, aligned_to_p1656, output_file)
    
    os.remove(temp_results)
    print("Script completed successfully.")

if __name__ == "__main__":
    main()
