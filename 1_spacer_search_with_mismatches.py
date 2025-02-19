import os
import csv

def reverse_complement(sequence):
    complement = str.maketrans("acgt", "tgca")
    return sequence.translate(complement)[::-1]

def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance between two sequences."""
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

def allowed_mismatches(sequence, query, max_mismatches):
    """Find all positions where the query matches the sequence with allowed mismatches."""
    matches = []
    query_len = len(query)
    
    for i in range(len(sequence) - query_len + 1):
        subseq = sequence[i:i + query_len]
        mismatches = hamming_distance(subseq, query)
        if mismatches <= max_mismatches:
            matches.append((i, subseq, mismatches))
    return matches

def process_fastq(input_fastq, output_csv, max_mismatches):
    input_fastq = input_fastq if input_fastq.endswith(".fastq") else input_fastq + ".fastq"
    output_csv = output_csv if output_csv.endswith(".csv") else output_csv + ".csv"

    if not os.path.exists(input_fastq):
        print(f"Error: The input file '{input_fastq}' does not exist.")
        return

    search_sequence = "tatgtttagagtgttccccgcgccagcggggataaacc".lower()

    results_found = False
    all_results = [] 

    try:
        with open(input_fastq, 'r') as fastq_file:
            while True:
                header = fastq_file.readline().strip() 
                sequence = fastq_file.readline().strip().lower()
                fastq_file.readline()
                fastq_file.readline() 

                if not sequence:
                    break

                forward_matches = allowed_mismatches(sequence, search_sequence, max_mismatches)
                for position, mismatched_seq, mismatches in forward_matches:
                    following_33bp = sequence[position + len(search_sequence):position + len(search_sequence) + 33]
                    concatenated = mismatched_seq + following_33bp
                    all_results.append([concatenated, mismatches, "Forward", following_33bp, header])
                    results_found = True

                rev_comp_sequence = reverse_complement(sequence)
                reverse_matches = allowed_mismatches(rev_comp_sequence, search_sequence, max_mismatches)
                for position, mismatched_seq, mismatches in reverse_matches:
                    following_33bp = rev_comp_sequence[position + len(search_sequence):position + len(search_sequence) + 33]
                    concatenated = mismatched_seq + following_33bp
                    all_results.append([concatenated, mismatches, "Reverse", following_33bp, header])
                    results_found = True

        if results_found:
            all_results.sort(key=lambda x: x[2]) 

            for count, result in enumerate(all_results, start=1):
                result.insert(0, count) 
            try:
                with open(output_csv, 'w', newline='') as csv_file:
                    writer = csv.writer(csv_file)
                    writer.writerow(["Count", "Leader-Repeat-Spacer Unit", "Number of Mismatches", "Fastq Strand", "Spacer", "Metadata (Fastq Header)"])
                    writer.writerows(all_results)

                print(f"Results successfully written to '{output_csv}'.")

            except IOError as e:
                print(f"Error: Unable to write to the file '{output_csv}'. {e}")
        else:
            print(f"No matches found for '{search_sequence}'.")

    except FileNotFoundError:
        print(f"Error: The file '{input_fastq}' was not found.")
    except IOError as e:
        print(f"Error: There was an issue with reading the file '{input_fastq}'. {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    try:
        input_fastq = input("Enter the input FASTQ file name (without extension): ").strip()
        output_csv = input("Enter the output CSV file name (without extension): ").strip()
        max_mismatches = int(input("Enter the number of allowed mismatches: ").strip())

        if not input_fastq or not output_csv or max_mismatches < 0:
            print("Error: All inputs must be provided, and the number of mismatches must be non-negative.")
            exit()

        process_fastq(input_fastq, output_csv, max_mismatches)

    except ValueError:
        print("Error: Invalid number of mismatches provided. It must be an integer.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
