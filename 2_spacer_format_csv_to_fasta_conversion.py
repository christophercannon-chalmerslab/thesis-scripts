import pandas as pd

def csv_to_fasta():
    csv_file = input("Enter the name of the CSV file (without suffix): ") + '.csv'
    fasta_file = input("Enter the name for the output FASTA file (without suffix): ") + '.fasta'
    
    df = pd.read_csv(csv_file)
    
    df.columns = df.columns.str.strip().str.lower()
    
    print("Columns in CSV file:", df.columns)
    
    base_name = csv_file.rsplit('.', 1)[0]
    
    with open(fasta_file, 'w') as fasta_out:
        for index, row in df.iterrows():
            sequence = row['spacer'] 
            count = row['count'] 
            
            header = f">{base_name}_{count}"
            
            fasta_out.write(f"{header}\n{sequence}\n")

csv_to_fasta()
