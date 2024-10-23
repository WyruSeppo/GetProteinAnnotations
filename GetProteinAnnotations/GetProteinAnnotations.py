from Bio import SeqIO
import requests
import pandas as pd

def read_tsv(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        return df
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return None

def get_uniprot_annotation(protein_id):
    url = f"https://www.uniprot.org/uniprotkb/{protein_id}.txt"
    response = requests.get(url)
    if response.status_code == 200:
        text_data = response.text
        protein_name = None
        function = None
        
        for line in text_data.splitlines():
            if line.startswith("ID"):
                protein_name = line.split()[1]
            elif line.startswith("CC   -!- FUNCTION:"):
                function = line[len("CC   -!- FUNCTION:"):].strip()
        
        return protein_name, function
    else:
        return None, None

def annotate_tsv(tsv_file):
    df = read_tsv(tsv_file)
    annotated_data = []

    if df is not None:
        for protein_id in df.iloc[:, 1]:
            protein_name, function = get_uniprot_annotation(protein_id)
            annotated_data.append({
                'Protein_ID': protein_id,
                'Protein_Name': protein_name,
                'Function': function
            })
            print(function)
    return pd.DataFrame(annotated_data)


#test it
tsv_file_path = "spore_formers_proteinIds"  # Update with your file path
annotated_df = annotate_tsv(tsv_file_path)

print(annotated_df)
