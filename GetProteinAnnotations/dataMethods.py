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
    try:
        response = requests.get(url)
    except Exception as e:
        print(f"An error occurred while calling uniprot API: {protein_id}")
        return None, None, None

    if response.status_code == 200:
        text_data = response.text
        protein_name = None
        function = None
        pfam = None
        
        for line in text_data.splitlines():
            if line.startswith("ID"):
                protein_name = line.split()[1]
            elif line.startswith("CC   -!- FUNCTION:"):
                function = line[len("CC   -!- FUNCTION:"):].strip()
            elif line.startswith("DR   Pfam;"):
                pfam = line.split()[2].removesuffix(';')

        return protein_name, function, pfam
    else:
        return None, None, None

def get_pfam_annotation(protein_id):
    api_url = "https://www.ebi.ac.uk/interpro/api"
    url = f"{api_url}/entry/pfam/{protein_id}"
    url += "?page_size=200"

    try:
        response = requests.get(url)
    except Exception as e:
        print(f"An error occurred while calling pfam API: {protein_id}")
        return None

    if response.status_code == 200:
        description_text = None
        jsonData = response.json()

        if jsonData['metadata']['description'] and isinstance(jsonData['metadata']['description'], list) and len(jsonData['metadata']['description']) > 0:
            description_text = jsonData['metadata']['description'][0]['text']
            description_text = description_text.removesuffix('</p>')
            description_text = description_text.removeprefix('<p>')
            
        return description_text
    else:
        return None


def annotate_tsv(tsv_file):
    df = read_tsv(tsv_file)
    annotated_data = []
    counter = 0
    total = len(df)

    if df is not None:
        for protein_id in df.iloc[:, 1]:
            print("\r" + str(counter // total) + "% " + str(counter) + "/" + str(total));
            counter += 1
            protein_name, function, pfamID = get_uniprot_annotation(protein_id)
            pfamAnnotation = get_pfam_annotation(pfamID)

            annotated_data.append({
                'Protein_ID': protein_id,
                'Pfam_ID':pfamID,
                'Protein_Name': protein_name,
                'Unitprot_Function': function,
                'Pfam_Description': pfamAnnotation
            })
            #print((protein_id or 'NA') + " " + (pfamID or 'NA') + " " + (protein_name or 'NA') + " " + (function or 'NA') + " " + (pfamAnnotation or 'NA'))
    return pd.DataFrame(annotated_data)
