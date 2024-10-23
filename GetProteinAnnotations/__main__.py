from dataMethods import *

tsv_file_path = "spore_formers_proteinIds"
annotated_df = annotate_tsv(tsv_file_path)
annotated_df.to_csv('output.txt', sep='\t', index=False)