import pandas as pd	

def main():
	print("Test!")
	data = {
    "uniprot_id": ["Q04917", "P62258", "P31947", "P27348", "P63104"],
    "binding_region_seq": ["QQDEEAGEGN", "QGDGEEQNKEALQDVEDENQ", "NAGEEGGEAPQEPQS", "SAGEECDAAEGAEN", "QGDEAEAGEGGEN"],
    "start_pos": [236, 235, 233, 231, 232],
    "end_pos": [246, 255, 248, 245, 245],
    "label": ["disorder_to_disorder"] * 5,
    "id_col": ["Q04917_236_246", "P62258_235_255", "P31947_233_248", "P27348_231_245", "P63104_232_245"]
    }
	df = pd.DataFrame(data)
	dataframe_to_fasta(df, "id_col", "binding_region_seq", "binding_regions_test.fasta")

	

def dataframe_to_fasta(df: pd.DataFrame, id_col: str, seq_col: str, output_file: str) -> None:
	'''
	Convert a DataFrame to FASTA format and save it to a file.
	
	Parameters:
		df (pd.DataFrame): The DataFrame to convert.
		id_col (str): The column name for the sequence IDs.
		seq_col (str): The column name for the sequence data.
		output_file (str): The path to the output FASTA file.
		
	'''
	with open(output_file, 'w') as f:
		for _, row in df.iterrows():
			f.write(f">{row[id_col]}\n{row[seq_col]}\n")

def process_mmseqs_clusters(input_clusterdb, output_clusters):
	
    cluster_num = 0
    with open(input_clusterdb, "r") as file:
        for line in file:
            rep_id, member_id = line.strip().split("\t")
			if rep_id == member_id:
				
            
			

                


if __name__ == "__main__":
	main()