import pandas as pd
import re

def main():
    print("Running...")
    input_fasta = "mobidb_search_2026-05-05T08-28-46.fasta"
    score_patterns = ["derived-binding_mode_disorder_to_order-mobi"]
    data = get_mobidb_scores(input_fasta, score_patterns)
    cols_to_check = ["derived_binding_mode_disorder_to_order_mobi"]
    filtered_data = data.dropna(subset=cols_to_check, thresh=1)
    print(filtered_data.tail())



def get_mobidb_scores(input_fasta: str, score_patterns: list[str]) -> pd.DataFrame:
    '''
    Produces an array containing the protein sequence along with certain score values representing the region of interest, by procesing a mobidb FASTA input file containing protein sequences along with score values.

    Args:
        input_fasta (str): Path to the input mobidb FASTA file containg the score values
        score_patterns (list[str]): A list containing the patterns to identify the score headers.
    
    Returns:
        (pd.DataFrame): A pandas DateFrame containing columns containing uniprot_id, sequence, and the given patterns.

    '''
    results = []
    
    # open the input file
    with open(input_fasta, "r") as file:
        # set current entry as None
        current_entry = None

        # loop through each line of the file
        for line in file:
            line = line.strip()

            if "sequence" in line:

                # if there is already a entry, save it to results before starting new
                if current_entry:
                    results.append(current_entry)

                
                # get the uniprot id for the current entry
                current_id = line.split("|")[0].strip(">").strip()

                # initialize the current entry with the current uniprot_id, sequence, and with None for the score patterns
                current_entry = {
                    "uniprot_id": current_id,
                    "sequence": "",
                    **{p: None for p in score_patterns}
                }                 
            
            # get the protein sequence for the current entry
            elif line[:1].isalpha() and current_entry:
                current_entry["sequence"] += line

            # get the scores for the given patterns and append it to the list
            elif current_entry:
                for pattern in score_patterns:
                    pattern_match = current_entry["uniprot_id"] + "|" + pattern
                    if pattern_match in line:
                        next_line = next(file, "")
                        current_entry[pattern] = next_line.strip()
        
        # Append the last entry
        if current_entry:
            results.append(current_entry)
        

    df = pd.DataFrame(results)
    df.columns = df.columns.str.replace('-', '_')
    return df



if __name__ == "__main__":
	main()