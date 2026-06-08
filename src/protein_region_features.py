import math
import pandas as pd
import re
from collections import Counter

def main():
    print("Features")
    aa_1 = "TTAHLITGSAGSSLCSLNSSLHEGHSHSAGSPLFFNHVKIKLAD"
    aa_2 = "PGAEYKAKKAKGDVKKKGRPDPYAYIPLN"
    print(calculate_scd(aa_1), calculate_scd(aa_2))
    print(calculate_shd(aa_1), calculate_shd(aa_2))
    print(calculate_fcr(aa_1), calculate_fcr(aa_2))

def features_binding_region_df(input_df: pd.DataFrame, binding_mode: str, features: set[str]) -> pd.DataFrame:
    '''
    Extracts the specific derived binding region based on the score and return the corresponding features from that of the entire protein sequence.

    Args:
        input_df (pd.DataFrame): The input dataframe containing the uniprot_id, sequence, binding mode score, and the features for the entire protein sequence.
        binding_mode (str): The column name for the binding mode score to identify the binding region.
        features (set[str]): The set of column names for the features to be extracted for the binding region.
    
    Returns:
        pd.DataFrame: A dataframe containing the uniprot_id, binding region sequence, start and end position of the binding region, and the features for the binding region.


    '''
    # check if the required columns are present in the input dataframe
    input_df_col = set(input_df.columns)
    
    if binding_mode not in input_df_col:
        raise ValueError("The column for the binding mode score doesn't exist")
    if "sequence" not in input_df_col:
        raise ValueError("The sequence column doesn't exist")
    if "uniprot_id" not in input_df_col:
        raise ValueError("The uniprot_id column doesn't exist")
    if not features.issubset(input_df_col):
        raise ValueError("The column for the feature doesn't exist")
    
    extracted_regions = []
    # loop through each row of the input dataframe
    for row in input_df.itertuples():
        mask = getattr(row, binding_mode)
        mask = str(mask) if pd.notna(mask) else ""
        full_seq = getattr(row, "sequence")
        uniprot_id = getattr(row, "uniprot_id")
        # find the regions where the mask is 1 and extract the corresponding sequence and features for those regions
        for match in re.finditer(r"1+", mask):
            start = match.start()
            end = match.end()
            
            entry = {
                "uniprot_id": uniprot_id,
                "binding_region_seq": full_seq[start:end],
                "start_pos": start,
                "end_pos": end
            }
            # extract the features for the binding region and add it to the entry
            for f in features:
                feature_val = getattr(row, f)
                if isinstance(feature_val, str) or isinstance(feature_val, list):
                    entry[f] = feature_val[start:end]
                else:
                    entry[f] = None
            
            extracted_regions.append(entry)

    return pd.DataFrame(extracted_regions)

def get_average_feature_value(feature_values: list[float | str]) -> float | None:
    '''
    Calculates the average value of a feature for a given binding region. Considers only the non-None values for the calculation.
    
    Args:
        feature_values (list[float | str]): A list of feature values for a binding region.
    Returns:
        float | None: The average value of the feature for the binding region. Returns None if there are no non-None values.
    '''
    not_none_values = [float(val) for val in feature_values if val is not None]
    if len(not_none_values) == 0:
        return None
    else:
        return round(sum(not_none_values) / len(not_none_values), 4)
    
def get_sec_str_freq(sec_str: list[str], convert_3_state: bool = False) -> dict[str, float]:
    '''
    Calculates the frequency of each secondary structure type in a given binding region.
    
    Args:
        sec_str (list[str]): A list of secondary structure types for a binding region.
        covert_3_state (bool): A flag to indicate whether to convert the secondary structure types to 3-state (H, E, C) or not. If True, the secondary structure types will be converted to 3-state before calculating the frequency.
    Returns:
        dict[str, float]: A dictionary containing the frequency of each secondary structure type in the binding region.
    '''
    # Calculate the total count of secondary structure types
    total_count = len(sec_str)

    # If the total count is zero, return an empty dictionary to avoid division by zero
    if total_count == 0:
        return {}
    # Get the unique secondary structure types
    sec_str_types = set(sec_str)

    # If convert_3_state is True and all secondary structure types are in the 9-state set, convert them to 3-state (H, E, C)
    if convert_3_state and sec_str_types.issubset({"H", "G", "I", "E", "B", "T", "S", "P", "C"}):
        sec_str_types = {"H" if s in ["H", "G", "I"] else "E" if s in ["E", "B"] else "C" for s in sec_str_types}

    # Calculate the frequency of each secondary structure type and store it in a dictionary
    frequency_dict = {sec_str_type: round(sec_str.count(sec_str_type) / total_count, 4) for sec_str_type in sec_str_types}
    
    return frequency_dict

def get_amino_acid_composition(sequence: str) -> dict[str, float]:
    '''
    Calculates the normalized composition of each amino acid type in a given binding region.
    
    Args:
        sequence (str): The amino acid sequence for a binding region.
    Returns:
        dict[str, float]: A dictionary containing the normalized composition of each amino acid type in the binding region. The keys are the amino acid types and the values are the normalized frequencies of each amino acid type in the sequence.
    '''

    amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    # Calculate the total count of amino acids
    total_count = len(sequence)

    # If the total count is zero, return an empty dictionary to avoid division by zero
    if total_count == 0:
        return {}
    
    # Calculate the frequency of each amino acid type and store it in a dictionary
    counts = Counter(sequence.upper())

    # Create a dictionary for the normalized composition of each amino acid type
    composition_dict = {aa: round(counts.get(aa, 0) / total_count, 4) for aa in amino_acids}

    return composition_dict

def calculate_scd(sequence: str, need_rounding: bool = False, round_to: int = 4) -> float:
    '''
    Calculates the Sequence Charge Decoration (SCD) for a given amino acid sequence. SCD is a measure of the distribution of charged amino acids in a sequence and is calculated using the formula: SCD = (1 / N) * Sum_{m=2 to N} [ Sum_{n=1 to m-1} ( q_m * q_n * sqrt(m - n) ) ], where N is the length of the sequence, q_m and q_n are the charges of the amino acids at positions m and n, respectively, and (m - n) is the distance between the positions of the two amino acids in the sequence.

    Args:
        sequence (str): The amino acid sequence for a binding region.
        need_rounding (bool): A flag to indicate whether to round the result or not.
        round_to (int): The number of decimal places to round the result to.

    Returns:
        float: The calculated SCD value for the given sequence. 
    '''
    charge_map = {
        'R': 1,  # Arginine
        'K': 1,  # Lysine
        'D': -1, # Aspartic acid
        'E': -1  # Glutamic acid
    }

    # map the sequence to charges using the charge_map, non-charged amino acids will be mapped to 0
    seq = sequence.upper().strip()
    N = len(seq)

    if N == 0:
        return 0.0
    
    charges = [charge_map.get(aa, 0) for aa in seq]

    # double summation loop to calculate the SCD value
    scd_value = 0.0

    for m in range(1, N): # m goes from 1 to N-1 (0-based indexing)
        for n in range(m): # n goes from 0 to m-1 (0-based indexing)

            if charges[m] != 0 and charges[n] != 0:
                scd_value += charges[m] * charges[n] * math.sqrt(m - n)

    scd = scd_value / N
    if need_rounding:
        scd = round(scd, round_to)
    return scd

def calculate_shd(sequence: str, need_rounding: bool = False, round_to: int = 4) -> float:
    '''
    Calculates the Sequence Hydropathy Decoration (SHD) for a given amino acid sequence. SHD is a measure of the distribution of hydrophobic and hydrophilic amino acids in a sequence and is calculated using the formula: SHD = (1 / N) * Sum_{i=2 to N} [ Sum_{j=1 to i-1} ( (h_i + h_j) / (i - j) ) ], where N is the length of the sequence, h_m and h_n are the scaled hydropathy values of the amino acids at positions m and n, respectively, and (m - n) is the distance between the positions of the two amino acids in the sequence.
    
    Args:
        sequence (str): The amino acid sequence for a binding region.
        need_rounding (bool): A flag to indicate whether to round the result or not.
        round_to (int): The number of decimal places to round the result to.

    Returns:
        float: The calculated SHD value for the given sequence.
    
    '''
    # Normalized hydropathy scales used for the HPS-Urry 
    urry_hydropathy_map = {
    "W": 1.0,
    "Y": 0.897059,
    "F": 0.82353,
    "H": 0.764707,
    "P": 0.758824,
    "L": 0.720589,
    "I": 0.705883,
    "M": 0.676471,
    "V": 0.664707,
    "C": 0.64706,
    "A": 0.602942,
    "N": 0.588236,
    "S": 0.588236,
    "T": 0.588236,
    "G": 0.57353,
    "R": 0.558824,
    "K": 0.558824,
    "D": 0.382354,
    "E": 0.294119,
    "Q": 0.0
    }

    # map the sequence to hydropathy values using the urry_hydropathy_map, non-standard amino acids will be mapped to 0
    seq = sequence.upper().strip()
    N = len(seq)
    
    if N == 0:
        return 0.0
    
    hydropathy_values = [urry_hydropathy_map.get(aa, 0) for aa in seq]

    # double summation loop to calculate the SHD value
    shd_value = 0.0

    for i in range(1, N): # i goes from 1 to N-1 (0-based indexing)
        for j in range(i): # j goes from 0 to i-1 (0-based indexing)

            shd_value += (hydropathy_values[i] + hydropathy_values[j]) / (i - j)

    shd = shd_value / N
    if need_rounding:
        shd = round(shd, round_to)
    return shd

def calculate_fcr(sequence: str, need_rounding: bool = False, round_to: int = 4) -> float:
    '''
    Calculates the Fraction of Charged Residues (FCR) for a given amino acid sequence.

    Args:
        sequence (str): The amino acid sequence for a binding region.
        need_rounding (bool): A flag to indicate whether to round the result or not.
        round_to (int): The number of decimal places to round the result to.

    Returns:
        float: The calculated FCR value for the given sequence.
    
    '''
    
    charge_residues = {'R', 'K', 'D', 'E'}
    seq = sequence.upper().strip()
    N = len(seq)

    if N == 0:
        return 0.0
    
    count_charged = sum(1 for aa in seq if aa in charge_residues)

    fcr = count_charged / N
    if need_rounding:
        fcr = round(fcr, round_to)
    return fcr



if __name__ == "__main__":
	main()