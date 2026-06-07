import requests
import pandas as pd
from io import BytesIO
from pathlib import Path


def main():
    #url = "https://alphasync.stjude.org/api/v1/protein/P31946"
    #file = requests.get(url)
    #byte_data = file.content
    #df = pd.read_json(BytesIO(byte_data))
    #print(df.head(5)["relasa"])
    #get_parse_alphasync_data(["P31946", "P04637", "A0JP26"], "")

    id = "P31946"
    #file_path = f"{id}.csv"
    #print(file_path)
    #df = pd.read_csv(file_path)
    #if 'accn' in (list(df.columns)):
        #print("ok")
    #plddt = [round(score/100, 4) for score in df["plddt"]]
    #print(plddt)
    #print(parse_alphasync_data(id, 'plddt', ""))

def parse_alphasync_data(uniprot_id: str, data_type: str, input_data_folder: str, divide_by_100: bool = False) -> list[str | float]:
    '''
    Gets the data for a given uniprot id and data type from the AlphaSync data.
     Args:
        uniprot_id (str): The uniprot id of the protein
        data_type (str): The type of data to be retrieved. It can be one of 'acc', 'site', 'afdb', 'aa', 'plddt', 'plddt10', 'asa', 'asa10', 'relasa', 'relasa10', 'dis', 'dis10', 'surf', 'surf10', 'sec', 'iso', 'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5', 'tau'
        input_data_folder (str): The path to where the AlphaSync csv files are present
        divide_by_100 (bool): Whether to divide the plddt scores by 100 or not. Default is False.


    Returns:
        list[str | float]: Returns the data for the given uniprot id and data type as a list.
    '''
    
    file_path = Path(input_data_folder).joinpath(f"{uniprot_id}.csv")
    # check if the file for the given uniprot id exists in the input data folder, if not download the data from the AlphaSync API and save it as a csv file in the input data folder
    if not file_path.is_file():
        print(f"AlphaSync file not found. Downloading for {uniprot_id}")
        get_alphasync_data([uniprot_id], input_data_folder)
    
    # read the csv file for the given uniprot id and return the data for the given data type
    df = pd.read_csv(file_path)
    if not data_type in list(df.columns):
        raise ValueError("data type not found")
    elif data_type in ["plddt", "plddt10"] and divide_by_100:
        return [round(data/100, 4) for data in df[data_type]]
    elif data_type in ["relasa", "relasa10"]:
        return [round(data, 4) for data in df[data_type]]
    elif data_type in ["sec"]:
        return ["C" if pd.isna(data) else data for data in df[data_type]]
    else:
        return [data for data in df[data_type]]
    


def get_alphasync_data(uniprot_ids: list[str], output_folder: str, supress_file_exists: bool =  True) -> None:
    '''
    Gets the AlphaSync data for all the uniprot_ids and saves each of them as a csv file

    Args:
        uniprot_id (list): The list containing uniprot ids of proteins
        output_folder (str): The path to the output folder 
        supress_file_exists (bool): Whether to suppress the message if the file already exists
    Returns:
        None: Saves the data for each uniprot id as a seperate csv file
    '''
    # get the unique uniprot ids from the list
    uniprot_ids_unique = set(uniprot_ids)

    # within a session
    with requests.Session() as session:

        # range throgh all the uniprot ids
        for id in uniprot_ids_unique:

            # only download file if it doesn't exist in the output folder
            file_path = Path(output_folder).joinpath(f"{id}.csv")
            if file_path.is_file() and not supress_file_exists:
                print(f"Alphasync details for {id} already downloaded")
            
            url = f"https://alphasync.stjude.org/api/v1/protein/{id}"

            # check if file can be accessed
            with session.get(url) as response:
                if response.status_code == 200:
                    # convert the content to a dataframe and save it as csv
                    df = pd.read_json(BytesIO(response.content))
                    df.to_csv(file_path, index = False)

                else:
                    print(f"Failed to download {id} (Status: {response.status_code})")

if __name__ == "__main__":
	main()