import requests
from pathlib import Path
import torch
import pydssp
from typing import Literal

def main():
    print("Test")
    #download_alphafold_structures(["F4HVG8","Q9I1F6"], "")
    #line = "ATOM     20  CB  LEU A   3     -33.929  14.913   6.050  1.00 25.50           C  "
    #print(float(line.split()[10])/100)
    # Sample coordinates
    #batch, length, atoms, xyz = 10, 100, 4, 3
    ## atoms should be 4 (N, CA, C, O) or 5 (N, CA, C, O, H)
    #coord = torch.randn([batch, length, atoms, xyz]) # batch-dim is optional
    #print(coord.dim())
    #print(get_plddt(["F4HVG8","Q9I1F6"],""))
    #parsed_pdb = parse_pdb("F4HVG8", "")
    #residue_num_list = parsed_pdb["residue_num_list"]
    #atom_name_list = parsed_pdb["atom_name_list"]
    #x_cord_list = parsed_pdb["x_cord_list"]
    #y_cord_list = parsed_pdb["y_cord_list"]
    #z_cord_list = parsed_pdb["z_cord_list"]
    #raw_coord = get_atom_cord(atom_name_list, x_cord_list, y_cord_list, z_cord_list)[0:16]
    #print(raw_coord)
    #cord = torch.tensor(raw_coord)
    #structured_cord = cord.view(4, 4, 3)
    #print(structured_cord)
    #print(get_secondary_structure(["F4HVG8","Q9I1F6"], "", output_type="index"))

    




def download_alphafold_structures(uniprot_ids: list[str], output_folder: str) -> None:
    '''
    Gets all the protein strcutures from the AlphaFold Protein Structure Database based on the list of uniprot ids and saves it.

    Args:
        uniprot_id (list): The list containing uniprot ids of proteins
        output_folder (str): The path to the output folder
    '''

    with requests.Session() as session:

        # range throgh all the uniprot ids
        for id in uniprot_ids:

            # only download file if it doesn't exist in the output folder
            file_path = Path(output_folder).joinpath(f"AF-{id}-F1-model_v6.pdb")
            if file_path.is_file():
                print(f"Structure for {id} already exist")
            
            url = f"https://alphafold.ebi.ac.uk/files/AF-{id}-F1-model_v6.pdb"

            with session.get(url) as response:
                if response.status_code == 200:
                    with open(file_path, 'wb') as file:
                        file.write(response.content)
                else:
                    print(f"Failed to download {id} (Status: {response.status_code})")

def get_plddt(uniprot_ids: list[str], input_folder: str) -> dict[str, list[float]]:
    '''
    Gets the scaled plddt scores (plddt_score/100) for each of the protein sequence with the given Uniprot ID.

    Args:
        uniprot_ids (list): the list containing uniprot ids of proteins
        input_folder (str): the path to where the AlphaFold PDB files are present
    
    Returns:
        dict[str, list[float]]: Returns the the plddt scores for each protein sequence as,
          {'uniprot_id': [plddt score for each residue]}
    '''

    plddt_scores: dict[str,list[float]] = {}

    # range over all the protein structures
    for id in uniprot_ids:
        # parse the pdb file for each uniprot id
        parsed_pdb = parse_pdb(id, input_folder)

        # store the scaled plddt scores for a given protein with a uniprot id.
        plddt_score = [round(float(score)/100, 4) for score in parsed_pdb["plddt_list"]]

        # add to the dictionary the plddt scores for each of the protein.
        plddt_scores[id] = plddt_score
    
    return plddt_scores


def get_secondary_structure(uniprot_ids: list[str], input_folder: str, output_type: Literal['onehot', 'index', 'c3'] = 'c3') -> dict[str, list[str]]:
    '''
    Gets the secondary structure for each of the protein sequence with the given Uniprot ID.
    Args:
        uniprot_ids (list): the list containing uniprot ids of proteins
        input_folder (str): the path to where the AlphaFold PDB files are present
        output_type (Literal['onehot', 'index', 'c3']): the type of output format for the secondary structure, default is 'c3'
    Returns:
        dict[str, list[str]]: Returns the the secondary structure for each protein sequence as,
          {'uniprot_id': [secondary structure for each residue]}
    '''

    dssp = {}

    # range over all the protein structures
    for id in uniprot_ids:
        # parse the pdb file for each uniprot id and get the coordinates and the type of atoms and number of residues
        parsed_pdb = parse_pdb(id, input_folder)
        residue_num_list = parsed_pdb["residue_num_list"]
        atom_name_list = parsed_pdb["atom_name_list"]
        x_cord_list = parsed_pdb["x_cord_list"]
        y_cord_list = parsed_pdb["y_cord_list"]
        z_cord_list = parsed_pdb["z_cord_list"]
        # length of the protein sequence
        num_res = len(residue_num_list)
        num_atoms = 4
        # get the xyz coordinates for the required atoms i.e, N, CA, C, O for all the residues
        raw_cord: list[list[float]] = get_atom_cord(atom_name_list, x_cord_list, y_cord_list, z_cord_list)
        xyz_dim = 3
        # transform the coordinate into suitable format for input to pydssp
        cord = torch.tensor(raw_cord)
        structured_cord = cord.view(num_res,num_atoms, xyz_dim)
        # get the secondary structure
        dssp[id] = pydssp.assign(structured_cord, out_type=output_type)
    
    return dssp

def get_atom_cord(atom_name_list: list[str], x_cord_list: list[float], y_cord_list: list[float], z_cord_list: list[float]) -> list[list[float]]:
    '''
    Gets the coordinates in the x,y,z plane for the atoms - "N", "CA", "C", "O" for each residue of the protein

    Args:
        atom_name_list (list): the list containing the type of atoms for each residue
        x_cord_list (list): the list containing the x coordinates for each atom
        y_cord_list (list): the list containing the y coordinates for each atom
        z_cord_list (list): the list containing the z coordinates for each atom
    
    Returns:
        list[list[float]]: Returns the list of coordinates for the specified atoms as, [[x1,y1,z1], [x2,y2,z2], ...]
    '''

    cord_list = []
    for atom in range(len(atom_name_list)):
        # gets the [x,y,z] coordinates for the specified atoms i.e, N, CA, C, O for each residue
        if atom_name_list[atom] in ["N", "CA", "C", "O"]:
            x = x_cord_list[atom]
            y = y_cord_list[atom]
            z = z_cord_list[atom]
            cord_list.append([x, y, z])
    
    return cord_list

def parse_pdb(uniprot_id: str,input_folder: str) -> dict[str, list]:
    '''
    Parses the AlphaFold PDB file and extracts the atom information, residue information, coordinates and plddt scores.
    
    Args:
        input_folder (str): the path to where the AlphaFold PDB files are present
        uniprot_id (str): the uniprot id of the protein structure to be parsed
    Returns:
        dict[str, list[str]]: Returns a dictionary containing the atom information, residue information, coordinates and plddt scores as,
        {"atom_name_list": [list of atom names], "atom_num_list": [list of atom numbers], "chain_id_list": [list of chain ids], "x_cord_list": [list of x coordinates], "y_cord_list": [list of y coordinates], "z_cord_list": [list of z coordinates], "residue_name_list": [list of residue names], "residue_num_list": [list of residue numbers], "plddt_list": [list of plddt scores]}
    '''
    
    file_path = Path(input_folder).joinpath(f"AF-{uniprot_id}-F1-model_v6.pdb")
    
    
    # check if file exist, if not, download it
    if not file_path.is_file():
        print(f"AlphaFold PDB file not found. Downloading for {uniprot_id}")
        download_alphafold_structures([uniprot_id],input_folder)
    
    # Initialize empty lists for structured collection
    data = {
        "atom_name_list": [], "atom_num_list": [],
        "x_cord_list": [], "y_cord_list": [], "z_cord_list": [],
        "residue_name_list": [], "residue_num_list": [], "plddt_list": []
    }


    # open the pdb file read through each line
    with open(file_path, "r") as file:
        for line in file:
            # get all the atom information
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                
                # Append all-atom data
                data["atom_name_list"].append(atom_name)
                data["atom_num_list"].append(int(line[6:11]))
                data["x_cord_list"].append(float(line[30:38]))
                data["y_cord_list"].append(float(line[38:46]))
                data["z_cord_list"].append(float(line[46:54]))
                
                # Extract residue info only if it matches AlphaFold's per-residue CA backbone
                if atom_name == "CA":
                    data["residue_name_list"].append(line[17:20].strip())
                    data["residue_num_list"].append(int(line[22:26]))
                    data["plddt_list"].append(float(line[60:66]))

    return data


def three_to_one_amino_acid_residue(three_letter):

    amino_acids = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
        "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
        "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
        "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
    }

    return amino_acids.get(three_letter.upper(), "X")


if __name__ == "__main__":
	main()