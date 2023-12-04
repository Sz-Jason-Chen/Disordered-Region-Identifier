import json
import os
import requests
from retrying import retry


def get_identifiers(protein):
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/{protein}'
    response = requests.get(url)

    o_dict = json.loads(response.text)
    identifiers = list(o_dict[protein]['PDB'].keys())
    return identifiers


def get_amino_sequence(protein):
    url = f'https://alphafold.ebi.ac.uk/api/prediction/{protein}'
    response = requests.get(url)
    o_dict = json.loads(response.text)

    sequence = o_dict[0]['uniprotSequence']
    return sequence

@retry(wait_fixed=1000)
def download_pdb(protein, identifier):
    file_path = f'./pdbs/{protein}_{identifier}.pdb'
    if os.path.exists(file_path):
        print(f'Checked {protein}_{identifier}.pdb')
    else:
        print(f'Downloading {protein}_{identifier}.pdb')
        url = f'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{identifier}.ent'
        response = requests.get(url)
        with open(file_path, 'wb') as f:
            f.write(response.content)

    return file_path

@retry(wait_fixed=1000)
def download_alphafold_pdb(protein):
    file_path = f'./pdbs/{protein}_alphafold.pdb'
    if not os.path.exists(file_path):
        print(f'Downloading {protein}_alphafold.pdb')
        data_url = f'https://alphafold.ebi.ac.uk/api/prediction/{protein}'
        data_response = requests.get(data_url)
        o_dict = json.loads(data_response.text)
        download_url = o_dict[0]['pdbUrl']
        download_response = requests.get(download_url)

        with open(file_path, 'wb') as f:
            f.write(download_response.content)
    return file_path


def access_protein_files(protein):
    file_paths = []
    identifiers = get_identifiers(protein=protein)
    for identifier in identifiers:
        file_paths.append(download_pdb(protein=protein, identifier=identifier))
    file_paths.append(download_alphafold_pdb(protein=protein))
    return file_paths


if __name__ == '__main__':
    access_protein_files('Q00535')
