import json
import os
import requests
from retrying import retry


@retry(wait_fixed=1000)
def get_pdb_params(protein):
    """

    :param protein:
    :return:
    """
    url = f'https://www.ebi.ac.uk/pdbe/api/mappings/{protein}'
    response = requests.get(url)

    o_dict = json.loads(response.text)
    params = {}
    for key, value in o_dict[protein]['PDB'].items():
        sections = []
        for section in o_dict[protein]['PDB'][key]:
            chain = section['chain_id']
            unp_start = section['unp_start']
            unp_end = section['unp_end']
            if section['start']['author_residue_number']:
                residue_start = section['start']['author_residue_number']
            else:
                residue_start = section['start']['residue_number']
            if section['end']['author_residue_number']:
                residue_end = section['end']['author_residue_number']
            else:
                residue_end = residue_start + unp_end - unp_start
            sections.append({'chain': chain, 'unp_start': unp_start, 'unp_end': unp_end, 'residue_start': residue_start,
                             'residue_end': residue_end})
        params[key] = sections

    # return pdb_paras{identifier:section[{chain, unp_start, unp_end, residue_start, residue_end}]}
    # eg. pdb_paras{'5i2n':['chain':'B', 'unp_start':394, 'unp_end':544, 'residue_start':3, 'residue_end':153]}
    return params


@retry(wait_fixed=1000)
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
        print(f'Downloading {protein}_{identifier}.pdb', end='   ')
        url = f'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{identifier}.ent'
        response = requests.get(url)
        with open(file_path, 'wb') as f:
            f.write(response.content)
        print('Done')

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
    identifiers = get_pdb_params(protein=protein).keys()
    for identifier in identifiers:
        file_paths.append(download_pdb(protein=protein, identifier=identifier))
    file_paths.append(download_alphafold_pdb(protein=protein))
    return file_paths


if __name__ == '__main__':
    # access_protein_files('Q00535')
    print(get_pdb_params('Q00535'))
