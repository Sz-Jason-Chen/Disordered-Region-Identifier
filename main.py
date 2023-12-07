from cores import *
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd


def main():
    protein = input('Enter protein: ')
    sequence = get_amino_sequence(protein)  # AA sequence
    print(f'Amino acid sequence: {sequence}')
    # get all structure files
    pdb_params = get_pdb_params(protein=protein)
    #print(pdb_params)
    file_paths = access_protein_files(protein=protein)

    identifiers = []  # identifier of one structure file
    gaps = []  # AA gaps
    resolutions = []  # structure measurement resolution
    valid_sections = []
    af_score = None  # AF score (100 - pLDDT)
    gaps_matrices = None
    # read every structure files
    for path in file_paths:
        # AF file
        if path.split('.')[-2][-9:] == 'alphafold':
            alphaFoldReader = AlphaFoldReader(file_path=path)
            af_score = 100 - np.array(alphaFoldReader.get_scores())  # 100 - pLDDT
        # real measurement files
        elif path.split('.')[-2][-9:] != 'alphafold':
            reader = Reader(file_path=path, params=pdb_params[path.split('.')[-2][-4:]])
            resolution = reader.get_resolution()
            # ignore structure without resolution (cannot be scored)
            if resolution:
                identifiers.append(reader.get_identifier())
                # gaps.append(reader.get_potential_gaps(sequence=sequence))
                if gaps_matrices is not None:
                    gaps_matrices = np.vstack((gaps_matrices, reader.get_potential_gaps(sequence=sequence)))
                else:
                    gaps_matrices = reader.get_potential_gaps(sequence=sequence)
                resolutions.append(resolution)
                valid_sections.append(reader.get_valid_sections())

    print(gaps_matrices)

    # testing outputs
    print('========================================================')

    # scoring
    # score_matrix = multi_chain_scoring(identifiers=identifiers, sequence=sequence, resolutions=resolutions, gaps=gaps, valid_sections=valid_sections, pdb_params=pdb_params)
    score_matrix = new_scoring(gaps_matrices=gaps_matrices)
    print(score_matrix)
    final_matrix = np.vstack((score_matrix, af_score))

    # merge the score with the AF prediction
    table = pd.DataFrame(final_matrix, index=['Score', 'AlphaFold'])


    # heatmap output
    sns.heatmap(table, cmap='inferno', vmin=0, vmax=100)
    plt.title(f'{protein}   (n={len(identifiers)})')
    plt.show()


if __name__ == '__main__':
    main()
