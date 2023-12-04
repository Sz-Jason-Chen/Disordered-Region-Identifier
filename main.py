from cores import *
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd


def main():
    protein = input('Enter protein: ')
    sequence = get_amino_sequence(protein)
    file_paths = access_protein_files(protein=protein)
    # print(file_paths)

    identifiers = []
    gaps = []
    resolutions = []
    af_score = None
    for path in file_paths:
        if path.split('.')[-2][-9:] == 'alphafold':
            alphaFoldReader = AlphaFoldReader(file_path=path)
            af_score = 100 - np.array(alphaFoldReader.get_scores())
        elif path.split('.')[-2][-9:] != 'alphafold':
            reader = Reader(file_path=path)
            resolution = reader.get_resolution()
            if resolution:
                identifiers.append(reader.get_identifier())
                gaps.append(reader.get_potential_gaps())
                resolutions.append(resolution)

    # print(identifiers)
    print(gaps)
    # print(resolutions)
    # print(af_score)

    print(f'Amino acid sequence: {sequence}')
    score_matrix = scoring(identifiers=identifiers, sequence=sequence, resolutions=resolutions, gaps=gaps)

    final_matrix = np.vstack((score_matrix, af_score))
    table = pd.DataFrame(final_matrix, index=['Score', 'Alphafold'])
    sns.heatmap(table, cmap='inferno', vmin=0, vmax=100)
    plt.show()


if __name__ == '__main__':
    main()
