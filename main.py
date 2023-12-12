from cores import *
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import seaborn as sns
import pandas as pd
from sklearn.metrics import roc_curve, auc



def main():
    print("█▓▒░ D I S O R D E R E D - R E G I O N - I D E N T I F I E R ░▒▓█")
    print("By Group 2")
    print('========================================')
    protein = input('Enter protein: ')

    print('----------------------------------------')
    print('Accessing sources')
    print('Reminder:')
    print('This section will access the Uniprot server and will automatically retry if access fails.')
    print('If it stuck for more than 1 minute, please restart the programme and check the internet connection.')
    sequence = get_amino_sequence(protein)  # AA sequence
    if sequence is None:
        print('No available protein sequence.')
        return -1
    print(f'Amino acid sequence length: {len(sequence)}')
    # get all structure files
    print('Accessing PDB files')
    pdb_params = get_pdb_params(protein=protein)
    file_paths = access_protein_files(protein=protein)
    print('File validation complete')

    print('----------------------------------------')
    print('Identifying gaps')
    valid_sections = []
    af_score = None  # AF score (100 - pLDDT)
    gaps_matrices = None
    # read every structure files
    for path in file_paths:
        # AF file
        if path.split('.')[-2][-9:] == 'alphafold':
            print(f'Resolving AlphaFold file')
            alpha_fold_reader = AlphaFoldReader(file_path=path)
            af_score = 100 - np.array(alpha_fold_reader.get_scores())  # 100 - pLDDT
        # real measurement files
        elif path.split('.')[-2][-9:] != 'alphafold':
            print(f'Examining {path[-8:-4].upper()} file')
            reader = Reader(file_path=path, params=pdb_params[path.split('.')[-2][-4:]])
            resolution = reader.get_resolution()
            # ignore structure without resolution (cannot be scored)
            if resolution:
                gap_array = reader.get_potential_gaps(sequence=sequence)
                if gaps_matrices is not None:
                    gaps_matrices = np.vstack((gaps_matrices, gap_array))
                else:
                    gaps_matrices = gap_array
                valid_sections.append(reader.get_valid_sections())

    #print(gaps_matrices)

    # testing outputs
    print('----------------------------------------')
    print('Scoring')
    if gaps_matrices is None:
        print('No available PDB file for gap identification.')
        return 1
    # scoring
    # score_matrix = multi_chain_scoring(identifiers=identifiers, sequence=sequence, resolutions=resolutions, gaps=gaps, valid_sections=valid_sections, pdb_params=pdb_params)
    score_matrix = scoring(gaps_matrices=gaps_matrices)
    #print(score_matrix)


    print('----------------------------------------')
    print('Result')
    """location = np.where(score_matrix > 50)[0]
    print(location)
    """
    # merge the score with the AF prediction
    final_matrix = np.vstack((score_matrix, af_score))
    table = pd.DataFrame(final_matrix, index=['Score', 'AlphaFold'])
    # Convert table from index 0 to index 1
    table = table.rename(columns=dict(zip(table.columns, table.columns + 1)))
    print(table)

    gap_positions = []
    start = 0
    avg_score = 0
    af_avg_score = 0
    for column_name, value in table.loc['Score'].items():
        if value >= 50 and start == 0:
            start = column_name
            avg_score = value
            af_avg_score = table.loc['AlphaFold', column_name]
        elif value >= 50 and start != 0 and column_name != len(sequence):
            avg_score += value
            af_avg_score += table.loc['AlphaFold', column_name]
        elif value < 50 and start != 0:
            end = column_name - 1
            avg_score /= end - start + 1
            af_avg_score /= end - start + 1
            gap_positions.append({'start': start, 'end': end, 'avg_score': avg_score, 'af_avg_score':af_avg_score})
            start = 0
            af_avg_score = 0
        elif value >= 50 and start != 0 and column_name == len(sequence):
            end = len(sequence)
            avg_score += value
            af_avg_score += table.loc['AlphaFold', column_name]
            avg_score /= end - start + 1
            af_avg_score /= end - start + 1
            gap_positions.append({'start': start, 'end': end, 'avg_score': avg_score, 'af_avg_score':af_avg_score})
        elif value >= 50 and start == 0 and column_name == len(sequence):
            gap_positions.append({'start': len(sequence), 'end': len(sequence), 'avg_score': value, 'af_avg_score':table.loc['AlphaFold', column_name]})
    # print(gap_positions)



    if len(gap_positions) > 0:
        print(f'{len(gap_positions)} gaps are identified.')
        print('The gaps and their average position scores are listed:')
        for gap in gap_positions:
            print(f'Start: {gap["start"]},\tEnd: {gap["end"]},\tAvg score: {round(gap["avg_score"], 2)},\tAF avg score: {round(gap["af_avg_score"], 2)}')
    else:
        print('No gap identified')

    print('The heatmap presents each position\'s score, and the gap regions are marked with red boxes.')
    # heatmap output
    sns.heatmap(table, cmap='viridis', vmin=0, vmax=100)
    for gap in gap_positions:
        rectangle = Rectangle((gap["start"], 0), width=(gap["end"]-gap["start"]+1), height=2, linewidth=1, edgecolor='red', facecolor='none')
        plt.gca().add_patch(rectangle)
    plt.title(f'{protein}   (n={len(gaps_matrices)})')
    plt.show()





    ###Draw ROC curve to estimate the efficiency of our model
    # extract Score and AlphaFold score
    score_column = table.loc['Score']
    alphafold_score_column = table.loc['AlphaFold']

    # use AlphaFold_Score as True Class Labels，Score as the input of models, 30 is determined by 1-70(70 is the threshold of high evidence in pLDDT)
    y_true = (alphafold_score_column > 30).astype(int)
    y_scores = score_column

    # use roc_curve to calculate the roc results
    fpr, tpr, thresholds = roc_curve(y_true, y_scores)
    roc_auc = auc(fpr, tpr)

    # draw ROC curve
    plt.figure(figsize=(8, 8))
    plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'AUC = {roc_auc:.2f}')
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Random', alpha=0.5)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    return 0

if __name__ == '__main__':
    main()
