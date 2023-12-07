from math import exp
from typing import List, Tuple
import numpy as np
import pandas as pd


def sigmoid(x):
    """
    sigmoid normalization
    :param x: val
    :return: sigmoid(val)
    """
    return 1 / (1 + exp(-x))


def tanh(x):
    """
    tanh normalization
    :param x: val
    :return: tanh(val)
    """
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x))


def min_max_scaler(x, min_v, max_v):
    """
    min max normalization
    :param x: val
    :param min_v: min
    :param max_v: max
    :return: (val - min) / (max - min)
    """
    return (x - min_v) / (max_v - min_v)


def single_chain_scoring(identifiers: List[str],
                         sequence: str,
                         resolutions: List[float],
                         valid_sections,
                         gaps,
                         pdb_params) -> np.ndarray:
    """
    Main scoring algorithm
    :param pdb_params:
    :param identifiers: all identifiers list
    :param sequence: protein AA sequence
    :param resolutions: all mesures' resolutions list
    :param gaps: potential gaps list, gaps[structure[gap(start, end)]]
    :return: gap score of all position
    """
    # load resolution and gap data in a matrix
    # row nums is the identifications nums, col nums is the aa sequence length
    reso_matrix = np.zeros((len(identifiers), len(sequence)))
    for row in range(len(identifiers)):
        # no gap position's score is -resolution
        for valid_section in valid_sections[row]:
            print(valid_section)
            for col in range(valid_section[0]-1, valid_section[1]):
                reso_matrix[row, col] = -resolutions[row]
        # gap position's score is resolution
        for gap in gaps[row]:
            # print(identifiers[row], gap)
            start = gap[0] - 1
            end = min(gap[1], len(sequence))  # P04439_6apn.pdb
            for col in range(start, end):
                reso_matrix[row, col] = resolutions[row]

    # init score matrix as zero matrix
    # row nums is 1, col nums is the aa sequence length (corresponding to reso matrix)
    score_matrix = np.zeros((1, len(sequence)))
    # position score = sum((-1 if no gap | 1 if gap) / resolution) = sum(reso matrix)
    for col in range(len(sequence)):
        score = 0
        for row in range(len(identifiers)):
            if reso_matrix[row, col] > 0:
                score += 1 / reso_matrix[row, col]
            elif reso_matrix[row, col] < 0:
                score += 1 / reso_matrix[row, col] - 1
        score_matrix[0, col] = score
    # min = score if a position is not in a gap in all structures
    # max = score if a position is in a gap in all structures
    min_score = 0
    for resolution in resolutions:
        min_score -= 1 / resolution
    max_score = -min_score
    # normalization
    for col in range(len(sequence)):
        score_matrix[0, col] = min_max_scaler(x=score_matrix[0, col], min_v=min_score, max_v=max_score)
        score_matrix[0, col] = tanh(score_matrix[0, col])
        score_matrix[0, col] = score_matrix[0, col] * 100
    return score_matrix

def new_scoring(gaps_matrices):
    for row in range(len(gaps_matrices)):
        for col in range(len(gaps_matrices[0])):
            if gaps_matrices[row, col] != 0:
                gaps_matrices[row, col] = 1 / gaps_matrices[row, col]
    print('gaps_matrices',gaps_matrices)
    score_matrix = np.zeros((1, len(gaps_matrices[0])))
    # position score = sum((-1 if no gap | 1 if gap) / resolution) = sum(reso matrix)
    for col in range(len(gaps_matrices[0])):
        score = 0
        valid_count = 0

        for row in range(len(gaps_matrices)):
            score += gaps_matrices[row, col]
            if gaps_matrices[row, col] != 0:
                valid_count += 1

        abs_scores = [abs(x) for x in gaps_matrices[:,col]]
        max_score = sum(abs_scores)
        min_score = -max_score
        valid_percentage = valid_count / len(gaps_matrices)

        if valid_count == 0:
            score_matrix[0, col] = 100
        else:
            scaled_score = min_max_scaler(x=score, min_v=min_score, max_v=max_score) * 100
            score_matrix[0, col] = (100 - scaled_score) * ((valid_percentage - 1) ** 8) * (exp(-0.5 * valid_count)) + scaled_score
    # min = score if a position is not in a gap in all structures
    # max = score if a position is in a gap in all structures

    return score_matrix


def best_scs(identifiers: List[str],
                         sequence: str,
                         resolutions: List[float],
                         valid_sections,
                         gaps,
                         pdb_params) -> np.ndarray:
    """
    Main scoring algorithm
    :param pdb_params:
    :param identifiers: all identifiers list
    :param sequence: protein AA sequence
    :param resolutions: all mesures' resolutions list
    :param gaps: potential gaps list, gaps[structure[gap(start, end)]]
    :return: gap score of all position
    """
    # load resolution and gap data in a matrix
    # row nums is the identifications nums, col nums is the aa sequence length
    reso_matrix = np.zeros((len(identifiers), len(sequence)))
    for row in range(len(identifiers)):
        # no gap position's score is -resolution
        for col in range(len(sequence)):
            reso_matrix[row, col] = -resolutions[row]
        # gap position's score is resolution
        for gap in gaps[row]:
            # print(identifiers[row], gap)
            start = gap[0] - 1
            end = min(gap[1], len(sequence))  # P04439_6apn.pdb
            for col in range(start, end):
                reso_matrix[row, col] = resolutions[row]

    # init score matrix as zero matrix
    # row nums is 1, col nums is the aa sequence length (corresponding to reso matrix)
    score_matrix = np.zeros((1, len(sequence)))
    # position score = sum((-1 if no gap | 1 if gap) / resolution) = sum(reso matrix)
    for col in range(len(sequence)):
        score = 0
        for row in range(len(identifiers)):
            if reso_matrix[row, col] > 0:
                score += 1 / reso_matrix[row, col]
            elif reso_matrix[row, col] < 0:
                score += 1 / reso_matrix[row, col]
        score_matrix[0, col] = score
    # min = score if a position is not in a gap in all structures
    # max = score if a position is in a gap in all structures
    min_score = 0
    for resolution in resolutions:
        min_score -= 1 / resolution
    max_score = -min_score
    # normalization
    for col in range(len(sequence)):
        score_matrix[0, col] = min_max_scaler(x=score_matrix[0, col], min_v=min_score, max_v=max_score)
        score_matrix[0, col] = tanh(score_matrix[0, col])
        score_matrix[0, col] = score_matrix[0, col] * 100
    return score_matrix


def multi_chain_scoring(identifiers: List[str],
                        sequence: str,
                        resolutions: List[float],
                        valid_sections,
                        gaps,
                        pdb_params):
    total_score = pd.DataFrame()

    chain_gaps = {}
    chain_valid_sections = {}
    for key, value in pdb_params.items():
        for section in value:
            print(section)
            if section['chain'] not in chain_gaps:
                chain_gaps[section['chain']] = []
                chain_valid_sections[section['chain']] = []

    chain_gaps = dict(sorted(chain_gaps.items()))
    chain_valid_sections = dict(sorted(chain_valid_sections.items()))

    for pdb in gaps:
        print('pdb',pdb)
        for key in chain_gaps.keys():
            if key in pdb.keys():
                chain_gaps[key].append(pdb[key])
            else:
                chain_gaps[key].append([])


    for pdb in valid_sections:
        print('pdb', pdb)
        for key in chain_valid_sections.keys():
            if key in pdb.keys():
                chain_valid_sections[key].append(pdb[key]['region'])
            else:
                chain_valid_sections[key].append([])


    print('gaps', gaps)
    print('chain_gaps', chain_gaps)
    print('chain_valid_sections', chain_valid_sections)


    for key, value in chain_gaps.items():
        chain_score = single_chain_scoring(identifiers, sequence, resolutions, chain_valid_sections[key], value, pdb_params)
        #chain_score = best_scs(identifiers, sequence, resolutions, chain_valid_sections[key], value, pdb_params)
        # print(key, chain_score)
        chain_score = pd.DataFrame(chain_score, index=[key])
        total_score = pd.concat([total_score, chain_score])
        print(total_score)

    return total_score

# test module
if __name__ == '__main__':
    pass

