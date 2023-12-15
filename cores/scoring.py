from math import exp
import numpy as np


def min_max_scaler(x: float, min_v: float, max_v: float) -> float:
    """
    Min-max scaling
    :param x: val (score)
    :param min_v: min
    :param max_v: max
    :return: scaled value (scaled score)
    """
    return (x - min_v) / (max_v - min_v)


def scoring(gaps_matrices: np.ndarray) -> np.ndarray:
    """
    Scoring by resolution, valid sequence count and percentage
    :param gaps_matrices: a matrix containing all gaps in all PDB sequences
    :return: a (1, len(protein_sequence)) matrix, row 0 is the score of each site
    """
    all_count = len(gaps_matrices)  # valid sequences number
    # take PDB's resolution into consideration, scored by sites
    for row in range(len(gaps_matrices)):
        for col in range(len(gaps_matrices[0])):
            if gaps_matrices[row, col] != 0:
                gaps_matrices[row, col] = 1 / gaps_matrices[row, col]

    # integrate sites in sequences
    score_matrix = np.zeros((1, len(gaps_matrices[0])))
    for col in range(len(gaps_matrices[0])):
        score = 0  # the site's score
        valid_count = 0
        for row in range(len(gaps_matrices)):
            score += gaps_matrices[row, col]
            if gaps_matrices[row, col] != 0:
                valid_count += 1

        # scaling and correction
        # calculate the possible min and max score of this site
        abs_scores = [abs(x) for x in gaps_matrices[:, col]]
        max_score = sum(abs_scores)
        min_score = -max_score
        valid_percentage = valid_count / len(gaps_matrices)  # the valid sequences percentage of this site
        # missing site
        if valid_count == 0:
            score_matrix[0, col] = (-50)*exp((-0.02) * all_count) + 100  # correction
        # not missing site
        else:
            scaled_score = min_max_scaler(x=score, min_v=min_score, max_v=max_score) * 100  # scaling
            score_matrix[0, col] = (100 - scaled_score) * ((valid_percentage - 1) ** 8) * (exp(-0.5 * valid_count)) + scaled_score  # correction

    return score_matrix
