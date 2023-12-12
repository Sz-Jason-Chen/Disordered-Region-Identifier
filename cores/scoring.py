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


def scoring(gaps_matrices):
    all_count = len(gaps_matrices)
    for row in range(len(gaps_matrices)):
        for col in range(len(gaps_matrices[0])):
            if gaps_matrices[row, col] != 0:
                gaps_matrices[row, col] = 1 / gaps_matrices[row, col]
    # print('gaps_matrices',gaps_matrices)
    score_matrix = np.zeros((1, len(gaps_matrices[0])))

    for col in range(len(gaps_matrices[0])):
        score = 0
        valid_count = 0
        for row in range(len(gaps_matrices)):
            score += gaps_matrices[row, col]
            if gaps_matrices[row, col] != 0:
                valid_count += 1

        abs_scores = [abs(x) for x in gaps_matrices[:, col]]
        max_score = sum(abs_scores)
        min_score = -max_score
        valid_percentage = valid_count / len(gaps_matrices)

        if valid_count == 0:
            #score_matrix[0, col] = 75
            score_matrix[0, col] = (-50)*exp((-0.02) * all_count) + 100
        else:
            scaled_score = min_max_scaler(x=score, min_v=min_score, max_v=max_score) * 100
            score_matrix[0, col] = (100 - scaled_score) * ((valid_percentage - 1) ** 8) * (exp(-0.5 * valid_count)) + scaled_score
            # score_matrix[0, col] = (100 - scaled_score) * ((valid_percentage - 1) ** 8) * (exp(-valid_count)) + scaled_score * (1 - exp(-2 * all_count))
            # min = score if a position is not in a gap in all structures
    # max = score if a position is in a gap in all structures

    return score_matrix


# test module
if __name__ == '__main__':
    pass
