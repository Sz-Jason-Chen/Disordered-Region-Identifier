from math import exp
import numpy as np

test_data = [{'start': 11, 'end': 14, 'resolution': 2.65},
             {'start': 11, 'end': 14, 'resolution': 2.30},
             {'start': 11, 'end': 14, 'resolution': 2.35},
             {'start': None, 'end': None, 'resolution': 2.20},
             {'start': 11, 'end': 15, 'resolution': 1.95},
             {'start': None, 'end': None, 'resolution': 1.90},
             {'start': 13, 'end': 15, 'resolution': 2.09},
             {'start': 9, 'end': 17, 'resolution': 2.91},
             {'start': 11, 'end': 17, 'resolution': 2.55},
             {'start': 8, 'end': 15, 'resolution': 3.05}]


def sigmoid(x):
    return 1 / (1 + exp(-x))


def tanh(x):
    return (exp(x) - exp(-x)) / (exp(x) + exp(-x))


def min_max_scaler(x, min_v, max_v):
    return (x - min_v) / (max_v - min_v)


def scoring(identifiers, sequence, resolutions, gaps):
    reso_matrix = np.zeros((len(identifiers), len(sequence)))
    for row in range(len(identifiers)):
        for col in range(len(sequence)):
            reso_matrix[row, col] = -resolutions[row]

        for gap in gaps[row]:
            # print(identifiers[row], gap)
            start = gap[0] - 1
            end = min(gap[1], len(sequence))  # P04439_6apn.pdb
            for col in range(start, end):
                reso_matrix[row, col] = resolutions[row]

    print(reso_matrix[:, 159])

    score_matrix = np.zeros((1, len(sequence)))
    for col in range(len(sequence)):
        score = 0
        for row in range(len(identifiers)):
            score += 1 / reso_matrix[row, col]
        score_matrix[0, col] = score
    min_score = 0
    for resolution in resolutions:
        min_score -= 1 / resolution
    max_score = -min_score
    print(max_score, min_score)
    for col in range(len(sequence)):
        score_matrix[0, col] = min_max_scaler(x=score_matrix[0, col], min_v=min_score, max_v=max_score) * 100
    return score_matrix


def legacy_scoring(test_data):
    start = {}
    end = {}
    start_c = {}
    end_c = {}

    length = len(test_data)
    for data in test_data:
        if data['start'] in start:
            start[data['start']] += 1 / data['resolution']
            start_c[data['start']] += 1
        else:
            start[data['start']] = 1 / data['resolution']
            start_c[data['start']] = 1

        if data['end'] in end:
            end[data['end']] += 1 / data['resolution']
            end_c[data['end']] += 1
        else:
            end[data['end']] = 1 / data['resolution']
            end_c[data['end']] = 1

    for key, value in start.items():
        start[key] = start[key] * start_c[key] / length
        start[key] = tanh(start[key])
    for key, value in end.items():
        end[key] = end[key] * end_c[key] / length
        end[key] = tanh(end[key])

    print(start)
    print(end)

    all = {}

    if None in start and None in end:
        all[None] = round((start[None] + end[None]) / 2, 3)
        start.pop(None)
        end.pop(None)
    elif None in start:
        start.pop(None)
    elif None in end:
        end.pop(None)

    for key_s, value_s in start.items():
        for key_e, value_e in end.items():
            all[f'{key_s}-{key_e}'] = round((value_s + value_e) / 2, 3)
    print(all)


if __name__ == '__main__':
    legacy_scoring(test_data)
