import csv
import matplotlib.pyplot as plt

scores = []
afs = []
valid_score = 0
valid_af = 0
protein = 0
with open(r'supp\results.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        protein += 1
        scores += eval(row[1])
        afs += eval(row[2])
        if len(eval(row[1])) > 0:
            valid_score += 1
        if len(eval(row[2])) > 0:
            valid_af += 1
print('proteins: ',protein)
print('valid score: ', valid_score)
print('valid af: ', valid_af)
print('amino acids: ', len(afs))
print('scores: ', len(scores))

plt.hist(scores, bins=100, edgecolor='black')
plt.title('Histogram of Scores')
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.show()

scores = [value for value in scores if value > 1 and value != 100]
# scores = [value for value in scores if value > 49.5 and value < 50.5]

plt.hist(scores, bins=100, edgecolor='black')
plt.title('Histogram of Scores (Without 0 and 100)')
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.show()

plt.hist(afs, bins=100, edgecolor='black')
plt.title('Histogram of AlphaFold')
plt.xlabel('Score')
plt.ylabel('Frequency')
plt.show()