# Disordered-Region-Identifier
The three-dimensional (3D) structures of proteins at the atomic level are crucial for understanding their functions and biological processes.  However, the dark proteome, which refers to regions that cannot be observed through traditional methods, presents a challenge in identifying protein structure. To address this issue, we have developed an algorithm that uses atomic coordinates to predict disordered regions within proteins. The algorithm produces a composite graph that includes a heatmap for score comparisons and a Receiver Operating Characteristic (ROC) curve. Through validation with AlphaFold scores and IUPred3 predictions, our algorithm exhibits relatively high efficiency and accuracy while implementing lightweight computations in Python.

## Installation and requirement

All functions are built on the Python 3.11.4 environment. Extra package requirements include matplotlib (3.7.1), numpy
(1.24.3), pandas (2.1.1), requests (2.31.0), retrying (1.3.4), scipy (1.10.1), seaborn (0.12.2), and scikit-learn (1.2.2). 
The sourcecode is updated to https://github.com/Sz-Jason-Chen/Disordered-Region-Identifier.

## Input and output
Input: a standard Uniprot Identifier (For example: Q00535)

Output: a composite graph presenting two distinctive panels. The left panel is the heatmap comparing the scores of input proteins in our algorithm with those obtained from AlphaFold. 

## Usage
Download source code and run main.py in the teminal or your integrated development environment (IDE). If you want to run the IUPRED module, please run iuPred.py.
