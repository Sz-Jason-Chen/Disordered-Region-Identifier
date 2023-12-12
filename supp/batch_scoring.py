import csv
import os
import sys
sys.path.append('.')
import threading
import time

from cores import *
import numpy as np
import pandas as pd

class BatchScoring:
    def __init__(self) -> None:
        pass

    def score_module(self, protein):
        sequence = get_amino_sequence(protein)  # AA sequence
        if sequence is None:
            return {'protein': protein, 'score': [], 'alphafold': []}
        pdb_params = get_pdb_params(protein=protein)
        file_paths = access_protein_files(protein=protein, folder='supp/pdbs/')
        # print(file_paths)
        
        valid_sections = []
        af_score = None  # AF score (100 - pLDDT)
        gaps_matrices = None
        # read every structure files
        for path in file_paths:
            # AF file
            if path.split('.')[-2][-9:] == 'alphafold':
                alpha_fold_reader = AlphaFoldReader(file_path=path)
                af_score = 100 - np.array(alpha_fold_reader.get_scores())  # 100 - pLDDT

            # real measurement files
            elif path.split('.')[-2][-9:] != 'alphafold':
                reader = Reader(file_path=path, params=pdb_params[path.split('.')[-2][-4:]])
                resolution = reader.get_resolution()
                # ignore structure without resolution (cannot be scored)
                if resolution:
                    if gaps_matrices is not None:
                        gaps_matrices = np.vstack((gaps_matrices, reader.get_potential_gaps(sequence=sequence)))
                    else:
                        gaps_matrices = reader.get_potential_gaps(sequence=sequence)
                    valid_sections.append(reader.get_valid_sections())
            os.remove(path)

        af_list = [round(s, 2) for s in af_score]
        if gaps_matrices is None:
            result = {'protein': protein, 'score': [], 'alphafold': af_list}
        else:
            score_matrix = scoring(gaps_matrices=gaps_matrices)
            score_list = [round(score_matrix[0, col], 2) for col in range(len(score_matrix[0]))]
            result = {'protein':protein, 'score':score_list, 'alphafold': af_list}

        
        return result

    def main(self):
        scored_proteins = []
        with open(r'supp\results.csv', 'r') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                scored_proteins.append(row[0])
        print(len(scored_proteins))

        not_scored_proteins = []
        with open(r'supp\human_proteins.csv', 'r') as f:
            reader = csv.reader(f)
            next(reader)
            for row in reader:
                if row[0] not in scored_proteins:
                    not_scored_proteins.append(row[0])
        print(len(not_scored_proteins))
        
        pool = 1
        pooled_proteins = [not_scored_proteins[i:i + pool] for i in range(0, len(not_scored_proteins), pool)]
        self.results = [None] * pool

        def single_thread(protein, index):
            result = self.score_module(protein)
            self.results[index] = result

        for p_proteins in pooled_proteins: 
            threads = [None] * pool
            for i in range(pool):
                threads[i] = threading.Thread(target=single_thread, kwargs={'protein': p_proteins[i], 'index':i})
                threads[i].start()
            while True:
                if None not in self.results:
                    with open(r'supp\results.csv', 'a', newline='') as f:
                        writer = csv.DictWriter(f, fieldnames=['protein', 'score', 'alphafold'])
                        for result in self.results:
                            writer.writerow(result)
                    self.results = [None] * pool
                    break
                else:
                    time.sleep(1)

                

if __name__=='__main__':
    batch_scoring = BatchScoring()
    batch_scoring.main()

