class Reader:
    def __init__(self, file_path):
        self.file_path = file_path

    def get_identifier(self):
        with open(self.file_path, 'r') as file:
            first_line = file.readline()
            first_line = first_line.rstrip()
            return first_line[-4:]

    def get_potential_gaps(self):
        atom_amino_acid_coordinates = []
        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    # extract the atom coordinate, amino acid and amino acid coordinate
                    atom_coordinate = int(line[8:12])
                    amino_acid_coordinate = int(line[23:26])
                    atom_amino_acid_coordinates.append((atom_coordinate, amino_acid_coordinate))

        prev_aa_coord = None
        potential_gaps = []

        for atom_coord, aa_coord in atom_amino_acid_coordinates:
            if prev_aa_coord is not None and aa_coord != prev_aa_coord:
                if aa_coord - prev_aa_coord > 1:
                    # Found a potential gap
                    potential_gaps.append((prev_aa_coord + 1, aa_coord - 1))
            prev_aa_coord = aa_coord


        return potential_gaps

    def get_resolution(self):
        resolution = None
        with open(self.file_path, 'r') as file:
            # print(self.file_path)
            for line in file:
                if line.startswith('REMARK   2 RESOLUTION'):
                    try:
                        resolution = float(line[25:30])
                    except ValueError:
                        resolution = None  # P04439_6mpp.pdb
                    finally:
                        break
        return resolution


class AlphaFoldReader:
    def __init__(self, file_path):
        self.file_path = file_path

    def get_scores(self):
        scores = []

        prev_aa_coordinate = 0
        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    aa_coordinate = int(line[23:26])
                    score = float(line[61:66])
                    if aa_coordinate > prev_aa_coordinate:
                        scores.append(score)
                        prev_aa_coordinate = aa_coordinate

        return scores
