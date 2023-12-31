import numpy as np


class Reader:
    """"
    PDB file reader
    """
    def __init__(self, file_path, params):
        """
        Initialize reader and align valid sections
        :param file_path: folder path saving reading files
        :param params: PDBs' identifiers, valid regions and the differences
        """
        self.file_path = file_path
        self.params = params
        self.resolution = self.get_resolution()
        self.valid_sections = {}  # Eg. {'A':{'region':[[2, 292]], 'transform':0}, 'B':{'region':[[2, 292]], 'transform':0}}
        for section in self.params:
            chain = section['chain']
            region = [section['unp_start'], section['unp_end']]
            transform = section['residue_start'] - section['unp_start']
            if chain not in self.valid_sections.keys():
                self.valid_sections[chain] = {'region': [region], 'transform': transform}
            else:
                self.valid_sections[chain]['region'].append(region)

    def get_valid_sections(self):
        """
        :return: valid sections
        """
        return self.valid_sections

    def get_potential_gaps(self, sequence):
        prev_aa_coord = 0  # the previous amino acid coordinate
        gaps = {}
        with open(self.file_path, 'r') as file:
            for line in file:  # exam line by line
                if line.startswith('ATOM'):
                    chain = line[21]
                    aa_coordinate = int(line[22:26])
                    if chain in self.valid_sections.keys():  # if valid chain
                        is_in_region = False
                        regions = self.valid_sections[chain]['region']
                        transform = self.valid_sections[chain]['transform']
                        for region in regions:
                            if region[0] + transform <= aa_coordinate <= region[1] + transform:
                                is_in_region = True
                        if is_in_region:  # if in valid region
                            prev_aa_coord = max(prev_aa_coord, region[0] + transform)
                            if prev_aa_coord is not None and aa_coordinate - prev_aa_coord > 1:
                                if chain not in gaps.keys():
                                    gaps[chain] = []
                                gaps[chain].append([prev_aa_coord + 1 - transform, aa_coordinate - 1 - transform])
                                
                            prev_aa_coord = aa_coordinate
        gap_array = np.zeros((len(self.valid_sections.keys()), len(sequence)))

        # align pdb and uniprot sequence
        try:
            for row in range(len(self.valid_sections.keys())):
                valid_regions = self.valid_sections[list(self.valid_sections.keys())[row]]['region']
                for v_region in valid_regions:
                    for col in range(v_region[0] - 1, v_region[1]):
                        gap_array[row, col] = -self.resolution
            for row in range(len(gaps.keys())):
                for section in gaps[list(gaps.keys())[row]]:
                    for col in range(section[0] - 1, section[1]):
                        gap_array[row, col] = self.resolution
        # The IndexError will only happen if the source of PDB is wrong, so the whole PDB will be considered invalid
        except IndexError:
            return np.zeros((len(self.valid_sections.keys()), len(sequence)))
        return gap_array

    def get_resolution(self):
        """
        Get the resolution of the PDB
        :return: resolution
        """
        resolution = None
        with open(self.file_path, 'r') as file:
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
    """
    AlphaFold PDB file reader
    """
    def __init__(self, file_path):
        self.file_path = file_path

    def get_scores(self):
        scores = []
        prev_aa_coordinate = 0
        with open(self.file_path, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    aa_coordinate = int(line[22:26])
                    score = float(line[61:66])
                    if aa_coordinate > prev_aa_coordinate:
                        scores.append(score)
                        prev_aa_coordinate = aa_coordinate

        return scores
