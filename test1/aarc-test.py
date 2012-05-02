from Bio import AlignIO

#Maximum (number of letter variations in a col)/(total # of seq)
MAX_VARIATION_RATIO = .3


class Homologue:
    '''A class to store all homologues'''
    def __init__(self, filename, fileFormat):
        self.alignment = AlignIO.read(open(filename), fileFormat)
        self.columns = []
        columns = self.columns
        for i in range(0, alignmnet.get_alignment_length()):
            columns[i] = {}
        for record in alignment:
            for letter in record.seq:
                if not col in columns[i]:
                    columns[i][letter] = {}
                    columns[i][letter]['count'] = 0
                    columns[i][letter]['records'] = []
                columns[i][letter]['count'] += 1
                columns[i][letter]['records'] += record
    def get_possible_coalesence_sites(self):
        '''Filters out columns/amino-acids that have either too much variation or no variation'''
        cols = []
        for column in self.columns:
            if len(column)/len(columns) < MAX_VARIATION_RATIO and len(column) > 1:
                cols += column
        return cols
    def compare(self, other):
        colsA = self.get_possible_coalesence_sites()
        colsB = other.get_possible_coalesence_sites()
        for colA in colsA:
           for letter in colA:
               
            
