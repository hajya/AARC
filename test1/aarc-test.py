from Bio import AlignIO

#Maximum (number of letter variations in a col)/(total # of seq)
MAX_VARIATION_RATIO = .3
LETTER_UNIQUENESS_WEIGHT = .3 #(Number of occurences of letter)/(Total # of letters)

def seq_column_check(seqA, seqApos, column, colBpos, colSize):
    '''Checks to see if a seq is a mutation in a column'''
    tuples = []
    for letter in column:
        if float(column[letter]['count'])/float(colSize) < LETTER_UNIQUENESS_WEIGHT and len(column) > 1:
            for seqB in column[letter]['records']:
                if column[letter]['records'][seqB].id == seqA.id:
                    ret = {}
                    ret['colApos'] = seqApos # basically X cord for A
                    ret['recordA'] = seqA.id # Which read it came from (Y cord for A)
                    ret['colBpos'] = colBpos # X cord for B
                    ret['recordB'] = column[letter]['records'][seqB].id # Y cord for B (which read it came from)
                    tuples += ret
    return tuples

def seq_columns_check(seq, seqPos, columns):
    tuples = []
    col_pos = 0;
    for column in columns:
        tuples += seq_column_check(seq, seqPos, column, col_pos, len(columns))
        col_pos += 1
    return tuples

def seqs_columns_check(seqs, seqsPos, columns):
    tuples = []
    for seq in seqs:
        tuples += seq_columns_check(seqs[seq], seqsPos, columns)
    return tuples

def column_columns_check(colA, colApos, colAlen, colsB):
    tuples = []
    for letter in colA:
        if float(colA[letter]['count'])/float(colAlen) < LETTER_UNIQUENESS_WEIGHT and len(colA) > 1:
            tuples += seqs_columns_check(colA[letter]['records'], colApos,colsB)
    return tuples


def columns_columns_check(A, B):
    colsA = A.get_possible_cols()
    colsB = B.get_possible_cols()
    tuples = []
    col_pos = 0
    # Each item in colsA represents a set of mutations
    for mutations in colsA:
        col_length = 0
        for letter in mutations:
            col_length += mutations[letter]['count']
        tuples += column_columns_check(mutations, col_pos, col_length, colsB)
        col_pos += 1
    return tuples
            
class Homologue:
    '''A class to store all homologues'''
    def __init__(self, filename, fileFormat):
        self.alignment = AlignIO.read(open(filename), fileFormat)
        self.columns = []
        columns = self.columns
        self.col_height = len(self.alignment._records)
        for i in range(0, self.alignment.get_alignment_length()):
            columns.append({}) 
        for record in self.alignment:
            i = 0
            for letter in record.seq:
                if not letter in columns[i]:
                    columns[i][letter] = {}
                    columns[i][letter]['count'] = 0
                    columns[i][letter]['records'] = {}
                columns[i][letter]['count'] += 1
                columns[i][letter]['records'][record.id] = record 
                i += 1
    def get_possible_cols(self):
        '''Filters out columns/amino-acids that have either too much variation or no variation'''
        cols = []
        col_len = len(self.columns)
        for column in self.columns:
            if float(len(column))/float(self.col_height) < MAX_VARIATION_RATIO and len(column) > 1:
                cols.append(column)
        return cols



A = Homologue('rhoa-aligned.fa', 'fasta')
B = Homologue('rock-aligned.fa', 'fasta')
print len(A.columns)
print len(B.columns)
print len(A.columns) * len(B.columns)
print len(B.get_possible_cols()) * len(A.get_possible_cols())
#for i in range(0,100):
#    columns_columns_check(A,B)
print len(columns_columns_check(A,B))
