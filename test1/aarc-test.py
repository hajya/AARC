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
                    ret['recordA'] = seqA.id.replace("|","")# Which read it came from (Y cord for A)
                    ret['colBpos'] = colBpos # X cord for B
                    ret['recordB'] = column[letter]['records'][seqB].id.replace("|","") # Y cord for B (which read it came from)
                    tuples.append(ret)
    return tuples

def seq_columns_check(seq, seqPos, columns):
    tuples = []
    col_pos = 0;
    for column in columns:
        #TODO: fix this, col_pos does not correspond to real column pos
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
        if (float(colA[letter]['count'])/float(colAlen) < LETTER_UNIQUENESS_WEIGHT) and (len(colA) > 1):
            tuples += seqs_columns_check(colA[letter]['records'], colApos,colsB)
    return tuples


def columns_columns_check(A, B):
    (colsA, colAmap) = A.get_possible_cols()
    (colsB, colBmap) = B.get_possible_cols()
    tuples = []
    col_pos = 0
    # Each item in colsA represents a set of mutations
    for mutations in colsA:
        col_length = 0
        for letter in mutations:
            col_length += mutations[letter]['count']
        tuples += column_columns_check(mutations, col_pos, col_length, colsB)
        col_pos += 1
    for element in tuples:
        element['colApos'] = colAmap[element['colApos']]
        element['colBpos'] = colBmap[element['colBpos']]
    return tuples

def filter_column(col, col_height):
    if "-" in col: #Filter out columns where there is an insertion in only a few homologues
        if float(col["-"]["count"])/float(col_height) < LETTER_UNIQUENESS_WEIGHT:
            return False
    if not (float(len(col))/float(col_height) < MAX_VARIATION_RATIO and len(col) > 1):
        return False 
    return True
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
        '''Filters out columns/amino-acids that have either too much variation or no variation
        Returns a list of the columns along with a list containing the map of the columns to
        their original positions'''
        cols = []
        col_len = len(self.columns)
        i = 0;
        column_map = []
        for column in self.columns:
            if filter_column(column,col_len):
                cols.append(column)
                column_map.append(i) # column_map[i] = original position of column in colum.self
            i += 1
        return (cols, column_map)

    
A = Homologue('rhoa-aligned-no-predicted.fa', 'fasta')
B = Homologue('rock-aligned-no-predicted.fa', 'fasta')
print "Number of columns in A: ", len(A.columns)
print "Number of columns in B: ", len(B.columns)
print "Number of Possible columns in A: ", len(A.get_possible_cols()[1])
print "Number of Possible columns in B: ", len(B.get_possible_cols()[1])
print "Total Number of possible combinations: ", len(B.get_possible_cols()[1]) * len(A.get_possible_cols()[1])
#for i in range(0,100):
#   columns_columns_check(A,B)
print len(columns_columns_check(A,B))
#for element in columns_columns_check(A,B):
#    print element
