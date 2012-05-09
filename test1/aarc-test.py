from Bio import AlignIO

#Maximum (number of letter variations in a col)/(total # of seq)
MAX_VARIATION_RATIO = .3
LETTER_UNIQUENESS_WEIGHT = .3 #(Number of occurences of letter)/(Total # of letters)

def is_mutation(column, letter, colHeight):
    '''Checks to see if a given letter is a mutation within a column'''
    if ((float(column[letter]['count'])/float(colHeight)) < LETTER_UNIQUENESS_WEIGHT) and len(column) > 1:
        return True
    return False

def letter_contains_record(column, letter, record):
    '''Checks to see if a given letter in a column contains is from a given record (species)
    If so it returns the record for the corresponding species in the column'''
    #For every record in the letter see if record is the same as
    #as the one in the letter
    for recordB in column[letter]['records']:
        if column[letter]['records'][recordB].id == record.id:
            return column[letter]['records'][recordB]
    return False

def seq_column_combine(seqA, seqApos, column, colBpos, colHeight):
    '''Checks to see if the sequence given by seqA, with amino acid at position
    seqApos has a mutation in column, with column position (index) colBpos for a
    given number of variations (colHeight)'''
    tuples = []
    for letter in column:
        #If the given letter is a mutation
        if is_mutation(column, letter, colHeight):
            #Check the records of that letter
            seqB = letter_contains_record(column, letter, seqA)
            if seqB:
                ret = {}
                ret['colApos'] = seqApos # basically X cord for A
                ret['recordA'] = seqA.id.replace("|","")# Which read it came from (Y cord for A)
                ret['colBpos'] = colBpos # X cord for B
                ret['recordB'] = seqB.id.replace("|","") # Y cord for B (which read it came from)
                tuples.append(ret)
    #Quick checksum to make sure that this only returns one pair (optional)
    if len(tuples) > 1:
        print "ERROR: multiple tuples returned for a single amino acid pairing"
    return tuples

def seq_columns_combine(seq, seqPos, columns, colHeight):
    '''Takes in a given sequence at position seqPos (aka an amino acid) and finds all the pairs
    in columns where there is a mutation in a column in the same variation'''
    tuples = []
    col_pos = 0;
    for column in columns:
        tuples += seq_column_combine(seq, seqPos, column, col_pos, colHeight)
        col_pos += 1
    return tuples

def seqs_columns_combine(seqs, seqsPos, columns, colHeight):
    '''Takes in a series of sequences (aka multiple mutation in a given column in an alignment)
    and creates a list of tuples where there was also a mutation in columns in the same variation'''
    tuples = []
    for seq in seqs:
        tuples += seq_columns_combine(seqs[seq], seqsPos, columns, colHeight)
    return tuples

def column_columns_combine(colA, colApos, colsB, colHeight):
    '''takes in a given column and finds all mutations. This is then compared to columns (another protien)
    and all of the possible pairs of corresponding mutations are returned'''
    tuples = []
    for letter in colA:
        if is_mutation(colA, letter, colHeight):
            tuples += seqs_columns_combine(colA[letter]['records'], colApos,colsB, colHeight)
    return tuples


def columns_columns_combine(A, B):
    '''Takes in two homologues and finds all of the mutations in each. Then it takes all of the
    mutations in one variation (protien sequence) and returns each mutation paired with all of the
    mutations in the same variation (usually species) in the other protien.'''
    (colsA, colAmap) = A.get_possible_cols()
    (colsB, colBmap) = B.get_possible_cols()
    tuples = []
    col_pos = 0
    # Each item in colsA represents a set of mutations
    for mutations in colsA:
        tuples += column_columns_combine(mutations, col_pos, colsB, A.col_height)
        col_pos += 1
    #Correct the column mapping since some columns were filtered out
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
#   columns_columns_combine(A,B)
print len(columns_columns_combine(A,B))
for element in columns_columns_combine(A,B):
    print element
