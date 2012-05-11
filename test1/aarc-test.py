from Bio import AlignIO
from scoring import get_amino_acid_score

#Maximum (number of letter variations in a col)/(total # of seq)
MAX_VARIATION_RATIO = .3
LETTER_UNIQUENESS_WEIGHT = .3 #(Number of occurences of letter)/(Total # of letters)


def print_possible_tuple(tup):
    print_str = "["
    print_str += '\'A-Pos\': ' + str(tup['colApos']) + ', '
    print_str += '\'B-Pos\': ' + str(tup['colBpos']) + ', '
    print_str += '\'deltaAScore\': ' + str(tup['deltaAScore']) + ', '
    print_str += '\'deltaBScore\': ' + str(tup['deltaBScore']) + ', '
    print_str += '\'deltaMutationScore\':' + str(tup['deltaMutationScore']) + ']'
    print print_str
    
def calc_mutation_score(column, letter):
    '''Uses a scoring matrix to calculate a score based
    upon what the amino acid changed from to'''
    #Select the most common letter in the column
    #As this will be what the letter mutated from
    max_letter_count = 0
    max_letter = ''
    for letter in column:
        if column[letter]['count'] > max_letter_count:
            max_letter_count = column[letter]['count']
            max_letter = letter
    return get_amino_acid_score(max_letter, letter);

def calc_delta_mutation_score(scoreA, scoreB):
    '''Calculates the similarity between two amino acid scores
    and assigns them a value. The theory is that the more similar
    the changes in the amino acids the more likely it is to be a
    covariance'''
    #Score is 1/(.1 + abs(delta(a,b)))
    #.1 is to prevent diviide by zero
    #highest possible score is 10 and increases more rapidly
    #as the distance approaches 0
    return float(1.0/(0.1 + abs(float(scoreA) - float(scoreB))))

def calc_mutation_count_score(column, letter):
    '''Uses the number of mutations in a column to get
    a score for the variation in the column and how unique
    a given letter is'''
    #TODO: finish this not sure quite how to score this
    return 0

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

def seq_column_combine(seq_a, seq_a_pos, b_col, col_b_pos, A, B):
    '''Checks to see if the sequence given by seqA, with amino acid at position
    seqApos has a mutation in column, with column position (index) colBpos for a
    given number of variations (colHeight)'''
    colHeight = A.col_height
    tuples = []
    for b_letter in b_col:
        #If the given letter is a mutation
        if is_mutation(b_col, b_letter, colHeight):
            #Check the records of that letter
            #To see if the sequence A is in it
            seq_b = letter_contains_record(b_col, b_letter, seq_a)
            #If so we have a match
            if seq_b:
                ret = {}
                ret['colApos'] = seq_a_pos # basically X cord for A
                ret['seqA'] = seq_a # Sequence from protien A
                ret['colBpos'] = col_b_pos # X cord for B
                ret['seqB'] = seq_b # Sequence from protien B
                ret['seqID'] = seq_a.id #Id of the sequence
                #Score for the Amino acid change in protien A
                ret['deltaAScore'] = calc_mutation_score(A.filtered_columns[seq_a_pos], seq_a[seq_a_pos])
                #Score for the Amino acid change in protien B
                ret['deltaBScore'] = calc_mutation_score(b_col, b_letter);
                ret['deltaMutationScore'] = calc_delta_mutation_score(ret['deltaAScore'], ret['deltaBScore'])
                tuples.append(ret)
    #Quick checksum to make sure that this only returns one pair (optional)
    if len(tuples) > 1:
        print "ERROR: multiple tuples returned for a single amino acid pairing"
    return tuples

def seq_columns_combine(seq_a, seq_a_pos, A, B):
    '''Takes in a given sequence at position seqPos (aka an amino acid) and finds all the pairs
    in columns where there is a mutation in a column in the same variation'''
    colHeight = A.col_height
    tuples = []
    for b_col_pos, b_col in enumerate(B.filtered_columns):
        tuples += seq_column_combine(seq_a, seq_a_pos, b_col, b_col_pos, A, B)
    return tuples

def seqs_columns_combine(seqs_a, seqs_a_pos, A, B):
    '''Takes in a series of sequences (aka multiple mutation in a given column in an alignment)
    and creates a list of tuples where there was also a mutation in columns in the same variation'''
    colHeight = A.col_height
    tuples = []
    for seq_a in seqs_a:
        tuples += seq_columns_combine(seqs_a[seq_a], seqs_a_pos, A, B)
    return tuples

def column_columns_combine(col_a, col_a_pos, A, B):
    '''takes in a given column and finds all mutations. This is then compared to columns (another protien)
    and all of the possible pairs of corresponding mutations are returned'''
    colHeight = A.col_height
    tuples = []
    for letter in col_a:
        if is_mutation(col_a, letter, colHeight):
            tuples += seqs_columns_combine(col_a[letter]['records'], col_a_pos, A, B)
    return tuples


def columns_columns_combine(A, B):
    '''Takes in two homologues and finds all of the mutations in each. Then it takes all of the
    mutations in one variation (protien sequence) and returns each mutation paired with all of the
    mutations in the same variation (usually species) in the other protien.'''
    (cols_a, col_a_map) = A.get_possible_cols()
    (cols_b, col_b_map) = B.get_possible_cols()
    tuples = []
    # Each item in colsA represents a set of mutations
    for col_a_pos, col_a in enumerate(cols_a):
        tuples += column_columns_combine(col_a, col_a_pos, A, B)
    #Correct the column mapping since some columns were filtered out
    for element in tuples:
        element['colApos'] = col_a_map[element['colApos']]
        element['colBpos'] = col_b_map[element['colBpos']]
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
        self.filtered_columns = []
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
        their original positions and sets the variable self.filtered_cols to the columns left over
        from the fileter.'''
        cols = []
        self.filtered_columns = []
        col_len = len(self.columns)
        column_map = []
        for i, column in enumerate(self.columns):
            if filter_column(column,col_len):
                cols.append(column)
                self.filtered_columns.append(column)
                column_map.append(i) # column_map[i] = original position of column in colum.self
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
    print_possible_tuple(element)
