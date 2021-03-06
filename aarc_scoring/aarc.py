from Bio import AlignIO
from scoring import get_amino_acid_score
import pickle
import sys
import json
import math
#Maximum (number of letter variations in a col)/(total # of seq)
MAX_VARIATION_RATIO = .3
LETTER_UNIQUENESS_WEIGHT = .3 #(Number of occurences of letter)/(Total # of letters)
DELTA_MUTATION_SCORE_MINIMUM = .2 #Used to filter out points whose conservation scores are not close
MUTATION_SCORE_FILTER_THRESHOLD = 0.0 #Any point with a matrix amino acid score lower than this will be filtered
VARIATION_AMINO_ACID_SCORE_PENALTY = 2.0 #The higher this number the more that columns with lots of variation
                                         #are penalized in their score

MAX_ADD_MATRIX_SCORE = 6.0
MIN_ADD_MATRIX_SCORE = -3.0


def print_possible_tuple(tup):
    print_str = ""
    print_str += str(tup['colApos']) + ","
    print_str += str(tup['colBpos']) + ","
    #print_str = "["
    #print_str += '\'A-Pos\': ' + str(tup['colApos']) + ', '
    #print_str += '\'B-Pos\': ' + str(tup['colBpos']) + ', '
    #print_str += '\'deltaAScore\': ' + str(tup['deltaAScore']) + ', '
    #print_str += '\'deltaBScore\': ' + str(tup['deltaBScore']) + ', '
    #print_str += '\'deltaMutationScore\':' + str(tup['deltaMutationScore']) + ']'
    print print_str
    
def calc_mutation_score(column, letter):
    '''Uses a scoring matrix to calculate a score based
    upon what the amino acid changed from to'''
    #Select the most common letter in the column
    #As this will be what the letter mutated from
    
    #TODO: add some sort of modifier so that is there is more than one mutation
    #the score is worse
    
    max_letter_count = 0
    max_letter = ''
    for letter in column:
        if column[letter]['count'] > max_letter_count:
            max_letter_count = column[letter]['count']
            max_letter = letter
    #Score is matrix_score - 1/number_of_variations in column
    score_penalty = 1.0 - float(len(column))/VARIATION_AMINO_ACID_SCORE_PENALTY;
    return float(get_amino_acid_score(max_letter, letter)) - score_penalty;



def calc_delta_mutation_score(scoreA, scoreB, col_height):
    '''Calculates the similarity between two amino acid scores
    and assigns them a value. The theory is that the more similar
    the changes in the amino acids the more likely it is to be a
    covariance'''
    #Score is 1.0/(1 + abs(delta(a,b)))
    #.1 is to prevent diviide by zero
    #highest possible score is 10 and increases more rapidly
    #as the distance approaches 0

    #TODO: Fix this score, its kinda weird
    
    return float(1.0/(1.0 + abs(float(scoreA) - float(scoreB))))

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
        if column[letter]['records'][recordB]['id'] == record['id']:
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
                #Note I have tried to keep each tuple dumpable as a json
                #This leads to a lot of redundancy since each tuple
                #contains the sequence of both protiens and their desciptions and IDs
                ret = {}
                ret['colApos'] = seq_a_pos # basically X cord for A
                ret['seqA'] = seq_a
                ret['colBpos'] = col_b_pos # X cord for B
                ret['seqB'] = seq_b
                ret['seqID'] = seq_a['id'] #Id of the sequence
                #Score for the Amino acid change in protien A
                ret['deltaAScore'] = calc_mutation_score(A.filtered_columns[seq_a_pos], seq_a['seq'][seq_a_pos])
                #Score for the Amino acid change in protien B
                ret['deltaBScore'] = calc_mutation_score(b_col, b_letter);
                ret['deltaMutationScore'] = calc_delta_mutation_score(ret['deltaAScore'], ret['deltaBScore'],colHeight)
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
        #if float(col["-"]["count"])/float(col_height) < LETTER_UNIQUENESS_WEIGHT:
        #    return False
        return False
    if not (float(len(col))/float(col_height) < MAX_VARIATION_RATIO and len(col) > 1):
        return False 
    return True

def alignments_to_json(filename, fileFormat):
    '''Opens a file with one alignment and converts it to a list that can
    be exported to a json string'''
    alignments = AlignIO.read(open(filename), fileFormat)
    homologue = {}
    #homologue["proteinName"] = alignments[0].description
    homologue["searchQuery"] = str(alignments[0].seq)
    homologue['species'] = {}
    for element in alignments:
        homologue['species'][element.id] = []
        temp = {}
        temp['id'] = element.id
        temp['description'] = element.description
        temp['seq'] = str(element.seq)
        homologue['species'][element.id].append(temp)
    return homologue

def seq_to_dict(seq):
    '''Takes in a Biopython seq and returns
    a dict with the important bits to represent
    a protien'''
    ret_dict = {}
    ret_dict['id'] = seq.id
    ret_dict['description'] = seq.description
    ret_dict['seq'] = str(seq.seq)

def read_homologue_pair(pair):
    '''Reads in a pair of homologues
    stored in a dict'''
    A = jsonHomologue(pair[0])
    B = jsonHomologue(pair[1])

class jsonHomologue():
    '''Homologues represented by json lines'''
    def __init__(self, homologue):
        if 'proteinName' in homologue:
            self.proteinName = homologue["proteinName"] #Metadata on search
        else:
            self.proteinName = "N/A"
        self.searchQuery = homologue["searchQuery"] # Query used to find homologues
        self.alignments = homologue["species"] #this is the individual protien seqs
        #check to make sure alignments ar same length
        #Other sanity checks should go here
        first_key = False
        for species_name in self.alignments:
            if not first_key:
                first_key = species_name
            if len(self.alignments[species_name][0]['seq']) != len(self.alignments[first_key][0]['seq']):
                raise Exception("Length of alignments didn't match")
        self.columns = []
        columns = self.columns
        self.filtered_columns = []
        self.col_height = len(self.alignments)
        for i in range(0,len(self.alignments[first_key][0]['seq'])):
            columns.append({})
        for species in self.alignments:
            for i,letter in enumerate(self.alignments[species][0]['seq']):
                if not letter in columns[i]:
                    columns[i][letter] = {}
                    columns[i][letter]['count'] = 0
                    columns[i][letter]['records'] = {}
                columns[i][letter]['count'] += 1
                columns[i][letter]['records'][species] = self.alignments[species][0]

    def get_possible_cols(self):
        cols = []
        self.filtered_columns = []
        col_len = len(self.columns)
        column_map = []
        for i,column in enumerate(self.columns):
            if filter_column(column,col_len):
                cols.append(column)
                self.filtered_columns.append(column)
                column_map.append(i)
        return (cols, column_map)
    def __str__(self, column_len= 80):
        species_names = []
        seqs = []
        for species_name,val in self.alignments.items():
            species_names.append(species_name)
            seqs.append(val[0]['seq'])
        i = 0
        retval = self.proteinName + "\n"
        stuff_to_print = True
        while stuff_to_print:
            for species_name,seq in zip(species_names, seqs):
                retval += "%30s  " % species_name[:30]
                if column_len < (len(seq) - i):
                    retval += seq[i:i+column_len]
                else:
                    retval += seq[i:]
                    stuff_to_print = False
                retval += "\n"
            i = i + column_len + 1
            retval += "\n"
        return retval

def get_common_species_count(A,B):
    '''Takes in two homologues and returns the number of common species between them'''
    count = 0
    for key in A.alignments:
        if key in B.alignments:
            count += 1
    return count

    
def filter_results(results):
    filtered_results = []
    global DELTA_MUTATION_SCORE_MINIMUM
    global MAX_ADD_MATRIX_SCORE
    global MIN_ADD_MATRIX_SCORE
    for element in results:
        delta_score_add = float(element['deltaAScore']) + float(element['deltaBScore'])
        if float(element['deltaMutationScore']) < DELTA_MUTATION_SCORE_MINIMUM:
            continue
        elif (delta_score_add <= MAX_ADD_MATRIX_SCORE) and (delta_score_add >= MIN_ADD_MATRIX_SCORE):
            continue
        else:
            filtered_results.append(element)
    return filtered_results

def calc_raw_total_delta_mutation(results):
    i = 0.0
    for element in results:
        i += float(element['deltaMutationScore'])
    return i

def calc_total_delta_mutation(results, A, B):
    i = calc_raw_total_delta_mutation(results)
    i = i * float(1000)
    temp = float(min((len(A.columns), len(B.columns))))
    #return i / float(get_common_species_count(A,B) * len(A.columns) * len(B.columns))
    return (i / float(get_common_species_count(A,B) * temp * temp)) * 100.0
    #try:
    #    return math.log(i)
    #except:
    #    return 0.0
    

def calc_and_print_average(results, column):
    total = 0.0
    for element in results:
        total += float(element[column])
    total = total/float(len(results))
    print "Average " + column + ": " + str(total)
    return total

def print_averages(results):
    for element in ['deltaAScore', 'deltaBScore', 'deltaMutationScore']:
        calc_and_print_average(results, element)

def map():
    '''Reads from stdin a series of alignments and then creates
    a series '''
    pass


if __name__ == "__main__":

    A = alignments_to_json('data/rhoa-aligned.fa', 'fasta')
    B = alignments_to_json('data/rock-aligned.fa', 'fasta')
    A = jsonHomologue(A)
    B = jsonHomologue(B)

    print "Number of Species common in A and B: ", get_common_species_count(A,B)
    
    print "Number of columns in A: ", len(A.columns)
    print "Number of columns in B: ", len(B.columns)
    print "Number of Possible columns in A: ", len(A.get_possible_cols()[1])
    print "Number of Possible columns in B: ", len(B.get_possible_cols()[1])
    print "Total Number of possible combinations: ", len(B.get_possible_cols()[1]) * len(A.get_possible_cols()[1])
    print "\n"

    
    results = columns_columns_combine(A,B)
    print "Total Combinations found before filtration: ", len(results)
    print "Total Combinations deltaMutationScore: ", calc_total_delta_mutation(results, A, B)
    print "Averages Before filtration: "
    print_averages(results)
    print "\n"

    filtered_results = filter_results(results)
    print "Filtered Results combinations: ", len(filtered_results)
    print "Filtered Results deltaMutationScore: ", calc_total_delta_mutation(filtered_results, A, B)
    print "Averages After filtration: "
    print_averages(filtered_results)
    for element in results:
        print_possible_tuple(element)
    #for line in filtered_results:
    #    print json.dumps(line)
