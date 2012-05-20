#!/usr/bin/python2.7

from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import re
import sys
import os
import StringIO
import tempfile
import json

E_VAL_THRESHOLD = 0.0001
titles = {}
theMAINlist = []
prots = {}
#outfile = 'out.phy'
formata = 'phylip'

'''Takes 0 or more arguments, which are names of fasta files to read in.  If 0 arguments, assumes files
will be given from stdin'''

class unique_identifier:
    def __init__(self, x):
        self.count = x
    def get_uid(self):
        self.count += 1
        return self.count
UID = unique_identifier(-1)

def parseXML():
    result_handle = open("my_blast.xml")
    b = NCBIXML.read(result_handle)  #read() is meant for a single query.  Should use parse()
    return b

def getHomologues(blast_records, searchQuery):
    '''Takes in a blast NCBIXML read and turns it into a
    homologue, returns a dictionary containing the following structure:
    {'searchQuery': <the function argument searchQuery>,
     'species':{ 'species-1': [{'id':'species-1', 'seq':<protien-seq-1>, 'desc':<blast-desc-1>, 'uid': 1},
                                {'id':'species-1', 'seq':<protien-seq-2>, 'desc':<blast-desc-2>, 'uid': 2},
                               ],
                  'species-2': [{'id':'species-2', 'seq':<protien-seq-3>, 'desc':<blast-desc-3>, 'uid': 3}
                               ],
                  . . .
                }
    }'''

    global E_VAL_THRESHOLD
    global UID
    homologue = {}
    homologue['searchQuery'] = searchQuery
    homologue['species'] = {}
    species = homologue['species']
    for alignment in blast_records.alignments:
        if len(alignment.hsps) != 1: #ERROR CONDITION
            print >> sys.stderr, "ERROR: unexpected blast results" 
            continue
        if alignment.hsps[0].expect > E_VAL_THRESHOLD:
            continue
        species_results = re.findall("\[.*?\]", alignment.title)
        species_set = set(species_results) # get only unique elements
        if len(species_set) == 1:
            species_name = list(species_set)[0]
            if species_name not in species:
                species[species_name] = []
            species_protien_list = species[species_name]
            temp = {}
            temp['uid'] = UID.get_uid() 
            temp['id'] = species_name
            temp['description'] = alignment.title
            temp['score'] = alignment.hsps[0].score
            temp['seq'] = alignment.hsps[0].sbjct
            species_protien_list.append(temp)
    return homologue

def align(hom):
    '''Takes in a homologue from getHomologues() and
    aligns all of the sequences that it contains'''
    with tempfile.NamedTemporaryFile() as temp_file:
        uid_map = {}
        for species in hom['species']:
            temp = hom['species'][species][0]
            temp_file.write('>' + str(temp['uid']) + '\n')
            temp_file.write(temp['seq'] + '\n')
            temp_file.write('\n')
            uid_map[temp['uid']] = species
        temp_file.flush()
        align_io_temp_file = StringIO.StringIO()
        cline = ClustalwCommandline("clustalw", infile=temp_file.name, align='true',output='PHYLIP')
        align_io_temp_file.write(cline)
        alignments = AlignIO.parse(align_io_temp_file, 'phylip')
    for alignment in alignments:
        #TODO: Don't throw out data here
        #get the first protien for the species
        temp = hom[uid_map[alignment.id]]['species'][0]
        #clear out all others since we only currently want one
        hom[uid_map[alignment.id]]['species'] = []
        #stick in the aligned sequence
        temp['seq'] = str(alignment.seq)
        #re append the protien
        hom[uid_map[alignment.id]]['species'].append(temp)
    seq_len = False
    for species in hom['species']:
        if not seq_len:
            seq_len = len(hom['species'][species][0]['seq'])
        if len(hom['species'][species][0]['seq']) != seq_len:
            raise Exception("ALIGNMENT ERROR")
    return hom

#if __name__ == "__main__":
if len(sys.argv) > 1:
    hom_list = []
    for filename in sys.argv[1:]:
        #format the seq from the file
#        print "filename is"+filename
        record = SeqIO.read(open(filename), format="fasta")
        #"create" the .xml file for NCBIXML
        result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
        with tempfile.NamedTemporaryFile() as temp_file:
            tempstr = result_handle.read()
            result_handle.close()
            temp_file.write(tempstr)
            temp_file.flush()
            blast_records = NCBIXML.read(open(temp_file.name))
        #temp_file.write(result_handle.read())     #Write this to a file we will immediately throw away
#        result_handle.close()
        #blast_records = NCBIXML.read(temp_file)    #Can't avoid using this function
        hom = getHomologues(blast_records, "FILL_IN")
        hom = align(hom)
        hom_list.append(hom)
else:
    pass
    #read from stdin
print json.dumps(hom_list)
