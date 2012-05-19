from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
import re
import sys

E_VAL_THRESHOLD = 0.0001
titles = {}
theMAINlist = []
prots = {}
outfile = 'out.phy'
formata = 'phylip'

def make_uid(x):
    def inc():
        return x + 1
    return inc

UID = make_uid(-1)

def makeXML(rec):
    #call to blast
    print(rec)
    result_handle = NCBIWWW.qblast("blastp", "nr", rec.format("fasta"))
    #put it into an xml file
    save_file = open("my_blast.xml", "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()


def parseXML():
    result_handle = open("my_blast.xml")
    b = NCBIXML.read(result_handle)  #read() is meant for a single query.  Should use parse()
    return b

def getHomologues(blast_records, searchQuery):
    '''Takes in a blast NCBIXML read and turns it into a
    homologue, returns a dictionary containing the following structure:
    {'searchQuery': <the function argument searchQuery>,
     'protiens':{ 'species-1': [{'id':'species-1', 'seq':<protien-seq-1>, 'desc':<blast-desc-1>, 'uid': 1},
                                {'id':'species-1', 'seq':<protien-seq-2>, 'desc':<blast-desc-2>, 'uid': 2},
                               ],
                  'species-2': [{'id':'species-2', 'seq':<protien-seq-3>, 'desc':<blast-desc-3>, 'uid': 3}
                               ],
                  . . .
                }
    }'''

    global E_VAL_THRESHOLD
    homologue = {}
    homologue['searchQuery'] = searchQuery
    homologue['protiens'] = {}
    protiens = homologue['protiens']
    for alignment in blast_records.alignments:
        if len(alignment.hsps) != 1: #ERROR CONDITION
            print >> sys.stderr, "ERROR: unexpected blast results" 
            continue
        if alignment.hsps[0].expect > E_VAL_THRESHOLD:
            continue
        species = re.findall("\[.*?\]", alignment.title)
        species = set(species) # get only unique elements
        if len(species) == 1:
            species_name = list(species)[0]
            if species_name not in protiens:
                protiens[species_name] = []
            species_protien_list = protiens[species_name]
            temp = {}
            temp['uid'] = "FOO" #UID() 
            temp['id'] = list(species)[0]
            temp['description'] = "TEMP"#alignment.title
            temp['score'] = alignment.hsps[0].score
            temp['seq'] = "TEMP" #alignment.hsps[0].sbjct
            species_protien_list.append(temp)
    return homologue

def align(homologue):
    '''Takes in a homologue from getHomologues() and
    aligns all of the sequences that it contains'''
    cline = ClustalwCommandline("clustalw", infile='toAlign.fasta',outfile=outfile,align='true',output='PHYLIP')
    cline()
    alignment = AlignIO.parse(outfile,formata).next()
    return alignment

def makeRec():
    thelist =[]

def alignments_to_json(filename, fileFormat):          #Experimental.  filename="out.phy", format="phylip"
    '''Opens a file with one alignment and converts it to a list that can
    be exported to a json string'''
    alignments = AlignIO.read(open(filename), fileFormat)
    homologue = {}
    homologue["protienID"] = alignments[0].description
    homologue["searchQuery"] = str(alignments[0].seq)
    protiens = []
    for element in alignments:
        temp = {}
        temp['id'] = element.id
        temp['description'] = element.description
        temp['seq'] = str(elemenet.seq)
        protiens.append(temp)
    homologue["protiens"] = protiens
    return homologue

if __name__ == "__main__":
    #format the seq from the file
    record = SeqIO.read(open("hemo.fasta"), format="fasta")
    #makeXML(record)
    blast_records = parseXML()
    unqList = getHomoLs(blast_records)
    print("-------------------------")
    print(unqList)
    print("-------------------------")
    #a = align(unqList)
    a = align()

    '''Desired Output:  A json string of the following format...
    A list of homologues.  Homologues are dictionaries

    [ { SearchQuery: " ... "                                             #Homologue number 1
        protiens: [ { "id":   unique-species-id
                      "seq":  aligned-sequence-of-uniform-length
                      "desc": " ... "
                    },
                    additional dictionaries for each unique-species-id...
                  ]
      },
      { ... }, #Homologue number 2...
      ...
    ]
    '''
