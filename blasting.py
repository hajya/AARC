result_handle = open("my_blast.xml")
from Bio import SeqIO
from Bio.Blast import NCBIXML
blast_records = NCBIXML.read(result_handle)
specCount = 1
E_VAL = 0.0001
titles = {}
theMAINlist = []
unqList = {}
thelist =[]
record = SeqIO.read(open("hemo.fasta"), format="fasta")
#for recs in record:
prots = {}
handler = open("toAlign.fasta",'w')
prots['proteinID'] = record.description
prots['SearchQuery'] = str(record.seq)
for b in blast_records.alignments:
    for hsp in b.hsps:
        if hsp.expect < E_VAL:
            c = 1
            start = 0;
            finish = 0
            for m in b.title:
                if(m == '['):
                    start = c
                if(m == ']'):
                    finish= c-1
                    this = b.title[start:finish]
                    
                    if(this not in titles 
                        and this != 'synthetic construct' 
                        and b.title[33:42] != 'PREDICTED'):
                        titles[this] = this
                        unqList['id'] = specCount 
                        unqList['description'] = b.title
                        handler.write('>')
                        handler.write(str(b.title))
                        handler.write('\n')
                        handler.write(str(hsp.sbjct))
                        handler.write('\n\n')
                        
                    #   print(this)O
                    #   print(hsp.sbjct)
                        specCount += 1      
                        thelist.append(unqList)
                        
                        break
                c += 1
prots['Proteins'] = thelist
theMAINlist.append(prots)
thelist = []


'''for m in theMAINlist:
    print(m['proteinID'])
    print(m['SearchQuery'])
    for prots in m['Proteins']:
        print(prots['id'])
        print(prots['description'])
        print(prots['seq'])'''
