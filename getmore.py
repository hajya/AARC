from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO


E_VAL = 0.0001
titles = {}
theMAINlist = []
prots = {}
outfile = 'out.phy'
formata = 'phylip'

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
	b = NCBIXML.read(result_handle)
	return b

def getHomoLs(bb):
	handler = open("toAlign.fasta",'w')
	unqList = {}
	specCount = 1
	for b in bb.alignments:
		for hsp in b.hsps:
			if hsp.expect < E_VAL:
				c = 1
				start = 0
				finish = 0
				count = 0
				endTit = 0
				for m in b.title:
					if(m == '|'):
						if(count == 1):
							endTit = c-1
							count = 2

						if(count == 0):
							count = 1
				
					if(m == '['):
						start = c
					if(m == ']'):
						finish= c-1
						this = b.title[start:finish]
						ida = b.title[0:endTit]
						if(this not in titles 
							and this != 'synthetic construct' 
							and b.title[33:42] != 'PREDICTED'):
							titles[this] = this
							unqList[b.title[0:10]] = this
							unqList[b.title[0:10]+'COUNT'] = specCount 
							unqList[b.title[0:10]+'descr'] = b.title
							handler.write('>')
							handler.write(str(b.title))
							handler.write('\n')
							handler.write(str(hsp.sbjct))
							handler.write('\n\n')
						#	print(this)O
						#	print(hsp.sbjct)
							specCount += 1		
							break
					c += 1
	return unqList
def align():
	
	cline = ClustalwCommandline("clustalw", infile='toAlign.fasta',outfile=outfile,align='true',output='PHYLIP')
	cline()
	alignment = AlignIO.parse(outfile,formata).next()
	return alignment
def makeRec():
	thelist =[]

if __name__ == "__main__":
	#format the seq from the file
	record = SeqIO.read(open("hemo.fasta"), format="fasta")
	#makeXML(record)
	blast_records = parseXML()
	unqList = getHomoLs(blast_records)
	a = align(unqList)
