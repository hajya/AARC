result_handle = open("my_blast.xml")

#parse the xml
from Bio.Blast import NCBIXML
blast_records = NCBIXML.read(result_handle)

E_VAL = 0.0001
i = 0
titles = []
for b in blast_records.alignments:
	for hsp in b.hsps:
		if hsp.expect < E_VAL:
			titles.append(b.title)			
for i in titles:
	print(i)
	print('\n')
