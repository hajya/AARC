result_handle = open("my_blast.xml")

#parse the xml
from Bio.Blast import NCBIXML
blast_records = NCBIXML.read(result_handle)

E_VAL = 0.0001
titles = []
for b in blast_records.alignments:
	for hsp in b.hsps:
		if hsp.expect < E_VAL:
			titles.append(b.title)			
#for i in titles:
	c = 0
	for m in titles[0]:
		c += 1
		print(m)
