from Bio.Blast import NCBIWWW
from Bio import SeqIO

#format the seq from the file
record = SeqIO.read(open("hemo.fasta"), format="fasta")

#call to blast
result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
#put it into an xml file
save_file = open("my_blast.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()

