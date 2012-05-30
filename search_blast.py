#!/usr/bin/python

from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from search import *
import re
import sys
import os
import StringIO
import tempfile
import json
import subprocess

if len(sys.argv) > 1:
    hom_list = []
    for filename in sys.argv[1:]:
        try:
            try:
                sys.stderr.write("OPENING: " + filename + "\n")
                record = SeqIO.read(open(filename), format="fasta")
            except:
                sys.stderr.write("ERROR reading file: " + filename + "\n")
            #"create" the .xml file for NCBIXML
            try:
                result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"))
                xml_file = open(filename + ".xml", "w")
                xml_file.write(result_handle.read())
                result_handle.seek(0)
                with tempfile.NamedTemporaryFile() as temp_file:
                    tempstr = result_handle.read()
                    result_handle.close()
                    temp_file.write(tempstr)
                    temp_file.flush()
                    blast_records = NCBIXML.read(open(temp_file.name))
            except:
                sys.stderr.write("ERROR parsing results from blast for file: " + filename + "\n")
            hom = getHomologues(blast_records, record.id + "\n" + str(record.seq), filename)
            hom = align(hom)
            hom_list.append(hom)
        except Exception as undefined_error:
            raise
            
        
else:
    sys.stderr.write("Reading from STDIN not supported yet\n")
    exit(0)
    #read from stdin
print json.dumps(hom_list)
