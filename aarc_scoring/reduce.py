#!/usr/bin/python

import sys
from aarc import *
import json

print "%55s %55s %6s %6s %7s %15s %15s" % ("Protien A Name", "Protien B Name", "lenA", "lenB", "comspec", "raw_score", "norm_score")
for line in sys.stdin:
    try:
        line = line.lstrip("0123456789 ")
        homologues = json.loads(line)
        A = jsonHomologue(homologues[0])
        B = jsonHomologue(homologues[1])
        A_len = len(A.columns)
        B_len = len(B.columns)
        num_common_species = get_common_species_count(A,B)
        results = filter_results(columns_columns_combine(A,B))
        raw_score = calc_raw_total_delta_mutation(results)
        normalized_score = calc_total_delta_mutation(results, A,B)
        print "%55s %55s %6i %6i %7i %15f %15f" % (A.proteinName.rstrip(".fasta"), B.proteinName.rstrip(".fasta"), A_len, B_len, num_common_species, raw_score, normalized_score)
        #print A.proteinName.rstrip(".fasta") + "*" + B.proteinName.rstrip(".fasta") + ", A_len:" + str(A_len) + ", B_len: " + str(B_len) + ", common_species: " + str(num_common_species) + ", Raw Score: " + str(raw_score) + ", Normalized Score: " + str(normalized_score) 
    except:
        raise
        print >> sys.stderr, "ERROR RUNNING reduce"
#        print >> sys.stderr, line
