from aarc import *
import json


A = alignments_to_json('data/rhoa-aligned-no-predicted.fa', 'fasta')
B = alignments_to_json('data/rock-aligned-no-predicted.fa', 'fasta')

print json.dumps(A)
print json.dumps(B)
