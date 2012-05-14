from aarc import *
import json


A = alignments_to_json('rhoa-aligned-no-predicted.fa', 'fasta')
B = alignments_to_json('rock-aligned-no-predicted.fa', 'fasta')

print json.dumps(A)
print json.dumps(B)
