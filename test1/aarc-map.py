from aarc import *
import json


A = alignments_to_json('rhoa-aligned-no-predicted.fa', 'fasta')
B = alignments_to_json('rock-aligned-no-predicted.fa', 'fasta')

foo = [json.loads(A),json.loads(B)]
for i in range(0,1000):
    print json.dumps(foo)
