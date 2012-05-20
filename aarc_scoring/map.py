#!/usr/bin/python
from aarc import *
import json
import sys

for line in sys.stdin:
    #print >> sys.stderr, "INPUT"
    #print >> sys.stderr, line
    temp = json.loads(line)
    combinations = []
    #print >> sys.stderr, "END INPUT"
    for i,homologue in enumerate(temp):
        for element in temp[i+1:]:
            combinations.append([homologue, element])
    try:
        #foo = combinations[0]
        #for i in range(1,100000):
        print str(i) + " " + json.dumps(foo)
    	    #print >> sys.stderr, "MAP:Generated output"
    except:
        print >> sys.stderr, "MAP:Didn't recive any input"
        pass
