#!/usr/bin/python
from aarc import *
import json
import sys

temp = []
for line in sys.stdin:
    temp.append(json.loads(line))
combinations = []
for i,homologue in enumerate(temp):
    for element in temp[i+1:]:
        combinations.append([homologue, element])
try:
    foo = combinations[0]
    for i in range(0,100000):
        print str(i) + " " + json.dumps(foo)
except:
    pass
