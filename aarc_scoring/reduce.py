#!/usr/bin/python

import sys
from aarc import *
import json


for line in sys.stdin:
    line = line.lstrip("0123456789 ")
    homologues = json.loads(line)
    A = jsonHomologue(homologues[0])
    B = jsonHomologue(homologues[1])
    results = filter_results(columns_columns_combine(A,B))
    print calc_total_delta_mutation(results)
