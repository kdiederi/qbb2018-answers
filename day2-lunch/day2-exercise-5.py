#!/usr/bin/env python3

import sys
import pandas as pd

if len(sys.argv) > 1:
    f = open( sys.argv[1] )
else:
    f = sys.stdin

count = 0
total = 0
   
for i, line in enumerate( f ):
    if line[0] == "@":
        continue
    count += 1
    fields = line.strip().split("\t")
    total += int(fields[4])
    mapfq_avg = float(total / count)
    
print(float(mapfq_avg))