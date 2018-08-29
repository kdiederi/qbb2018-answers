#!/usr/bin/env python3

import sys


count = 0
for line in open(sys.argv[1]):
    if line.startswith("#!"):
        continue
    fields = line.rstrip("\r\n").split() 
    if "gene" in fields[2]:
        if "protein_coding" in line:
            count += 1
        else:
            continue
print(count)
