#!/usr/bin/env python3
#we only want entries with DROME. there are many entries with empty columns, we are ignoring those.  

import sys

for line in open(sys.argv[1]):

    if "DROME" in line:
        fields = line.strip("\r\n").split()
        if len(fields) == 4:
            print(fields[3], fields[2])
