#!/usr/bin/env python3
#when running program enter program name then dictionary file (flyFull.out mapping file from HW question#1) then file we want to look at (SSR...)


import sys

if len(sys.argv) < 2:
    print("need to enter map file < SRR file after program name")
    quit()
if len(sys.argv) < 3:
    print("enter S or U after SRR file name. S=skip lines with no match, U=print unknown")
    quit()
#making a dictionary of uniprot and flybase ids. file being used is flyFull.out from exercise 1

mapping_dict = {}
for line in open(sys.argv[1]):
    fields = line.strip().split()
    FlyBaseID_key = fields[0]
    UniprotID_value = fields[1]
    mapping_dict[FlyBaseID_key] = UniprotID_value

#using flybase ID from SRR file to search dictionary. if match found, print SRR line with newly found uniprot id


for i, line in enumerate(sys.stdin):
    if i == 0:
        continue
    fields = line.rstrip("\r\n").split()
    FlyBaseID = fields[8]
    if FlyBaseID in mapping_dict:
        UniProtID = mapping_dict[FlyBaseID]
        print(line + "  " + UniProtID)
#for Unknown to be filled in for uniprotids not found enter U in command line after SRR file. if you want line to be skipped enter S
        if sys.argv[2] == "U":
            print(line + "  " + "Unknown")
        if sys.argv[2] == "S":
            continue
        