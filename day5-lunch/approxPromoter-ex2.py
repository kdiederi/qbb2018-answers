#!/usr/bin/env python3

"""
Usage: ./approxPromoter.py <t_data.ctab>

create a new file that approximates promoter region (p_start, p_end).
output is a tab separated file .bed
contains t_name, chrom, p_start, p_end
"""



import sys


for i, line in enumerate(open(sys.argv[1])):
    if i == 0:
        continue
    fields = line.strip().split()
    chrom = fields[1]
    t_name = fields[5]
    start = int(fields[3])
    end = int(fields[4])
    
    if fields[2] == "+":
        p_start = start-500
        p_end = start+500
    if fields[2] == "-":
        p_start = end-500
        p_end = end+500
    if p_start < 0 or p_end < 0:
        continue
    print(chrom+"\t"+str(p_start)+"\t"+str(p_end)+"\t"+t_name)
    # print ("\t".join([chrom, t_name, sdf]))
    

