#!/usr/bin/env python3

"""
calculate  the number of contigs, minimum/maximum/average contig length, and N50
 
./contig_compute.py <contigs.fasta> > <output>

contigs.fasta generated from velvet or SPAdes
"""

import sys
import fasta
import numpy as np


contig_list = []
contig_length = []
sorted_contig_length = []

contig_count = 0
for ident, sequence in fasta.FASTAReader(open(sys.argv[1])):
    contig_count += 1
    contig_list.append(sequence)
# print(contig_list)

for i in contig_list:
    length = len(i)
    contig_length.append(length)

print("#contigs= " + str(contig_count))
print("min contig length= " + str(min(contig_length)))
print("max contig length= " + str(max(contig_length)))
print("avg contig length= " + str(np.mean(contig_length)))

#to calculate N50
#sorts list and re-stores it as the sorted one
contig_length.sort()
total_size = sum(contig_length)
fiftypercent = total_size / 2
# print("total size= " + str(total_size))
# print("50% = " + str(fiftypercent))
contig_sum = 0
for k in contig_length:
    contig_sum += int(k)
    if contig_sum < fiftypercent:
        continue
    else:
        print("N50= " + str(k))
        break


