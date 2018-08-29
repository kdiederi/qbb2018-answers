#!/usr/bin/env python3

import sys
import fasta
#reader 1 is for subset.fa, reader 2 for droyak2.fa (the query)
reader1 = fasta.FASTAReader(open(sys.argv[2]))
reader2 = fasta.FASTAReader(open(sys.argv[1]))

kmers = {}
k = int(sys.argv[3])
target_position = []

for ident, sequence in reader1:
    for i in range(0, len(sequence)-k):
        kmer = sequence[i:i+k]
        if kmer not in kmers:
            kmers[kmer] = [i]
        else:
            kmers[kmer].append(i)



for ident, sequence in reader2:
    for i in range(0, len(sequence)-k):
        kmer = sequence[i:i+k]
        if kmer not in kmers:
            continue
        else:
            print("target sequence name: " + str(ident) + "   " + "target start: " + str(i) + "   " + "query start: " + str(kmers[kmer]) + "   " + "kmer: " + str(kmer))
