#!/usr/bin/env python3
import sys

gene_type_count = {}

for line in open(sys.argv[1]):
    if line.startswith("#!"):
        continue
    fields = line.rstrip("\r\n").split() 
    if "gene" in fields[2]:
        for i, val in enumerate(fields):
            if val == "gene_biotype":
                gene_name = fields[i+1]
                if gene_name not in gene_type_count:
                    gene_type_count[gene_name] = 1
                else:
                    gene_type_count[gene_name] += 1

for gene_type, count in gene_type_count.items():
    print(gene_type, count)