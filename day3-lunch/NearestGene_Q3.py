#!/usr/bin/env python3

import sys

new_closest_proteinCoding_pos = 30000000000
new_closest_proteinCoding_gene = "x"
new_closest_noncoding_pos = 30000000000
new_closest_noncoding_gene = "x"


for line in open(sys.argv[1]):
    if line.startswith("#!"):
        continue
    fields = line.rstrip("\r\n").split() 
    if "3R" in fields[0] and "gene" in fields[2]:
        for i, val in enumerate(fields):
            if val == "gene_biotype":
                if fields[i+1] == '"protein_coding";':
        
                    my_dist = 0
                    find_pos = 21378950
                    gene_start = int(fields[3])
                    gene_end = int(fields[4])
        
        
                    if find_pos < gene_start:
                        my_dist = gene_start - find_pos
                    elif find_pos > gene_end:
                        my_dist = find_pos - gene_end
                    if my_dist < new_closest_proteinCoding_pos:
                        new_closest_proteinCoding_pos = my_dist
                        for i, val in enumerate(fields):
                            if val == "gene_name":
                                gene_name = fields[i+1]
                                new_closest_proteinCoding_gene = gene_name
        for i, val in enumerate(fields):
            if val == "gene_biotype":
                if fields[i+1] != '"protein_coding";':
        
                    my_dist = 0
                    find_pos = 21378950
                    gene_start = int(fields[3])
                    gene_end = int(fields[4])
        
        
                    if find_pos < gene_start:
                        my_dist = gene_start - find_pos
                    elif find_pos > gene_end:
                        my_dist = find_pos - gene_end
                    if my_dist < new_closest_noncoding_pos:
                        new_closest_noncoding_pos = my_dist
                        for i, val in enumerate(fields):
                            if val == "gene_name":
                                gene_name = fields[i+1]
                                new_closest_noncoding_gene = gene_name
        
print("Closest Protein Coding Gene and distance: " + " " + str(new_closest_proteinCoding_gene) + " " + str(new_closest_proteinCoding_pos))
print("Closest Non-Protein Coding Gene and distance: " + " " + str(new_closest_noncoding_gene) + " " + str(new_closest_noncoding_pos))
