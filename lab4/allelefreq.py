#!/usr/bin/env python3

"""
Use allele frequency from vcf file and make histogram
 
./allelefreq.py <BYxRM_segs_saccer3.bam.simplified.vcf>
"""

import sys
import matplotlib.pyplot as plt

vcf_file = open(sys.argv[1])

allele_freq = []
for line in vcf_file:
    if line.startswith("#"):
        continue
    fields = line.split()
    info = fields[7]
    for id_val in info.split(";"):
        id, val = id_val.split("=")
        if id == "AF":
            allele_freq.append(float(val.split(",")[0]))

fig1, ax1 = plt.subplots()


ax1.hist(allele_freq, bins=60)
ax1.set_xlabel('Variant')
ax1.set_ylabel('Count')
ax1.set_title("Allele frequencies")
fig1.savefig("allelefreqPlot.png")
plt.close(fig1)