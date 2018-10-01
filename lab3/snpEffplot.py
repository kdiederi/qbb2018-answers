#!/usr/bin/env python3

"""
parse snpEff output file,
create plots for:
a) read depth distribution across each variant
b) genotype quality distribution
c) allele frequency spectrum of identified varients
d) summary of predicted effect of each variant (barplot?)
 
./snpEffplot.py <snpEffoutput.vcf> <predicted_effects.csv>
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#parse through snpEff output file (which is tab delimited for first columns, then ; delimited)

read_depth = []
allele_frequency = []
genotype_quality = []
i=1

for line in open(sys.argv[1]):
    if line.startswith("#"):
        continue
    else:
        fields = line.strip().split()
        quality = fields[9].split(":")
        try:
            gq = quality[1]
        except IndexError:
            i+=1
        genotype_quality.append(float(gq))
        info = fields[7]
        data = info.strip(",").split(";")
        readdepth = data[7].split("=")
        read_depth.append(readdepth[1])
        allelefreq = data[3].split("=")
        allele_frequency.append(allelefreq[1])

read_depth_new = []
for line in read_depth:
    rd_new = line.split(",")[0]
    read_depth_new.append(float(rd_new))
predicted_effect_count = []
predicted_effect_name = []

allele_freq_new = []
for line in allele_frequency:
    afreq_new = line.split(",")[0]
    allele_freq_new.append(float(afreq_new))

file2 = open(sys.argv[2])
for line in file2:
    col = line.rstrip().split()
    name = col[0]
    predicted_effect_name.append(name)
    count = col[1]
    predicted_effect_count.append(int(count))


fig, axes = plt.subplots(nrows=2,ncols=2)
axes = axes.flatten()
plt.tight_layout(pad=4.5,w_pad=4.5,h_pad=4.5)
fig.set_size_inches(15,15)

#plot histogram of varient effect
axes[0].bar(predicted_effect_name, predicted_effect_count)
axes[0].set_xlabel('Predicted Effect')
axes[0].set_ylabel('Frequency')
axes[0].set_title("snpEff predicted variant effects")
#sets axis labels to the variant name, rotates label
axes[0].set_xticklabels(predicted_effect_name, rotation=90)

#plot read depth
axes[1].hist(np.log10(read_depth_new), bins=100)
axes[1].set_xlabel('Variant')
axes[1].set_ylabel('log10 Read Depth')
axes[1].set_title("Read depth distribution")


#plot allele freq
axes[2].hist(allele_freq_new, bins=50)
axes[2].set_xlabel('Variant')
axes[2].set_ylabel('Allele Frequency')
axes[2].set_title("Allele frequency spectrum")

#plot quality

axes[3].hist(genotype_quality, bins=50)
axes[3].set_xlabel('Variant')
axes[3].set_ylabel('Quality Score')
axes[3].set_title("Genotype quality distribution")

# plt.tight_layout()
fig.savefig("multipanelPlot.png")
plt.close(fig)
