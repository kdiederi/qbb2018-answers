#!/usr/bin/env python3

"""
make manhattan plots for all plink files
 
./manhattanplot.py <plink.*.qassoc>
"""

import sys
import matplotlib.pyplot as plt
import numpy as np


for fname in sys.argv[1:]:
    treatment = fname.strip("plink."".qassoc")
    data = {}
    with open(fname) as f:
        for line in f:
            fields = line.strip().split()
            if "BP" in line:
                continue
            if "NA" in line:
                continue
            chrom = fields[0]
            pos = int(fields[2])
            pval = float(fields[8])
            pval_log10 = -np.log10(pval)
            if chrom not in data:
                data[chrom] = {'positions':[], 'logpvals':[]}
            data[chrom]['positions'].append(pos)
            data[chrom]['logpvals'].append(pval_log10)

    fig,ax = plt.subplots(figsize=(20,5))

    colors = ['skyblue', 'sandybrown']
    highlights = ['steelblue', 'coral']

    offset = 0
    tick_pos = []
    tick_labels = []
    
    for i, chrom in enumerate(data.keys()):
        x = np.array(data[chrom]['positions'])
        y = np.array(data[chrom]['logpvals'])
    
        sig = (y > 5)
    
        ax.scatter(x[sig] + offset,y[sig], marker='.', c=highlights[i%2])
        ax.scatter(x[~sig] + offset,y[~sig], marker='.', c=colors[i%2])
    
        tick_labels.append(chrom)
        maxx = max(x)
        tick_pos.append(offset + maxx/2)
        offset += maxx
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels)
    ax.axhline(5, c='k', ls=':', label="Significance Cutoff")
    ax.legend()
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("-log10 P-value")
    ax.set_title("Manhattan Plot\n"+str(treatment))
    plt.savefig(str(treatment) + ".png")
    plt.close()
