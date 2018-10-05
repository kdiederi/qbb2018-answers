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
    chromosome = []
    position =[]
    pvalue_log10 = []
    with open(fname) as f:
        for line in f:
            fields = line.strip().split()
            if "BP" in line:
                continue
            if "NA" in line:
                continue
            chrom = fields[0]
            pos = float(fields[2])
            pval = float(fields[8])
            pval_log10 = np.log10(pval)
            chromosome.append(chrom)
            position.append(pos)
            pvalue_log10.append(-pval_log10)

    plt.scatter(position,pvalue_log10)
    #instructions said higlight p values below 10E-5. -log(10E-5)=4
    plt.axhline(y=4, color='r', linestyle='--')
    plt.text(1700000,5,'threshold=\np<10E-5',fontsize= 'xx-small',horizontalalignment='center',verticalalignment='top')
    plt.xlabel("Position")
    plt.ylabel("-log(p)")
    plt.title(str(treatment) + " Manhattan Plot")
    plt.savefig(str(treatment) + ".png")
    plt.close()
