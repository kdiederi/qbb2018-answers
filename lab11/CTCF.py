#!/usr/bin/env python2

"""
CTCF interactions
 
./CTCF.py <Chr17_GSM2418860_WT_CTCF_peaks.txt>
"""

import sys
import hifive
import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1], sep="\t", index_col="chr", header=0)
df['CTCFmidpt'] = (df["start"] + df["end"])/2

#get binds for CTCF peaks
CTCF_list = df['CTCFmidpt'].values.tolist()
CTCF_subtract = np.subtract(CTCF_list, 15000000)

#divide by 10000 for bin size
CTCF_mid = np.divide(CTCF_subtract, 10000).astype(np.int32)
CTCF_uniq = np.unique(CTCF_mid)



hic = hifive.HiC('normalizing.hcp')
data = hic.cis_heatmap(chrom='chr17', start=15000000, stop=17500000, binsize=10000, datatype='fend', arraytype='full')

data[:, :, 1] *=np.sum(data[:, :, 0]) / np.sum(data[:, :, 1])
where = np.where(data[:, :, 1] > 0)
data[where[0], where[1], 0] /= data[where[0], where[1], 1]
data = data[:, :, 0]

# for all pairwise combintations of bins, which have both CTCF and if so is there interaction above 1. which cells in matrix have value greater than one?
for i in range(len(CTCF_uniq)):
    jcount = 0
    for j in range((i+1),len(CTCF_uniq)):
        row = CTCF_uniq[i]
        column = CTCF_uniq[j]
        enrichment = data[row,column]
        if enrichment >= 1:
        # print bin coordinates and enrichment value
            print("bin coordinate= " + str([i,j]) + "   " + "enrichment= " + str(enrichment))
        jcount = jcount + 1
        if jcount >= len(CTCF_uniq):
            break

