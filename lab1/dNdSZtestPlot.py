#!/usr/bin/env python3

"""
calculate standard error of dN-dS across all codon positions
compare dN and dS at each codon with Z test. Plot Z score versus codon position 

./dNdSZtestPlot.py <NucleotideAlightment output file> 
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#set up dataframe from nucleotide alignment file (which is tab delimited)
df = pd.read_table(sys.argv[1], sep='\t')

#to calculate Z test for each codon
    # z = (Dc - Dn) / SE where Dc is the difference between dN and dS at that codon,
  #   Dn is the null hypothesis, and SE is the standard error of the differences between dN and dS across all the codon positions.
dNdSdiff = df.iloc[:,5]
dNdSdiff_stddev = dNdSdiff.std()
# print("dNdSdiff_stddev : " + str(dNdSdiff_stddev))
Dc = dNdSdiff
Dn = 0
#n is sample size, for this exercise it was 1247
n = 1247
stderror = dNdSdiff_stddev / np.sqrt(n)
# print("stderror: " + str(stderror))
ztest = (Dc - Dn) / stderror
# print(ztest)

newdf = pd.concat([df,ztest.rename('ztest')], axis=1)

# print(newdf)

#plot ztest versus codon position
plt.scatter(newdf.iloc[:,0],newdf.iloc[:,6])
plt.xlabel("Codon Number")
plt.ylabel("Z score")
plt.title("Z-test of dN-dS at every codon position")
plt.savefig("ZtestPlot.png")
plt.close()