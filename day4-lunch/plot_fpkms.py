#!/usr/bin/env python3

"""
Usage: plot_fpkms.py <sample1/t_data.ctab> <sample2/t_data.ctab>

make a scatter plot of FPKM values of two samples
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df1 = pd.read_csv(sys.argv[1], sep="\t", index_col="t_name")
df2 = pd.read_csv(sys.argv[2], sep="\t", index_col="t_name")

fpkm1 = df1.loc[:,"FPKM"]
fpkm2 = df2.loc[:,"FPKM"]

#to find log values to make plot log scale
fpkm1_log = np.log(fpkm1 + 1)
fpkm2_log = np.log(fpkm2 + 1)

#to print scatter plot
fig, ax = plt.subplots()
ax.scatter(fpkm1_log, fpkm2_log, alpha=0.05)

#title and axis labels
ax.set_title("FPKM comparison")
ax.set_xlabel("log fpkm1")
ax.set_ylabel("log fpkm2")

#for fit line
p = np.polyfit(fpkm1_log,fpkm2_log,1) #returns coefficients for mx+b (since deg 1 was used)
poly = np.poly1d(p) #returns function that you can use
x = np.linspace(min(fpkm1_log), max(fpkm1_log)) #determine numbers by eye (look for what values xaxis covers)
plt.plot(x, poly(x), 'k') #'k' changes line to be bold. poly(x) gives y values



#these are the numbers for polyfit


fig.savefig("FinalPlot.png")
plt.close(fig)
