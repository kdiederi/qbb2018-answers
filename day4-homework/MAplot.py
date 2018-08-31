#!/usr/bin/env python3

"""
Usage: plot_fpkms.py <sample1/t_data.ctab> <sample2/t_data.ctab>

make a MA plot of FPKM values of two samples
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df1 = pd.read_csv(sys.argv[1], sep="\t", index_col="t_name")
df2 = pd.read_csv(sys.argv[2], sep="\t", index_col="t_name")

fpkm1 = df1.loc[:,"FPKM"]
fpkm2 = df2.loc[:,"FPKM"]


m_val = np.log2(fpkm1 + 0.1) - np.log2(fpkm2 + 0.1)
a_val = 0.5*(np.log2(fpkm1 + 0.1)+np.log2(fpkm2 + 0.1))




#to print scatter plot
fig, ax = plt.subplots()
ax.scatter(a_val, m_val, alpha=0.05)

#title and axis labels
ax.set_title("FPKM MA plot")
ax.set_xlabel("Average")
ax.set_ylabel("Ratio")
fig.savefig("MAPlot1.png")
plt.close(fig)

