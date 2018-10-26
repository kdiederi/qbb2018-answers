#!/usr/bin/env python3

"""
density plot of motif starts
 
./motifdensityplot.py <motifpositons.txt>
"""

import sys
import matplotlib.pyplot as plt
import pandas as pd


dictionary = {}

for line in open(sys.argv[1]):
    fields = line.strip().split()
    position = int(fields[2])
    if position not in dictionary:
        dictionary[position] = 1
    else:
        dictionary[position] += 1

fig,ax = plt.subplots()
for key in dictionary:
    plt.scatter([key], dictionary[key])

ax.set_title("Location of motif matches in input sequences")
ax.set_xlabel("Position in sequence")
ax.set_ylabel("Counts")
fig.savefig("Plot.png")
plt.close()

