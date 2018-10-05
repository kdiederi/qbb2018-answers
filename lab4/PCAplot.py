#!/usr/bin/env python3

"""
plot PCA data from plink
 
./PCAplot.py <plink_output.eigenvec>
"""

import sys
import matplotlib.pyplot as plt

x = []
y = []
plink_file = open(sys.argv[1])
for line in plink_file:
    fields = line.strip("/r/n").split(" ")
    x.append(float(fields[2]))
    y.append(float(fields[3]))
    
    
fig1, ax1 = plt.subplots()
ax1.scatter(x,y)
ax1.set_xlabel("PCA1")
ax1.set_ylabel("PCA2")
ax1.set_title("PCA of S.cerevsiae strains")
fig1.savefig("PCAplot.png")
plt.close()