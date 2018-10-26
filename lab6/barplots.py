#!/usr/bin/env python3

"""
create a two panel plot. one panel for number of sites gained and lost upon differentiation, the other for number of sites in each type of region
 
./barplots.py <peaksgained.bed> <peakslost.bed> <Mus_musculus_features.bed> <G1E.narrowpeaks> <ER4.narrowpeaks>
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

#count number of peaks gained and peaks lost upon cell differation
peaks_gained = 0
peaks_lost = 0
for line in open(sys.argv[1]):
    peaks_gained += 1
for line in open(sys.argv[2]):
    peaks_lost += 1
variables = ['Peaks Gained', 'Peaks Lost']
x_pos = np.arange(len(variables))
data = [peaks_gained, peaks_lost]
print("change in peaks calculated")

#count number of features overlapping with G1E and ER4
features = open(sys.argv[3])
G1E_cells = open(sys.argv[4])
ER4_cells = open(sys.argv[5])

G1E_exon = 0
G1E_intron = 0
G1E_promoter = 0
ER4_exon = 0
ER4_intron = 0
ER4_promoter = 0

class G1E:
    def __init__(self,start,end):
        self.start = start
        self.end = end
    def getstart(self):
        return self.start
    def getend(self):
        return self.end
        
    def __str__(self):
        #formatted this way to make the printing look good. and tab delimited
        return (self.start+ "\t" + self.end)
        
G1E_positions = []
    
for line in G1E_cells:
    rfields = line.strip().split()
    G1E_positions.append(G1E(rfields[1],rfields[2]))
    
G1E_cells.close()

for each in G1E_positions:
    features.seek(0,0)
    for line2 in features:
        fields2 = line2.strip().split()
        if each.getstart() >= fields2[1] and each.getend() <= fields2[2]:
            if fields2[3] == 'promoter':
                G1E_promoter += 1
            if fields2[3] == 'exon':
                G1E_exon += 1
            if fields2[3] == 'intron':
                G1E_intron += 1

print('G1E exons: ' + str(G1E_exon))
print('G1E introns: ' + str(G1E_intron))
print('G1E promoters: ' + str(G1E_promoter))

class ER4:
    def __init__(self,start,end):
        self.start = start
        self.end = end
    def getstart(self):
        return self.start
    def getend(self):
        return self.end
        
    def __str__(self):
        #formatted this way to make the printing look good. and tab delimited
        return (self.start+ "\t" + self.end)
        
ER4_positions = []
    
for line3 in ER4_cells:
    fields3 = line3.strip().split()
    ER4_positions.append(ER4(fields3[1],fields3[2]))
    
ER4_cells.close()

for each2 in ER4_positions:
    features.seek(0,0)
    for line4 in features:
        fields4 = line4.strip().split()
        if each2.getstart() >= fields4[1] and each2.getend() <= fields4[2]:
            if fields4[3] == 'promoter':
                ER4_promoter += 1
            if fields4[3] == 'exon':
                ER4_exon += 1
            if fields4[3] == 'intron':
                ER4_intron += 1

print('ER4 exons: ' + str(ER4_exon))
print('ER4 introns: ' + str(ER4_intron))
print('ER4 promoters: ' + str(ER4_promoter))



features.close()

#plot data
fig, (ax2, ax1) = plt.subplots(ncols=2, figsize=(20,10))
ax1.bar(x_pos, data, align='center', color='royalblue', capsize=15)
ax1.set_ylabel('Number of peaks')
ax1.set_xticks(x_pos)
ax1.set_xticklabels(variables)
ax1.set_title('Changes in CTCF binding \n upon cell differentiation')
ax1.yaxis.grid(True, alpha=0.5)
ax1.set_axisbelow(True)

site_types = ['Intron', 'Exon', 'Promoter']
G1E_data = [G1E_intron, G1E_exon, G1E_promoter]
ER4_data = [ER4_intron, ER4_exon, ER4_promoter]
index = np.arange(len(site_types))
bar_width = 0.35

ER4data = ax2.bar(index, ER4_data, bar_width, color='royalblue', capsize=15, label='ER4 cells')
G1Edata = ax2.bar(index + bar_width, G1E_data, bar_width, color='magenta', capsize=15, label='G1E cells')
# ax2.bar(x_pos, data, align='center', color='royalblue', capsize=15)
ax2.set_ylabel('Count')
ax2.set_xlabel('Site type')
ax2.set_xticks(index + bar_width/2)
ax2.set_xticklabels(site_types)
ax2.set_title('Number of exons, introns, and promoters \n before and after cell differentiation')
ax2.yaxis.grid(True, alpha=0.5)
ax2.legend()
ax2.set_axisbelow(True)

##  Save the figure
# plt.tight_layout()
plt.savefig('plot1.png')
plt.close()