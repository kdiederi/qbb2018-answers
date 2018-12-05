#!/usr/bin/env python3

"""
sort genes, create dendogram and heatmap of gene expression
 
./clustering.py <hema_data.txt>
"""

import sys
from matplotlib import pyplot as plt
plt.style.use('ggplot')
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from sklearn.cluster import KMeans
from scipy import stats
import numpy as np
import pandas as pd
import csv


file = open(sys.argv[1])

header = file.readline()
celltype = header.split()
celltype_label = celltype[1:]


genes = []
data = []
for line in file:
    fields = line.split()
    genes.append(fields[0])
    data_row = []
    for fpkm in fields[1:]:
        data_row.append(float(fpkm))
    data.append(data_row)
    
X = np.array(data)
genes = np.array(genes)

Z = linkage(X, method='average', metric='euclidean', optimal_ordering=True) 
gene_index = leaves_list(Z)

Z_T = linkage(X.T, method='average', metric='euclidean', optimal_ordering=True) 
cell_index = leaves_list(Z_T)

fig = plt.figure(figsize=(25, 10))
dn = dendrogram(Z)
ax = plt.gca()
ax.set_xticklabels(genes[gene_index])
ax.set_xlabel("gene")
fig.savefig("dendrogram.png")


labels = np.array(celltype_label)
species = genes

sorted_X = X[gene_index]
sorted_X = sorted_X[:, cell_index]

fig2, ax2 = plt.subplots()          # Open a blank canvas, 8 inches x 6 inches
ax2.set_title("Heatmap of gene expression") # Add a title to the top
im = ax2.pcolor(                              # Treat the values like pixel intensities in a picture
	sorted_X,                                      # ... Using X as the values
	cmap="RdBu",                             # ... Use the Red-white-blue colormap to assign colors to your pixel values
)

ax2.grid(False)                      # Turn of the grid lines (a feature added automatically by ggplot)
ax2.set_xticks(                      # Edit the xticks being shown
	np.arange(0.5, X.shape[1]+0.5), # ... use the values centered on each column of pixels
	)
ax2.set_xticklabels(                 # Label the ticks
	labels[cell_index],                         # ... at position which correspond to the indices of our labels
	rotation=50,                    # ... and rotate the labels 50 degrees counter-clockwise
	)
ax2.set_yticks([])                   # Edit the ticks on the y-axis to show....NOTHING

cbar = fig.colorbar(im, ax=ax2)      # Add a bar to the right side of the plot which shows the scale correlating the colors to the pixel values

fig2.subplots_adjust( # Adjust the spacing of the subplots, to help make everything fit
    left = 0.05,     # ... the left edge of the left-most plot will be this percent of the way across the width of the plot
    bottom = 0.15,   # ... the bottom edge of the bottom-most plot will be this percent of the way up the canvas
    right = 1.0,     # ... the right edge of the right-most plot will be this percent of the way across the width
    top = 0.95,      # ... the top edge of the top-most plot will be this percent of the way from the bottom
)

fig2.savefig("heatmap.png") # Save the image
plt.close(fig2) # Close the canvas

#looked at heat map and saw 6 clusters. used this value in KMeans calculations

#kmeans for CFU and poly

cluster_kmeans = KMeans(n_clusters=6).fit(sorted_X)
index = cluster_kmeans.predict(sorted_X)
fig3, ax3 = plt.subplots()          # Open a blank canvas, 8 inches x 6 inches
ax3.scatter(sorted_X[:,0],sorted_X[:,5], c = index)
ax3.set_xlabel("CFU expression")
ax3.set_ylabel("poly expression");
ax3.set_title("K-means clustering") # Add a title to the top

fig3.savefig("kmeans2.png") # Save the image
plt.close(fig3) # Close the canvas


#t-test to identify genes differentially expressed between the two earliest stages in differentiation as compared to the two latest stages. 
#early stage cells= poly & unk
#late stage cells= mys &CFU

df = pd.read_csv(sys.argv[1], sep="\t", index_col="gene", header=0).dropna()
#average gene expression between cell types in each stage
df['early_avg'] = df[["poly","unk"]].mean(axis=1)
df['late_avg'] = df[["mys","CFU"]].mean(axis=1)

early = df[["poly","unk"]]
late = df[["mys","CFU"]]
stats, pval = stats.ttest_ind(early,late, axis=1)
df['pval'] = pval
#to check how many NaN pvalues (answer=2) 
# print(df.to_string())

diff_genes = []
panthergene_expression = 0
panthergene_name = ''


for index, row in df.iterrows():
    if (pd.isnull(row["pval"])):
        continue
    if row['pval'] <= 0.05 and row['early_avg'] < row['late_avg']:
        gene = row.name
        diff_genes.append(gene)
        if row['early_avg'] != 0:
            increase_expression = float(row['late_avg'])/float(row['early_avg'])
            if increase_expression > panthergene_expression:
                panthergene_expression = increase_expression
                panthergene_name = row.name
    else:
        continue

np.savetxt("diff_genes.csv", diff_genes, delimiter=",", fmt='%s', header=header)
print("The gene that is the most upregulated in the late differentiation stage is: " + panthergene_name) 
print("Expression level for this gene in late differentiation stage= " + str(panthergene_expression))