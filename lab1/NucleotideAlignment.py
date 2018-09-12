#!/usr/bin/env python3

"""
align nucleotides guided by protein alignment. considers both AA alignment
and the original sequence. gaps notated by -. then each codon is compared to query
and dN,dS values are calculated for each codon. dN/dS and dN-dS also calculated for each position.
plot made for dN/dS vs codon position
pipe to output file of choice (needed for dNdSZtestPlot.py)
./NucleotideAlignment.py <blast file> <MAFFT file> > <outputfile>
for this lab command is: ./NucleotideAlignmnet.py alignmentseq aaalignment > nucalign_output
"""

import sys
import fasta
import itertools
import pandas as pd
import matplotlib.pyplot as plt

#used during testing to slow down printing
#import time


dna_reader = fasta.FASTAReader(open(sys.argv[1]))
#dna imput is blast output file (this exercise its called alignmentseq)
aa_reader = fasta.FASTAReader(open(sys.argv[2]))
#aa imput is MAFFT output file (this exercise its called aaalignment)


nuc_seqs = []
aa_lists = []

for (dna_id, dna), (aa_id, aa) in zip(dna_reader,aa_reader):
    nuc_seq = []
    aa_list = []
    j = 0
    for i in range(len(aa)):
        a = aa[i]
        nuc = dna[j:j+3]
        aa_list.append(a)
        if a == "-":
            nuc = "---"
            nuc_seq.append(nuc)
        else:
            j += 3
            nuc_seq.append(nuc)
    nuc_seqs.append(nuc_seq)
    aa_lists.append(aa_list)



#separate query sequences from alignment sequences so you don't compare query to query for the first iteration
query_nuc = nuc_seqs[0]
query_aa = aa_lists[0]
nuc_match = nuc_seqs[1:]
aa_match = aa_lists[1:]

#keep track of dN/dS/noChange counts per codon
class CodonComparison:
    #define the subgroups in Mutation class (variables i am keeping count for, for each codon). keep three counts for each codon separate from other codons (why not to use count method)
    def __init__(self,dN,dS,no_change,codonpos,dNdSratio,dNdSdiff):
        self.dN = dN
        self.dS = dS
        self.no_change = no_change
        self.codonpos = codonpos
        self.dNdSratio = float(dNdSratio)
        self.dNdSdiff = float(dNdSdiff)
        
        #make function to increase count by one when called        
    def increasedN(self):
        self.dN = self.dN + 1
    def increasedS(self):
        self.dS = self.dS + 1
    def noChange(self):
        self.no_change = self.no_change + 1
    #to calculate dN/dS for codons with dS > 0
    def ratiocalc(self):
        if self.dS > 0:
            self.dNdSratio = float(self.dN) / float(self.dS)
    # #to include codons where dS= 0 ended up not ploting this
  #   def ratiocalc(self):
  #       self.dNdSratio = (float(self.dN) + 0.00001) / (float(self.dS) + 0.00001)
    def __str__(self):
        #formatted this way to make the printing look good. and tab delimited
        return (str(self.codonpos)+ "\t" + str(self.dN) + "\t" + str(self.dS) + "\t"
        + str(self.no_change) + "\t" + str(self.dNdSratio) + "\t" + str(self.dNdSdiff))
    #to make dictonary for dN/dSplot
    def to_dict(self):
        return {
            'x':self.codonpos,
            'y':self.dNdSratio
                }
    def diff(self):
        self.dNdSdiff = float(self.dN) - float(self.dS)

#iterate through matched sequences and compare to query sequence.
#codoncomparison will be aded to results (results= list containing codoncomparison, which has the dN/dS/no change counts for that specific codon)
results = []
#start with first position having zeros
results.append(CodonComparison(0,0,0,1,0,0))

for codon,aa in zip(nuc_match,aa_match):
    for r in range(len(aa)):
        #if this is the first time a codon is being looked at, we need to set each value to zero (initialize it).
        if len(results) < r + 1 :
             results.append(CodonComparison(0,0,0,r+1,0,0))
        #double check that the list is getting initialized properly
        # print(len(results))
        # print(r)
        # time.sleep(10)
        
        #set variables to compare    
        a = aa[r]
        c = codon[r]
        #to see if the variables are what I think they are
        # print(str(r))
        # print(a)
        # print(c)
        # print(query_aa[r])
        # time.sleep(10)
        
#if the aa in test is different from query it is for sure a non-synonomous mutation
        if a != query_aa[r]:
            results[r].increasedN()
#if the aa in test is the same from quary it might be a nt mutation (which is a synonomus mut) or there might not be a mutation (no change)
        else:
            if c == query_nuc[r]:
                results[r].noChange()
            else:
                results[r].increasedS()
   

print("codon_pos" + "\t" + "dN" + "\t" + "dS" + "\t" + "no_change" + "\t" + "dN/dS" +  "\t" + "dNdS_diff")
for result in results:
    result.ratiocalc()
    result.diff()
    print(str(result))


#make dataframe from dictionary of codon pos and dN/dS https://stackoverflow.com/questions/34997174/how-to-convert-list-of-model-objects-to-pandas-dataframe/41762270#41762270

df = pd.DataFrame.from_records([res.to_dict() for res in results])

#plot dN/dS versus codon location

# Create a blank canvas
fig1, ax1 = plt.subplots()         
# Plot x vs y as points, from dataframe https://chrisalbon.com/python/data_visualization/matplotlib_scatterplot_from_pandas/
ax1.scatter(df.x,df.y)             
plt.ylabel("dN/dS")
plt.xlabel("Codon Number")
plt.title("Ratio of dN/dS for each codon")
# Save the figure
fig1.savefig("dNdSscatter.png") 
# Close the canvas
plt.close()      
