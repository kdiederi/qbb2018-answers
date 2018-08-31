#!/usr/bin/env python3

"""
Usage: ./timecourse.py <gene_name> <samples.csv> <ctab_dir> 

create a timecourse of a given gene for females and males. average FPKM
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



def timecourse(gender):
    df = pd.read_csv(sys.argv[2])
    #soi=samples of interest
    soi = df.loc[:,"sex"] == gender
    df = df.loc[soi,:] 
    stages = []
    fpkms = []
    for index, sample, sex, stage in df.itertuples():
        filename = os.path.join(sys.argv[3],sample,"t_data.ctab",)
        ctab_df = pd.read_table(filename,index_col="t_name")
        roi = ctab_df.loc[:,"gene_name"] == sys.argv[1]
        
        fpkm = ctab_df.loc[roi,"FPKM"]
        fpkms.append(ctab_df.loc[roi,"FPKM"].mean())
        stages.append(stage)
    return fpkms,stages
        

fpkms_female,stages = timecourse("female")
fpkms_male,stages = timecourse("male")



fig, ax = plt.subplots()
ax.plot(fpkms_female, label="female", c="red")
ax.plot(fpkms_male, label="male", c="blue")

ax.set_xlabel('developmental stage')
ax.set_ylabel('Mean mRNA abundance (RPKM)')
ax.set_title(str(sys.argv[1]))
#sets axis labels to the stages, rotates label
ax.set_xticklabels(stages, rotation=90)
#legend. need to have labels specified in ax.plot
plt.legend(bbox_to_anchor=(1.05,0.55),loc=2, borderaxespad=0)
plt.tight_layout()
fig.savefig("geneSearch.png", bbox_inches='tight')
plt.close(fig)

