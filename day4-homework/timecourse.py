#!/usr/bin/env python3

"""
Usage: ./timecourse.py <t_name> <samples.csv> <ctab_dir> <replicates.csv>

create a timecourse of a given transcript (FBtr0331261 for females and males)
"""

import sys
import os
import pandas as pd
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
    #subset out rows and columns we care about
        fpkms.append(ctab_df.loc[sys.argv[1],"FPKM"])
        stages.append(stage)
    return fpkms,stages
def timecourse_replicate(gender):
    df = pd.read_csv(sys.argv[4])
    #soi=samples of interest
    soi = df.loc[:,"sex"] == gender
    df = df.loc[soi,:]
    
    stages2 = []
    fpkms = []
    for index, sample, sex, stage in df.itertuples():
        filename = os.path.join(sys.argv[3],sample,"t_data.ctab",)
        ctab_df = pd.read_table(filename,index_col="t_name")
    #subset out rows and columns we care about
        fpkms.append(ctab_df.loc[sys.argv[1],"FPKM"])
        stages2.append(stage)
    return fpkms,stages2
    
fpkms_female,stages = timecourse("female")
fpkms_male,stages = timecourse("male")

fpkms_female_rep,stages2 = timecourse_replicate("female")
fpkms_male_rep,stages2 = timecourse_replicate("male")



fig, ax = plt.subplots()
ax.plot(fpkms_female, label="female", c="red")
ax.plot(fpkms_male, label="male", c="blue")
ax.plot([4,5,6,7],fpkms_female_rep, label="female replicate", c="green")
ax.plot([4,5,6,7],fpkms_male_rep, label="male replicate", c="orange")
ax.set_xlabel('developmental stage')
ax.set_ylabel('mRNA abundance (RPKM)')
ax.set_title('Sxl')
#sets axis labels to the stages, rotates label
ax.set_xticklabels(stages, rotation=90)
#legend. need to have labels specified in ax.plot
plt.legend(bbox_to_anchor=(1.05,0.55),loc=2, borderaxespad=0)
plt.tight_layout()
fig.savefig("timecourse1.png", bbox_inches='tight')
plt.close(fig)

