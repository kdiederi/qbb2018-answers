#!/usr/bin/env python3

"""
Usage: ./timecourse.py <samples.csv> <ctab_dir> <gene_name1> <gene_name2> ...<gene_namei> 

create a plot of FPKMs for many different genes for females and males
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


for search in sys.argv[3:]:
    
    def timecourse(gender,gene):
        df = pd.read_csv(sys.argv[1])
        #soi=samples of interest
        soi = df.loc[:,"sex"] == gender
        df = df.loc[soi,:]
        stages = []
        fpkms = []
        for index, sample, sex, stage in df.itertuples():
            filename = os.path.join(sys.argv[2],sample,"t_data.ctab",)
            ctab_df = pd.read_table(filename,index_col="t_name")
            
            
            roi = ctab_df.loc[:,"gene_name"] == gene

            fpkms.append(ctab_df.loc[roi,"FPKM"].mean())
            stages.append(stage)
        return fpkms,stages


    fpkms_female,stages = timecourse("female",search)
    fpkms_male,stages = timecourse("male",search)

    fig, ax = plt.subplots()
    ax.plot(fpkms_female, label="female", c="red")
    ax.plot(fpkms_male, label="male", c="blue")

    ax.set_xlabel('developmental stage')
    ax.set_ylabel('Mean mRNA abundance (RPKM)')
    ax.set_title(search)
    #sets axis labels to the stages, rotates label
    ax.set_xticklabels(stages, rotation=90)
    #legend. need to have labels specified in ax.plot
    plt.legend(bbox_to_anchor=(1.05,0.55),loc=2, borderaxespad=0)
    plt.tight_layout()
    fig.savefig(search, bbox_inches='tight')
    plt.close(fig)

