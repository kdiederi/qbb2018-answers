#!/usr/bin/env python3

"""
Usage: merge.py <threshold> <sample1/t_data.ctab> <sample2/t_data.ctab> ... <samplei/t_data.ctab> 

Create csv file with FPKMs from i samples. then sum across rows and check if sum is greater than input threshold.
print if sum fpkm is greater than threshold
need to specify in input what directory files are in
"""

import sys
import os
import pandas as pd
import numpy as np

list_of_files = sys.argv[2:len(sys.argv)]
fpkm_dict = {} 

#this iterates through the files and adds sample name and all fpkms to the dictionary
for file in list_of_files:
    name = file.split(os.sep)[-2]
    fpkm = pd.read_csv(file, sep="\t", index_col="t_name").loc[:,"FPKM"]
    fpkm_dict[name] = fpkm
#makes the dictionary into a dataframe (joins each individual sample column together with others)    

fpkms_all = pd.DataFrame(fpkm_dict)

#need to sum accross rows
fpkm_sum = fpkms_all.sum(axis=1)
#to add sums column to data frame
fpkm_sum_col = fpkms_all.assign(sums=fpkm_sum)

over_threshold = fpkm_sum_col.loc[:,"sums"] > float(sys.argv[1])
print(fpkm_sum_col.loc[over_threshold,:])
