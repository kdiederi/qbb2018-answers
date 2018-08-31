#!/usr/bin/env python3

"""
Usage: ./LinearRegressionbigWigAvg.py <bigWig.tab1><bigWig.tab2><bigWig.tab3>
<bigWig.tab4><bigWig.tab5> <ctab_file>

create a linear regression for all of the five big wig files generated previously
reports R2, p-value, coefficients
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pandas as pd

#list of files and what columns/values we want from them
mean1 = pd.read_table(sys.argv[1], index_col=0).iloc[:,5-1]
mean2 = pd.read_table(sys.argv[2], index_col=0).iloc[:,5-1]
mean3 = pd.read_table(sys.argv[3], index_col=0).iloc[:,5-1]
mean4 = pd.read_table(sys.argv[4], index_col=0).iloc[:,5-1]
mean5 = pd.read_table(sys.argv[5], index_col=0).iloc[:,5-1]
ctab = pd.read_table(sys.argv[6], index_col='t_name').loc[:,"FPKM"]

#make the keys the title of the file name (keys that relate to columns of interest)
name1 = sys.argv[1].split(os.sep)[-1]
name2 = sys.argv[2].split(os.sep)[-1]
name3 = sys.argv[3].split(os.sep)[-1]
name4 = sys.argv[4].split(os.sep)[-1]
name5 = sys.argv[5].split(os.sep)[-1]
ctab_name = sys.argv[6].split(os.sep)[-1]

d = {name1:mean1,
     name2:mean2,
     name3:mean3,
     name4:mean4,
     name5:mean5,
     ctab_name:ctab}

df = pd.DataFrame(d)
#look through dataframe and find rows that have means of not a number (promoter start or end were <0 so ignored when making mean files. ctab still has those FPKMS and so pandas but NA into the dataframe). for roi look at one of the mean columns and get rid of those that dont have values
roi = df.loc[:,name1].notna()

small_df = df.loc[roi,:]

y = small_df.loc[:,ctab_name]
x = small_df.loc[:,[name1,name2,name3,name4,name5]]

small_df.columns = small_df.columns.str.replace(".","_")
# to check column names to use in formula command
# print(list(small_df.columns.values))
#formula= response variable ~ imput variable 1 + imput variable 2#
mod = smf.ols(formula="t_data_ctab ~ H3k27ac_tab + H3k27me3_tab + H3k4me1_tab + H3k4me3_tab + H3k9ac_tab", data=small_df)
res = mod.fit()
print(res.summary())