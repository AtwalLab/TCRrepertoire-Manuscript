"""
ADD CLONE ID TO EACH SEQUENCE.
This code takes the clone id assigned to each sequence from mixcr clonality
analyses and adds it to the final dataframe. It is a way of identifying sequences
that are different in nucleotides, but considered the same clone.

If sequences get dropped during clonality analysis, the code finds sequences with
no assigned ID and drops them from the dataframe.

"""

from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

#import files
TRAs_prod=pd.read_csv('TRAsFinal_productive.txt', sep='\t')
TRAs_prod_clones=pd.read_csv('TRAs_final_prod_clones.txt', sep='\t')

TRBs_prod=pd.read_csv('TRBsFinal_productive.txt', sep='\t')
TRBs_prod_clones=pd.read_csv('TRBs_final_prod_clones.txt', sep='\t')


#Split description column in files from mixcr
DescAs_prod=TRAs_prod_clones['descrR1'].str.split('_', 3, expand=True)
DescBs_prod=TRBs_prod_clones['descrR1'].str.split('_', 3, expand=True)

#Extract only wanted columns
DescAs_prod = DescAs_prod[DescAs_prod.columns.tolist()[1:4]]
DescBs_prod = DescBs_prod[DescBs_prod.columns.tolist()[1:4]]


#Rename columns
DescAs_prod.columns=['droplet', 'molecular','number']
DescBs_prod.columns=['droplet', 'molecular','number']


#Add droplet column to clone ID.
TRAs_prod_clones=pd.concat([DescAs_prod[['droplet']],TRAs_prod_clones.drop('descrR1', 1)],axis=1)
TRBs_prod_clones=pd.concat([DescBs_prod[['droplet']],TRBs_prod_clones.drop('descrR1', 1)],axis=1)


#Find dropped sequences if there are any and remove them from original dataframe.
while True:
    restart = False
    for i in range(len(TRAs_prod_clones)):
        if TRAs_prod_clones.droplet[i]!=TRAs_prod.droplet[i]:
            TRAs_prod=TRAs_prod.drop(TRAs_prod.index[i])
            TRAs_prod=TRAs_prod.reset_index(drop=True)
            restart = True
            break
    if not restart:
        break

while True:
    restart = False
    for i in range(len(TRBs_prod_clones)):
        if TRBs_prod_clones.droplet[i]!=TRBs_prod.droplet[i]:
            TRBs_prod=TRBs_prod.drop(TRBs_prod.index[i])
            TRBs_prod=TRBs_prod.reset_index(drop=True)
            restart = True
            break
    if not restart:
        break


#Add clone ID to full dataframe as a first column.
TRAsFinal_prod=pd.concat([TRAs_prod_clones[['cloneId']],TRAs_prod], axis=1)
TRBsFinal_prod=pd.concat([TRBs_prod_clones[['cloneId']],TRBs_prod], axis=1)

#Save new dataframes.
TRAsFinal_prod.to_csv('TRAs_Final_productive.txt', sep='\t', index=False)
TRBsFinal_prod.to_csv('TRBs_Final_productive.txt', sep='\t', index=False)
