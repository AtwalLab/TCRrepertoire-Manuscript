"""
FILTERING. POST CLONALITY.
This code takes files after clonality analysis with mixcr, finds all identical
clones, collapses the sequences and adds to the original unique droplets
dataframes.

"""

from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd


#Import all necessary dataframes.
TRAm_prod=pd.read_csv('TRAm_productive.txt', sep='\t')
TRAsingle_prod=pd.read_csv('TRAsingle_productive.txt',sep='\t')
TRAm_prod_clones=pd.read_csv('TRAm_prod_clones.txt', sep='\t')

TRBm_prod=pd.read_csv('TRBm_productive.txt', sep='\t')
TRBsingle_prod=pd.read_csv('TRBsingle_productive.txt',sep='\t')
TRBm_prod_clones=pd.read_csv('TRBm_prod_clones.txt', sep='\t')

#Find dropped sequences if there are any and remove them from original dataframe.
while True:
    restart = False
    for i in range(len(TRAm_prod_clones)):
        if TRAm_prod_clones.droplet[i]!=TRAm_prod.droplet[i]:
            TRAm_prod=TRAm_prod.drop(TRAm_prod.index[i])
            TRAm_prod=TRAm_prod.reset_index(drop=True)
            restart = True
            break
    if not restart:
        break

while True:
    restart = False
    for i in range(len(TRBm_prod_clones)):
        if TRBm_prod_clones.droplet[i]!=TRBm_prod.droplet[i]:
            TRBm_prod=TRBm_prod.drop(TRBm_prod.index[i])
            TRBm_prod=TRBm_prod.reset_index(drop=True)
            restart = True
            break
    if not restart:
        break

#Insert clone ID into the main dataframes.
TRAm_prod.insert(0,'clone',TRAm_prod_clones.cloneId)
TRBm_prod.insert(0,'clone',TRBm_prod_clones.cloneId)


#Identify identical clones with the same droplet barcode.
#Get a list of unique droplet barcodes for each dataframe.
dropletA_prod=TRAm_prod.droplet.unique().tolist()
dropletB_prod=TRBm_prod.droplet.unique().tolist()

#Identify droplets that are associated with more than one clone ID.
#Make a new dataframe of those droplets by dropping droplets that are
#associated with unique clone ID only.
print 'collapsing identical clone sequences'
TRAm_new_prod = TRAm_prod
for i in dropletA_prod:
    a=TRAm_prod.loc[TRAm_prod['droplet'] == i]
    if len(a.clone.value_counts())==1:
        TRAm_new_prod = TRAm_new_prod[TRAm_new_prod.droplet != i]

TRBm_new_prod = TRBm_prod
for i in dropletB_prod:
    a=TRBm_prod.loc[TRBm_prod['droplet'] == i]
    if len(a.clone.value_counts())==1:
        TRBm_new_prod = TRBm_new_prod[TRBm_new_prod.droplet != i]


#Drop the extracted indeces from the original dataframe to acquire a dataframe
#of droplets associated with a single clone.
TRAs_prod=TRAm_prod.drop(TRAm_prod.index[TRAm_new_prod.index.tolist()])
TRBs_prod=TRBm_prod.drop(TRBm_prod.index[TRBm_new_prod.index.tolist()])

#Now take those single droplet dataframes and collapse the sequences.
#The new count is the resulting sum of those sequences.
#The sequence that goes into the new dataframe is the one that is most common.

#Make a list of unique droplets.
TRAs_droplet_prod=TRAs_prod.droplet.unique().tolist()
TRBs_droplet_prod=TRBs_prod.droplet.unique().tolist()

#Create empty dataframes to fill in with collapsed sequences.
TRAs_new_prod=pd.DataFrame(index=range(0,len(TRAs_droplet_prod)),columns=TRAm_prod.columns.tolist())
TRBs_new_prod=pd.DataFrame(index=range(0,len(TRBs_droplet_prod)),columns=TRBm_prod.columns.tolist())

print 'collapsing identical clone sequences - filling in dataframes'
#Fill in the dataframe with collapsed sequences. Sum the counts of all sequences.
for i in range(0,len(TRAs_droplet_prod)):
    a=TRAs_prod.loc[TRAs_prod['droplet'] == TRAs_droplet_prod[i]]
    TRAs_new_prod.iloc[i]=a[a.number == max(a.number)].iloc[0]
    TRAs_new_prod.number.iloc[i]=sum(a.number)

for i in range(0,len(TRBs_droplet_prod)):
    a=TRBs_prod.loc[TRBs_prod['droplet'] == TRBs_droplet_prod[i]]
    TRBs_new_prod.iloc[i]=a[a.number == max(a.number)].iloc[0]
    TRBs_new_prod.number.iloc[i]=sum(a.number)


#Drop the clone column from unique droplet dataframes.
TRAs_new_prod=TRAs_new_prod.drop('clone', 1)
TRBs_new_prod=TRBs_new_prod.drop('clone', 1)

#Add the resulting unique droplet dataframe to the single droplet dataframe from FILTER3.
TRA_single_prod=pd.concat([TRAsingle_prod,TRAs_new_prod],axis=0)
TRB_single_prod=pd.concat([TRBsingle_prod,TRBs_new_prod],axis=0)

#Save the dataframes as .txt file.
TRA_single_prod.to_csv('TRAs_final_productive.txt', sep='\t', index=False)
TRB_single_prod.to_csv('TRBs_final_productive.txt', sep='\t', index=False)
TRA_single_unprod.to_csv('TRAs_final_unproductive.txt', sep='\t', index=False)
TRB_single_unprod.to_csv('TRBs_final_unproductive.txt', sep='\t', index=False)
