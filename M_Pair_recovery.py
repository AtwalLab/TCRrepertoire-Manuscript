"""
PAIRING OF ALPHA/BETA CHAINS.

This code inputs final dataframes acquired from filtering of raw data, and
matches alpha and beta chains via their droplets. Some droplets won't have
a matched pair, and will not be included in the final paired dataframe.
Make sure to be in correct location before running the code.
"""

from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

#Import single chain files.
TRA_prod=pd.read_csv('TRAs_Final_productive.txt', sep='\t')
TRB_prod=pd.read_csv('TRBs_Final_productive.txt', sep='\t')

#Find all the droplets that appear in both alpha and beta dataframes.
#It's done by concatenating alpha and beta dataframes and taking only the droplets
#that appear twice in the resulting concatenated dataframe.
TRAB_prod = pd.concat([TRA_prod,TRB_prod],axis=0)

Droplet_prod = pd.DataFrame(TRAB_prod.droplet.value_counts())

Droplet_list_prod=Droplet_prod[Droplet_prod['droplet']==2].index.tolist()


#Make empty dataframes to fill in with alpha/beta chains for pairing later.
TRA_forpairing_prod=pd.DataFrame(index=range(len(Droplet_list_prod)),columns=TRA_prod.columns.tolist())
TRB_forpairing_prod=pd.DataFrame(index=range(len(Droplet_list_prod)),columns=TRB_prod.columns.tolist())


#Fill in the dataframes.
for i in range(len(Droplet_list_prod)):
    TRA_forpairing_prod.iloc[i]=TRA_prod[TRA_prod['droplet']==Droplet_list_prod[i]].iloc[0]

for i in range(len(Droplet_list_prod)):
    TRB_forpairing_prod.iloc[i]=TRB_prod[TRB_prod['droplet']==Droplet_list_prod[i]].iloc[0]


#Make columns names for alpha and beta chains, for distinguishing in the paired dataframe.
Alpha_columns = ['dropletA','molecularA','numberA','sequenceA',
'BestVhitA','BestJhitA','NSeqCDR1A','LengthCDR1A','NSeqCDR2A','LengthCDR2A',
'NSeqCDR3A','LengthCDR3A','NSeqFR1A','LengthFR1A','NSeqFR2A','LengthFR2A','NSeqFR3A',
'LengthFR3A','NSeqFR4A','LengthFR4A','AASeqCDR1A','Length_CDR1A','AASeqCDR2A',
'Length_CDR2A','AASeqCDR3A','Length_CDR3A','AASeqFR1A','Length_FR1A','AASeqFR2A',
'Length_FR2A','AASeqFR3A','Length_FR3A','AASeqFR4A','Length_FR4A']

Beta_columns = ['dropletB','molecularB','numberB','sequenceB',
'BestVhitB','BestDhitB','BestJhitB','NSeqCDR1B','LengthCDR1B','NSeqCDR2B','LengthCDR2B',
'NSeqCDR3B','LengthCDR3B','NSeqFR1B','LengthFR1B','NSeqFR2B','LengthFR2B','NSeqFR3B',
'LengthFR3B','NSeqFR4B','LengthFR4B','AASeqCDR1B','Length_CDR1B','AASeqCDR2B',
'Length_CDR2B','AASeqCDR3B','Length_CDR3B','AASeqFR1B','Length_FR1B','AASeqFR2B',
'Length_FR2B','AASeqFR3B','Length_FR3B','AASeqFR4B','Length_FR4B']

#Rename the columns of resulting dataframes before concatenating.
TRA_forpairing_prod.columns=Alpha_columns
TRB_forpairing_prod.columns=Beta_columns


#Concatenate alpha and beta dataframes into one dataframe of paired information.
TR_pairs_prod=pd.concat([TRA_forpairing_prod,TRB_forpairing_prod], axis=1)

#Save files as .txt.
TR_pairs_prod.to_csv('TRpairs_productive.txt', sep='\t', index=False)
