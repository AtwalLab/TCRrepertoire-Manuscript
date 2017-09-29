"""
ASSIGN CD4 OR CD8 TYPE.
This code takes CD4/CD8 tables, extracts all the droplets that have an assigned
type and makes new separate CD4/CD8 dataframes.

Make sure to check what CD type table look like in pandas dataframe before running
the code. May vary subject to subject.
"""

from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
import glob


#Make a list of CD type files.
file_list=glob.glob('*.tab')

#Initialize CD type dataframe by importing the first file.
CD=pd.read_table(file_list[0])
#Loop over the remaining files to add to the initialized dataframe.
for i in range(1,len(file_list)):
    a=pd.read_table(file_list[i])
    CD=pd.concat([CD,a], axis=0)

#Import all files for stratification.
TR_prod=pd.read_csv('TRpairs_productive.txt', sep='\t')

TRA_prod=pd.read_csv('TRAs_Final_productive.txt', sep='\t')
TRB_prod=pd.read_csv('TRBs_Final_productive.txt', sep='\t')

#Rename the columns to more convenient names.
CD=CD[['DB','CD4','CD8']]

#Select droplets that correspond to CD4 or CD8.
CD4=((CD[CD.CD4>0])[(CD[CD.CD4>0]).CD8==0]).reset_index(drop=True)
CD8=((CD[CD.CD8>0])[(CD[CD.CD8>0]).CD4==0]).reset_index(drop=True)

#Get a unique list of droplets assigned to CD type.
CD4=CD4.DB.unique().tolist()
CD8=CD8.DB.unique().tolist()

#Create empty dataframes cd4/cd8
TRcd4_prod=pd.DataFrame(index=range(0,len(CD4)),columns=TR_prod.columns.tolist())
TRcd8_prod=pd.DataFrame(index=range(0,len(CD8)),columns=TR_prod.columns.tolist())

TRAcd4_prod=pd.DataFrame(index=range(0,len(CD4)),columns=TRA_prod.columns.tolist())
TRAcd8_prod=pd.DataFrame(index=range(0,len(CD8)),columns=TRA_prod.columns.tolist())

TRBcd4_prod=pd.DataFrame(index=range(0,len(CD4)),columns=TRB_prod.columns.tolist())
TRBcd8_prod=pd.DataFrame(index=range(0,len(CD8)),columns=TRB_prod.columns.tolist())

#Fill in the dataframes
print 'Paired dataframe'
for i in range(0,len(CD4)):
    a=TR_prod[TR_prod.dropletA==CD4[i]]
    if len(a)==1:
        TRcd4_prod.iloc[i]=a.iloc[0]
for i in range(0,len(CD8)):
    a=TR_prod[TR_prod.dropletA==CD8[i]]
    if len(a)==1:
        TRcd8_prod.iloc[i]=a.iloc[0]


print 'Alpha dataframe'
for i in range(0,len(CD4)):
    a=TRA_prod[TRA_prod.droplet==CD4[i]]
    if len(a)==1:
        TRAcd4_prod.iloc[i]=a.iloc[0]
for i in range(0,len(CD8)):
    a=TRA_prod[TRA_prod.droplet==CD8[i]]
    if len(a)==1:
        TRAcd8_prod.iloc[i]=a.iloc[0]


print 'Beta dataframe'
for i in range(0,len(CD4)):
    a=TRB_prod[TRB_prod.droplet==CD4[i]]
    if len(a)==1:
        TRBcd4_prod.iloc[i]=a.iloc[0]
for i in range(0,len(CD8)):
    a=TRB_prod[TRB_prod.droplet==CD8[i]]
    if len(a)==1:
        TRBcd8_prod.iloc[i]=a.iloc[0]


#Drop empty rows.
TRcd4_prod = TRcd4_prod.dropna(subset=['dropletA'])
TRcd8_prod = TRcd8_prod.dropna(subset=['dropletA'])

TRAcd4_prod = TRAcd4_prod.dropna(subset=['droplet'])
TRAcd8_prod = TRAcd8_prod.dropna(subset=['droplet'])

TRBcd4_prod = TRBcd4_prod.dropna(subset=['droplet'])
TRBcd8_prod = TRBcd8_prod.dropna(subset=['droplet'])

TRAcd4_prod = TRAcd4_prod.dropna(subset=['droplet'])
TRAcd8_prod = TRAcd8_prod.dropna(subset=['droplet'])


#Save the dataframes as txt files.
TRcd4_prod.to_csv('TRcd4_productive.txt', sep='\t', index=False)
TRcd8_prod.to_csv('TRcd8_productive.txt', sep='\t', index=False)

TRAcd4_prod.to_csv('TRAcd4_productive.txt', sep='\t', index=False)
TRAcd8_prod.to_csv('TRAcd8_productive.txt', sep='\t', index=False)

TRBcd4_prod.to_csv('TRBcd4_productive.txt', sep='\t', index=False)
TRBcd8_prod.to_csv('TRBcd8_productive.txt', sep='\t', index=False)
