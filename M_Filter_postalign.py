"""

FILTERING. Post alignment.

PART A.
This code takes the MIXCR aligned sequence .txt files and creates productive.txt for each chain.
Productive = has all V(D)J hits, has all CDRs and FRs, has a stop codon only in FR4.
Make sure to be in the correct file location.

PART B.
This code creates two dataframes:
1. All droplets associated with only one sequence.
2. All droplets associated with more than one sequence.
This code is in preparation for clonality analysis to find sequences that
are identical under a given sequencing error rate.
"""

from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd

#Read aligned sequence files (made with Filter1.)
TRA = pd.read_csv('TRAaligned.txt', sep='\t')
TRB = pd.read_csv('TRBaligned.txt', sep='\t')

#Handle alpha chain.
TRA = TRA.drop(['bestDHit'], 1) #alpha doesn't have D gene.
#Removes rows that have any missing data (NaN)
TRA = TRA.dropna(subset=TRA.columns.tolist()).reset_index(drop=True)

#Handle beta chain.
TRBsubset=TRB.columns.tolist()
TRBsubset.remove('bestDHit') #beta chain may have D region mutated. Keep all.
#Removes rows with missing data (NaN), except for D gene.
TRB = TRB.dropna(subset=TRBsubset).reset_index(drop=True)

# Select sequences that don't have stop codons in forbidden regions.
# Also make dataframes for each forbidden region that have stop codons.
# alpha chain
a_CDR1=TRA[TRA.aaSeqCDR1.str.find('_')!=-1] #CDR1 for unproductive
TRAp=TRA[TRA.aaSeqCDR1.str.find('_')==-1]
a_CDR2=TRAp[TRAp.aaSeqCDR2.str.find('_')!=-1] #CDR2 for unproductive
TRAp=TRAp[TRAp.aaSeqCDR2.str.find('_')==-1]
a_CDR3=TRAp[TRAp.aaSeqCDR3.str.find('_')!=-1] #CDR3 for unproductive
TRAp=TRAp[TRAp.aaSeqCDR3.str.find('_')==-1]
a_FR1=TRAp[TRAp.aaSeqFR1.str.find('_')!=-1] #FR1 for unproductive
TRAp=TRAp[TRAp.aaSeqFR1.str.find('_')==-1]
a_FR2=TRAp[TRAp.aaSeqFR2.str.find('_')!=-1] #FR2 for unproductive
TRAp=TRAp[TRAp.aaSeqFR2.str.find('_')==-1]
a_FR3=TRAp[TRAp.aaSeqFR3.str.find('_')!=-1] #FR3 for unproductive
TRA_prod=TRAp[TRAp.aaSeqFR3.str.find('_')==-1] # Final productive alpha

#beta chain
b_CDR1=TRB[TRB.aaSeqCDR1.str.find('_')!=-1] #CDR1 for unproductive
TRBp=TRB[TRB.aaSeqCDR1.str.find('_')==-1]
b_CDR2=TRBp[TRBp.aaSeqCDR2.str.find('_')!=-1] #CDR2 for unproductive
TRBp=TRBp[TRBp.aaSeqCDR2.str.find('_')==-1]
b_CDR3=TRBp[TRBp.aaSeqCDR3.str.find('_')!=-1] #CDR3 for unproductive
TRBp=TRBp[TRBp.aaSeqCDR3.str.find('_')==-1]
b_FR1=TRBp[TRBp.aaSeqFR1.str.find('_')!=-1] #FR1 for unproductive
TRBp=TRBp[TRBp.aaSeqFR1.str.find('_')==-1]
b_FR2=TRBp[TRBp.aaSeqFR2.str.find('_')!=-1] #FR2 for unproductive
TRBp=TRBp[TRBp.aaSeqFR2.str.find('_')==-1]
b_FR3=TRBp[TRBp.aaSeqFR3.str.find('_')!=-1] #FR3 for unproductive
TRB_prod=TRBp[TRBp.aaSeqFR3.str.find('_')==-1] #Final productive beta

#Save productive dataframes into .txt files.
TRA_prod.to_csv('TRAproductive.txt', sep='\t', index=False)
TRB_prod.to_csv('TRBproductive.txt', sep='\t', index=False)

#Load productive and unproductive .txt files.
TRA_prod=pd.read_csv('TRAproductive.txt', sep='\t')
TRB_prod=pd.read_csv('TRBproductive.txt', sep='\t')

#Split the discription column into droplet barcode, molecular barcode and cell type.
DescA_prod=TRA_prod['descrR1'].str.split('_', 3, expand=True)
DescB_prod=TRB_prod['descrR1'].str.split('_', 3, expand=True)

#Select the columns with barcode, molecular and count information.
DescA_prod = DescA_prod[DescA_prod.columns.tolist()[1:4]]
DescB_prod = DescB_prod[DescB_prod.columns.tolist()[1:4]]

#Rename the columns to 'barcode', 'molecular', 'number'.
DescA_prod.columns=['droplet', 'molecular','number']
DescB_prod.columns=['droplet', 'molecular','number']

#Overwrite original dataframes with a separated Description column.
TRA_prod=pd.concat([DescA_prod,TRA_prod.drop('descrR1', 1)],axis=1)
TRB_prod=pd.concat([DescB_prod,TRB_prod.drop('descrR1', 1)],axis=1)

#We are now ready to start analysis.
print 'Ready for analysis'

#Get a list of droplets that appear only once.
DropletA_prod_count = pd.DataFrame(TRA_prod.droplet.value_counts())
UniqueA_prod_droplet = DropletA_prod_count[DropletA_prod_count['droplet']==1].index.tolist()

DropletB_prod_count = pd.DataFrame(TRB_prod.droplet.value_counts())
UniqueB_prod_droplet = DropletB_prod_count[DropletB_prod_count['droplet']==1].index.tolist()

#Create dataframe of droplets that appear multiple times: Copy original dataframe,
#and remove rows containing unique droplets.
print 'starting productive alpha chain'
TRAm_prod=TRA_prod
for i in UniqueA_prod_droplet:
    TRAm_prod = TRAm_prod[TRAm_prod['droplet'] != i]

print 'starting productive beta chain'
TRBm_prod=TRB_prod
for i in UniqueB_prod_droplet:
    TRBm_prod = TRBm_prod[TRBm_prod['droplet'] != i]

#Now extract index numbers from multiple droplet dataframe and drop those
#rows from original dataframe to create unique droplet dataframe.
TRA1_prod=TRA_prod.drop(TRA_prod.index[TRAm_prod.index.tolist()])
TRB1_prod=TRB_prod.drop(TRB_prod.index[TRBm_prod.index.tolist()])

#Some droplets appear more than once, but are associated with the same sequence.
#Next step is to find those droplets, collapse the sequences by adding the counts
#and add it to unique droplets dataframe.

#Reset indeces of multiple droplet dataframe.
TRAm_prod=TRAm_prod.reset_index(drop=True)
TRBm_prod=TRBm_prod.reset_index(drop=True)

#Get a list of droplets from multiple droplet dataframe
DropletA_prod_list = TRAm_prod.droplet.unique().tolist()
DropletB_prod_list = TRBm_prod.droplet.unique().tolist()

#Make a dataframe of droplets that are associated with different sequences.
print 'starting productive alpha chain 2'
TRAm_prod_final=TRAm_prod
for i in DropletA_prod_list:
    a=TRAm_prod[TRAm_prod['droplet']==i]
    if len(a.readSequence.value_counts())==1:
        TRAm_prod_final = TRAm_prod_final[TRAm_prod_final.droplet != i]

print 'starting productive beta chain 2'
TRBm_prod_final=TRBm_prod
for i in DropletB_prod_list:
    a=TRBm_prod[TRBm_prod['droplet']==i]
    if len(a.readSequence.value_counts())==1:
        TRBm_prod_final = TRBm_prod_final[TRBm_prod_final.droplet != i]

#Make a dataframe of droplets that have associated with identical sequences multiple times.
TRA1_prod_uncollapsed=TRAm_prod.drop(TRAm_prod.index[TRAm_prod_final.index.tolist()])
TRB1_prod_uncollapsed=TRBm_prod.drop(TRBm_prod.index[TRBm_prod_final.index.tolist()])

#Collapse identical sequences into one so that resulting dataframe has each droplet
#associated with one sequence.

#Make list of unique droplets.
TRA1_prod_droplet=TRA1_prod_uncollapsed.droplet.unique().tolist()
TRB1_prod_droplet=TRB1_prod_uncollapsed.droplet.unique().tolist()

#Make a collapsed dataframe to fill in with collapsed sequences.
TRA1_prod_collapsed=pd.DataFrame(index=range(0,len(TRA1_prod_droplet)),columns=TRA1_prod_uncollapsed.columns.tolist())
TRB1_prod_collapsed=pd.DataFrame(index=range(0,len(TRB1_prod_droplet)),columns=TRB1_prod_uncollapsed.columns.tolist())

#Convert number string into float.
TRA1_prod_uncollapsed.number=TRA1_prod_uncollapsed.number.astype(float)
TRB1_prod_uncollapsed.number=TRB1_prod_uncollapsed.number.astype(float)

#Fill the dataframes with collapsed sequences and added counts from number column.
print 'starting productive single alpha chain'
for i in range(len(TRA1_prod_droplet)):
    a=TRA1_prod_uncollapsed[TRA1_prod_uncollapsed['droplet'] == TRA1_prod_droplet[i]]
    TRA1_prod_collapsed.iloc[i]=a.iloc[0]
    TRA1_prod_collapsed.number.iloc[i]=sum(a.number)

print 'starting productive single beta chain'
for i in range(len(TRB1_prod_droplet)):
    a=TRB1_prod_uncollapsed[TRB1_prod_uncollapsed['droplet'] == TRB1_prod_droplet[i]]
    TRB1_prod_collapsed.iloc[i]=a.iloc[0]
    TRB1_prod_collapsed.number.iloc[i]=sum(a.number)

# Now there are two separate dataframes of droplets associated with a single sequence.
# and one dataframe of droplets associated with more than one sequence.

# Now concatenate the two unique droplet dataframes.
TRA_prod_single=pd.concat([TRA1_prod,TRA1_prod_collapsed],axis=0)
TRB_prod_single=pd.concat([TRB1_prod,TRB1_prod_collapsed],axis=0)

# Save these dataframes in a txt file.
print 'saving unique droplet files'
TRA_prod_single.to_csv('TRAsingle_productive.txt', sep='\t', index=False)
TRB_prod_single.to_csv('TRBsingle_productive.txt', sep='\t', index=False)

TRAm_prod_final.to_csv('TRAm_productive.txt', sep='\t', index=False)
TRBm_prod_final.to_csv('TRBm_productive.txt', sep='\t', index=False)

# Convert the remaining dataframe of droplets with muliple sequences into fasta file
# for clonality MIXCR analysis.
print 'making fasta files'
#FASTA file preparation:
#Insert > symbol for fasta file.
TRAm_prod_final.insert(0,'fasta','>')
TRBm_prod_final.insert(0,'fasta','>')

#Take only a sequence column
TRAm_prod_seq=TRAm_prod_final[['readSequence']]
TRBm_prod_seq=TRBm_prod_final[['readSequence']]

#Put all sequence info into one string.
TRAm_prod_info=pd.DataFrame(TRAm_prod_final.fasta+'_'+TRAm_prod_final.droplet+'_'+TRAm_prod_final.molecular+'_'+TRAm_prod_final.number, columns=['info'])
TRBm_prod_info=pd.DataFrame(TRBm_prod_final.fasta+'_'+TRBm_prod_final.droplet+'_'+TRBm_prod_final.molecular+'_'+TRBm_prod_final.number, columns=['info'])

#Concatenate info and sequence dataframes.
TRAm_prod_fasta=pd.concat([TRAm_prod_info,TRAm_prod_seq], axis=1)
TRBm_prod_fasta=pd.concat([TRBm_prod_info,TRBm_prod_seq], axis=1)

#Create fasta files so that >info of each sequence is the 1st row and each sequence is the 2nd row.
TRAm_prod_fasta.to_csv('TRAm_productive.fasta', sep='\n', index=False, header=False)
TRBm_prod_fasta.to_csv('TRBm_productive.fasta', sep='\n', index=False, header=False)
