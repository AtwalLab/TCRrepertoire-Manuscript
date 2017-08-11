"""

FILTERING. STEP 1.
This code takes the raw sequencing fasta files and prepares for aligning with MIXCR.
Removes all sequences with count=1.

REQUIREMENTS BEFORE RUNNING:
1. Go to location where the files are.
2. Edit the code adjusting for each subject when raw file format is different.
The script is set for Subjects 1-4.

!For Subject 5, change LINE 43 into: ID_df=pd.DataFrame(ID,columns=['droplet','molecular','number','celltype'])
!Also, change celltype (LINES 50,51) to: 'PRCONS=P5trPRC2-TRA' and 'PRCONS=P5trPRC2-TRB'.

4. Make sure the string for count is the same as the default in the script, LINE 54,55.

"""

from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import stats
import pandas as pd
from Bio import SeqIO
import glob

file_list=glob.glob('*.fasta') #grab all fasta files

#initialize sequence and sequence ID lists to fill from raw file.
ID=[]
sequence=[]

for file_name in file_list:
    handle = open(file_name)
    for record in SeqIO.parse(handle, "fasta"):
        a=sum([i.split('|') for i in [record.id]],[])
        ID.append(sum([i.split('_') for i in a],[]))
        sequence.append([record.seq])
    handle.close()

#concatenate ID and sequence lists into one dataframe.
ID_df=pd.DataFrame(ID,columns=['droplet','molecular','number','celltype','orig','count','assembly'])
ID_df=ID_df[['droplet','molecular','number','celltype']]

sequence_df=pd.DataFrame(sequence,columns=['sequence'])
df = pd.concat([ID_df, sequence_df], axis=1)

#call TRA/TRB cell types. Make sure to adjust the strings by looking at the raw file.
TRA=df[df['celltype']=='PRCONS=p5trpcr2-tra'].reset_index(drop=True)
TRB=df[df['celltype']=='PRCONS=p5trpcr2-trb'].reset_index(drop=True)

#Remove sequences with count of 1, and extract the count integer from the string.
TRA = TRA[TRA.number != 'CONSCOUNT=1']
TRB = TRB[TRB.number != 'CONSCOUNT=1']
TRA.number=TRA.number.str.extract('(\d+)')
TRB.number=TRB.number.str.extract('(\d+)')

#FASTA file preparation:
#Insert > symbol for fasta file.
TRA.insert(0,'fasta','>')
TRB.insert(0,'fasta','>')

#Take only a sequence column
TRA_seq=TRA[['sequence']]
TRB_seq=TRB[['sequence']]

#Put all sequence info into one string.
TRA_info=pd.DataFrame(TRA.fasta+'_'+TRA.droplet+'_'+TRA.molecular+'_'+TRA.number, columns=['info'])
TRB_info=pd.DataFrame(TRB.fasta+'_'+TRB.droplet+'_'+TRB.molecular+'_'+TRB.number, columns=['info'])

#Concatenate info and sequence dataframes.
TRAfasta=pd.concat([TRA_info,TRA_seq], axis=1)
TRBfasta=pd.concat([TRB_info,TRB_seq], axis=1)

#Create fasta files so that >info of each sequence is the 1st row and each sequence is the 2nd row.
TRAfasta.to_csv('TRAtoalign.fasta', sep='\n', index=False, header=False)
TRBfasta.to_csv('TRBtoalign.fasta', sep='\n', index=False, header=False)
