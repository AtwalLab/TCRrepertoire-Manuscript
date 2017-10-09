"""
Heatmaps of gene usage as well as CDR3 lengh frequencies:
p(state1,state2)/sample size
"""

from __future__ import division
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math

#Define lists for each gene strings.
TRAV=['TRAV1-1*00','TRAV1-2*00','TRAV2*00','TRAV3*00','TRAV4*00','TRAV5*00',
      'TRAV6*00','TRAV7*00','TRAV8-1*00','TRAV9-1*00','TRAV10*00','TRAV11*00',
      'TRAV12-1*00','TRAV8-2*00','TRAV8-3*00','TRAV13-1*00','TRAV12-2*00',
      'TRAV8-4*00','TRAV8-5*00','TRAV13-2*00','TRAV14DV4*00','TRAV9-2*00',
      'TRAV15*00','TRAV12-3*00','TRAV8-6*00','TRAV16*00','TRAV17*00','TRAV18*00',
      'TRAV19*00','TRAV20*00','TRAV21*00','TRAV22*00','TRAV23DV6*00','TRAV24*00',
      'TRAV25*00','TRAV26-1*00','TRAV8-7*00','TRAV27*00','TRAV28*00','TRAV29DV5*00',
      'TRAV30*00','TRAV31*00','TRAV32*00','TRAV33*00','TRAV26-2*00','TRAV34*00',
      'TRAV35*00','TRAV36DV7*00','TRAV37*00','TRAV38-1*00','TRAV38-2DV8*00',
      'TRAV39*00','TRAV40*00','TRAV41*00']

TRBV=['TRBV1*00','TRBV2*00','TRBV3-1*00','TRBV4-1*00','TRBV5-1*00','TRBV6-1*00',
      'TRBV7-1*00','TRBV4-2*00','TRBV6-2*00','TRBV3-2*00','TRBV4-3*00',
      'TRBV6-3*00','TRBV7-2*00','TRBV8-1*00','TRBV5-2*00','TRBV6-4*00',
      'TRBV7-3*00','TRBV8-2*00','TRBV5-3*00','TRBV9*00','TRBV10-1*00',
      'TRBV11-1*00','TRBV12-1*00','TRBV10-2*00','TRBV11-2*00','TRBV12-2*00',
      'TRBV6-5*00','TRBV7-4*00','TRBV5-4*00','TRBV6-6*00','TRBV7-5*00',
      'TRBV5-5*00','TRBV6-7*00','TRBV7-6*00','TRBV5-6*00','TRBV6-8*00',
      'TRBV7-7*00','TRBV5-7*00','TRBV6-9*00','TRBV7-8*00','TRBV5-8*00',
      'TRBV7-9*00','TRBV13*00','TRBV10-3*00','TRBV11-3*00','TRBV12-3*00',
      'TRBV12-4*00','TRBV12-5*00','TRBV14*00','TRBV15*00','TRBV16*00',
      'TRBV17*00','TRBV18*00','TRBV19*00','TRBV20-1*00','TRBV21-1*00',
      'TRBV22-1*00','TRBV23-1*00','TRBV24-1*00','TRBV25-1*00','TRBV26*00',
      'TRBV27*00','TRBV28*00','TRBV29-1*00','TRBV30*00']

TRBJ=['TRBJ1-1*00','TRBJ1-2*00','TRBJ1-3*00','TRBJ1-4*00','TRBJ1-5*00',
    'TRBJ1-6*00','TRBJ2-1*00','TRBJ2-2*00','TRBJ2-3*00','TRBJ2-4*00',
    'TRBJ2-5*00','TRBJ2-6*00','TRBJ2-7*00']

TRAJ=['TRAJ61*00','TRAJ60*00','TRAJ59*00','TRAJ58*00','TRAJ57*00','TRAJ56*00',
      'TRAJ55*00','TRAJ54*00','TRAJ53*00','TRAJ52*00','TRAJ51*00','TRAJ50*00',
      'TRAJ49*00','TRAJ48*00','TRAJ47*00','TRAJ46*00','TRAJ45*00','TRAJ44*00',
      'TRAJ43*00','TRAJ42*00','TRAJ41*00','TRAJ40*00','TRAJ39*00','TRAJ38*00',
      'TRAJ37*00','TRAJ36*00','TRAJ35*00','TRAJ34*00','TRAJ33*00','TRAJ32*00',
      'TRAJ31*00','TRAJ30*00','TRAJ29*00','TRAJ28*00','TRAJ27*00','TRAJ26*00',
      'TRAJ25*00','TRAJ24*00','TRAJ23*00','TRAJ22*00','TRAJ21*00','TRAJ20*00',
      'TRAJ19*00','TRAJ18*00','TRAJ17*00','TRAJ16*00','TRAJ15*00','TRAJ14*00',
      'TRAJ13*00','TRAJ12*00','TRAJ11*00','TRAJ10*00','TRAJ9*00','TRAJ8*00',
      'TRAJ7*00','TRAJ6*00','TRAJ5*00','TRAJ4*00','TRAJ3*00','TRAJ2*00','TRAJ1*00']

TRBD = ['TRBD1*00','TRBD2*00','nan']

def frequency_heatmap(df, gene1, gene2, gene1_list, gene2_list, fig_title=''):
    """This function calculates genes usage frequences: frequency=p(gene1,gene2)/sample size
    Inputs: df - dataframe, gene1, gene2 - genes of interest,
    gene1_list, gene2_list - list of all genes, fig_title - title for saved figure.
    Returns: saved heatmap figure."""

    #create a joint gene column
    gene1_2=gene1+'_'+gene2
    df[gene1_2]=df[gene1].astype(str)+'_'+df[gene2].astype(str)
    #create dataframes of counted genes
    g1=pd.DataFrame(df[gene1].astype(str).value_counts())
    g2=pd.DataFrame(df[gene2].astype(str).value_counts())
    g1_2=pd.DataFrame(df[gene1_2].astype(str).value_counts())
    #Make a list of all possibilities of joint genes.
    DJ_list=[]
    for i in gene2_list:
        for j in gene1_list:
            DJ_list.append(j+'_'+i)
    #Make and empty dataframe for gene counts.
    counts=pd.DataFrame(index=DJ_list, columns = [gene1_2,gene1,gene2])
    #Fill in the dataframe. If the joint genes is not in the dataset, fill in with 0.
    for i in range(len(counts)):
        if counts.index[i] in g1_2.index.tolist():
            genes=counts.index[i].split('_',1)
            counts[gene1_2].iloc[i]=g1_2[gene1_2].ix[counts.index[i]]
            counts[gene1].iloc[i]=g1[gene1].ix[genes[0]]
            counts[gene2].iloc[i]=g2[gene2].ix[genes[1]]
        else:
            counts[gene1_2].iloc[i]=0
            counts[gene1].iloc[i]=0
            counts[gene2].iloc[i]=0
    #Create a dataframe for frequences.
    frequencies=pd.DataFrame(index=gene1_list, columns=gene2_list)
    #Fill in the dataframe with frequency=p(gene1,gene2)/sample size
    for i in gene2_list:
        for j in gene1_list:
            genes=j+'_'+i
            if genes in g1_2.index.tolist():
                frequencies.ix[j,i]=(counts.ix[genes][gene1_2]/len(df))
            else:
                frequencies.ix[j,i]=0
    #Make a nunpy array of frequencies converting it into float.
    freq=np.array(frequencies).astype(float)
    #Make a heatmap.
    plt.figure(figsize=(10,10))
    plt.pcolor(freq,cmap=plt.cm.seismic)
    plt.yticks(np.arange(0.5, len(frequencies.index), 1), frequencies.index)
    plt.xticks(np.arange(0.5, len(frequencies.columns), 1), frequencies.columns, rotation=90)
    plt.colorbar()
    plt.tight_layout()

    return plt.savefig(fig_title+'.png')

#Modified function for CDR3 lengh usage.
def frequency_heatmap_CDR3(df, gene1, gene2, gene1_list, gene2_list, fig_title=''):
    """This function calculates CDR3 lenght usage frequences: frequency=p(gene1,gene2)/sample size
    Inputs: df - dataframe, gene1, gene2 - CDR3 lenghts of interest,
    gene1_list, gene2_list - list of all CDR3 lenghts, fig_title - title for saved figure.
    Returns: saved heatmap figure."""

    #create a joint gene column
    gene1_2=gene1+'_'+gene2
    df[gene1_2]=df[gene1].astype(str)+'_'+df[gene2].astype(str)
    #create dataframes of counted genes
    g1=pd.DataFrame(df[gene1].value_counts())
    g2=pd.DataFrame(df[gene2].value_counts())
    g1_2=pd.DataFrame(df[gene1_2].value_counts())
    #Make a list of all possibilities of joint genes.
    DJ_list=[]
    for i in gene2_list:
        for j in gene1_list:
            DJ_list.append(str(j)+'_'+str(i))
    #Make and empty dataframe for gene counts.
    counts=pd.DataFrame(index=DJ_list, columns = [gene1_2,gene1,gene2])
    #Fill in the dataframe. I the joint genes is not in the dataset, fill in with 0.
    for i in range(len(counts)):
        if counts.index[i] in g1_2.index.tolist():
            genes=counts.index[i].split('_',1)
            counts[gene1_2].iloc[i]=g1_2[gene1_2].ix[counts.index[i]]
            counts[gene1].iloc[i]=g1[gene1].ix[float(genes[0])]
            counts[gene2].iloc[i]=g2[gene2].ix[float(genes[1])]
        else:
            counts[gene1_2].iloc[i]=0
            counts[gene1].iloc[i]=0
            counts[gene2].iloc[i]=0
    #Create a dataframe for frequences.
    frequencies=pd.DataFrame(index=gene1_list, columns=gene2_list)
    #Fill in the dataframe with frequency=p(gene1,gene2)/sample size
    for i in gene2_list:
        for j in gene1_list:
            genes=str(j)+'_'+str(i)
            if genes in g1_2.index.tolist():
                frequencies.ix[j,i]=(counts.ix[genes][gene1_2]/len(df))
            else:
                frequencies.ix[j,i]=0
    #Make a nunpy array of frequencies converting it into float.
    freq=np.array(frequencies).astype(float)
    #Make a heatmap.
    plt.figure(figsize=(10,10))
    plt.pcolor(freq,cmap=plt.cm.seismic)
    plt.yticks(np.arange(0.5, len(frequencies.index), 1), frequencies.index)
    plt.xticks(np.arange(0.5, len(frequencies.columns), 1), frequencies.columns, rotation=90)
    plt.colorbar()
    plt.tight_layout()

    return plt.savefig(fig_title+'.png')

#Load data
TRp_1=pd.read_csv('S1/TRpairs_productive.txt', sep='\t')
TRp_2=pd.read_csv('S2/TRpairs_productive.txt', sep='\t')
TRp_3=pd.read_csv('S3/TRpairs_productive.txt', sep='\t')
TRp_4=pd.read_csv('S4/TRpairs_productive.txt', sep='\t')
TRp_5=pd.read_csv('S5/TRpairs_productive.txt', sep='\t')

#Concatenate to pool all subjects
TRp_all=pd.concat([TRp_1,TRp_2,TRp_3,TRp_4,TRp_5],axis=0)

#SUPPLEMENTARY FIGURE 4
#Make heatmaps for each gene pair across the chains
frequency_heatmap(TRp_1,'BestVhitA','BestVhitB', TRAV, TRBV, fig_title='aVbV_heatmap_S1')
frequency_heatmap(TRp_1,'BestVhitA','BestJhitB', TRAV, TRBJ, fig_title='aVbJ_heatmap_S1')
frequency_heatmap(TRp_1,'BestJhitA','BestVhitB', TRAJ, TRBV, fig_title='aJbV_heatmap_S1')
frequency_heatmap(TRp_1,'BestJhitA','BestJhitB', TRAJ, TRBJ, fig_title='aJbJ_heatmap_S1')
frequency_heatmap(TRp_1,'BestVhitA','BestDhitB', TRAV, TRBD, fig_title='aVbD_heatmap_S1')
frequency_heatmap(TRp_1,'BestJhitA','BestDhitB', TRAJ, TRBD, fig_title='aJbD_heatmap_S1')

frequency_heatmap(TRp_2,'BestVhitA','BestVhitB', TRAV, TRBV, fig_title='aVbV_heatmap_S2')
frequency_heatmap(TRp_2,'BestVhitA','BestJhitB', TRAV, TRBJ, fig_title='aVbJ_heatmap_S2')
frequency_heatmap(TRp_2,'BestJhitA','BestVhitB', TRAJ, TRBV, fig_title='aJbV_heatmap_S2')
frequency_heatmap(TRp_2,'BestJhitA','BestJhitB', TRAJ, TRBJ, fig_title='aJbJ_heatmap_S2')
frequency_heatmap(TRp_2,'BestVhitA','BestDhitB', TRAV, TRBD, fig_title='aVbD_heatmap_S2')
frequency_heatmap(TRp_2,'BestJhitA','BestDhitB', TRAJ, TRBD, fig_title='aJbD_heatmap_S2')

frequency_heatmap(TRp_3,'BestVhitA','BestVhitB', TRAV, TRBV, fig_title='aVbV_heatmap_S3')
frequency_heatmap(TRp_3,'BestVhitA','BestJhitB', TRAV, TRBJ, fig_title='aVbJ_heatmap_S3')
frequency_heatmap(TRp_3,'BestJhitA','BestVhitB', TRAJ, TRBV, fig_title='aJbV_heatmap_S3')
frequency_heatmap(TRp_3,'BestJhitA','BestJhitB', TRAJ, TRBJ, fig_title='aJbJ_heatmap_S3')
frequency_heatmap(TRp_3,'BestVhitA','BestDhitB', TRAV, TRBD, fig_title='aVbD_heatmap_S3')
frequency_heatmap(TRp_3,'BestJhitA','BestDhitB', TRAJ, TRBD, fig_title='aJbD_heatmap_S3')

frequency_heatmap(TRp_4,'BestVhitA','BestVhitB', TRAV, TRBV, fig_title='aVbV_heatmap_S4')
frequency_heatmap(TRp_4,'BestVhitA','BestJhitB', TRAV, TRBJ, fig_title='aVbJ_heatmap_S4')
frequency_heatmap(TRp_4,'BestJhitA','BestVhitB', TRAJ, TRBV, fig_title='aJbV_heatmap_S4')
frequency_heatmap(TRp_4,'BestJhitA','BestJhitB', TRAJ, TRBJ, fig_title='aJbJ_heatmap_S4')
frequency_heatmap(TRp_4,'BestVhitA','BestDhitB', TRAV, TRBD, fig_title='aVbD_heatmap_S4')
frequency_heatmap(TRp_4,'BestJhitA','BestDhitB', TRAJ, TRBD, fig_title='aJbD_heatmap_S4')

frequency_heatmap(TRp_5,'BestVhitA','BestVhitB', TRAV, TRBV, fig_title='aVbV_heatmap_S5')
frequency_heatmap(TRp_5,'BestVhitA','BestJhitB', TRAV, TRBJ, fig_title='aVbJ_heatmap_S5')
frequency_heatmap(TRp_5,'BestJhitA','BestVhitB', TRAJ, TRBV, fig_title='aJbV_heatmap_S5')
frequency_heatmap(TRp_5,'BestJhitA','BestJhitB', TRAJ, TRBJ, fig_title='aJbJ_heatmap_S5')
frequency_heatmap(TRp_5,'BestVhitA','BestDhitB', TRAV, TRBD, fig_title='aVbD_heatmap_S5')
frequency_heatmap(TRp_5,'BestJhitA','BestDhitB', TRAJ, TRBD, fig_title='aJbD_heatmap_S5')

#FIGURES 3B-G
frequency_heatmap(TRp_all,'BestVhitA','BestVhitB', TRAV, TRBV, fig_title='aVbV_heatmap_all')
frequency_heatmap(TRp_all,'BestVhitA','BestJhitB', TRAV, TRBJ, fig_title='aVbJ_heatmap_all')
frequency_heatmap(TRp_all,'BestJhitA','BestVhitB', TRAJ, TRBV, fig_title='aJbV_heatmap_all')
frequency_heatmap(TRp_all,'BestJhitA','BestJhitB', TRAJ, TRBJ, fig_title='aJbJ_heatmap_all')
frequency_heatmap(TRp_all,'BestVhitA','BestDhitB', TRAV, TRBD, fig_title='aVbD_heatmap_all')
frequency_heatmap(TRp_all,'BestJhitA','BestDhitB', TRAJ, TRBD, fig_title='aJbD_heatmap_all')

#Make CDR3 length lists for each dataframe for alpha and beta chains.
A_lengths_list1=sorted(TRp_1.LengthCDR3A.unique().tolist())
B_lengths_list1=sorted(TRp_1.LengthCDR3B.unique().tolist())
A_lengths_list2=sorted(TRp_2.LengthCDR3A.unique().tolist())
B_lengths_list2=sorted(TRp_2.LengthCDR3B.unique().tolist())
A_lengths_list3=sorted(TRp_3.LengthCDR3A.unique().tolist())
B_lengths_list3=sorted(TRp_3.LengthCDR3B.unique().tolist())
A_lengths_list4=sorted(TRp_4.LengthCDR3A.unique().tolist())
B_lengths_list4=sorted(TRp_4.LengthCDR3B.unique().tolist())
A_lengths_list5=sorted(TRp_5.LengthCDR3A.unique().tolist())
B_lengths_list5=sorted(TRp_5.LengthCDR3B.unique().tolist())
A_lengths_list_all=sorted(TRp_all.LengthCDR3A.unique().tolist())
B_lengths_list_all=sorted(TRp_all.LengthCDR3B.unique().tolist())

#SUPPLEMENTARY FIGURE 5
#Make heatmaps for CDR3 lenghts.
frequency_heatmap_CDR3(TRp_1,'LengthCDR3A','LengthCDR3B', A_lengths_list1, B_lengths_list1, fig_title='CDR3_length_heatmap_1')
frequency_heatmap_CDR3(TRp_2,'LengthCDR3A','LengthCDR3B', A_lengths_list2, B_lengths_list2, fig_title='CDR3_length_heatmap_2')
frequency_heatmap_CDR3(TRp_3,'LengthCDR3A','LengthCDR3B', A_lengths_list3, B_lengths_list3, fig_title='CDR3_length_heatmap_3')
frequency_heatmap_CDR3(TRp_4,'LengthCDR3A','LengthCDR3B', A_lengths_list4, B_lengths_list4, fig_title='CDR3_length_heatmap_4')
frequency_heatmap_CDR3(TRp_5,'LengthCDR3A','LengthCDR3B', A_lengths_list5, B_lengths_list5, fig_title='CDR3_length_heatmap_5')

#FIGURE 3A
frequency_heatmap_CDR3(TRp_all,'LengthCDR3A','LengthCDR3B', A_lengths_list_all, B_lengths_list_all, fig_title='CDR3_length_heatmap_all')
