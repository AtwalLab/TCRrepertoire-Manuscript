"""
This script inputs single alpha and beta chain dataframes divided by CD4 and CD8.
Plots clone size distributions of CD4 and CD8 single chain alpha and beta repertoires
of each subject.
"""

from __future__ import division
from scipy.optimize import newton
from scipy.special import zeta
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


def powerlaw_single_toplot(df):
    """This function calculates clone sizes and their frequencies.
    Input: pandas dataframe of paired sequences.
    Returns: numpy array of clone sizes, numpy array of clone size frequencies."""

    clone_size_count=np.array(pd.DataFrame(df.nSeqCDR3.value_counts()).nSeqCDR3.value_counts())
    clone_size=np.array(pd.DataFrame(df.nSeqCDR3.value_counts()).nSeqCDR3.value_counts().index)
    clone_size_frequency=np.divide(clone_size_count, np.sum(clone_size_count).astype(float))

    return clone_size, clone_size_frequency

#Load data
TRA_cd4_1=pd.read_csv('S1/TRAcd4_productive.txt', sep='\t')
TRA_cd4_2=pd.read_csv('S2/TRAcd4_productive.txt', sep='\t')
TRA_cd4_3=pd.read_csv('S3/TRAcd4_productive.txt', sep='\t')
TRA_cd4_4=pd.read_csv('S4/TRAcd4_productive.txt', sep='\t')
TRA_cd4_5=pd.read_csv('S5/TRAcd4_productive.txt', sep='\t')

TRA_cd8_1=pd.read_csv('S1/TRAcd8_productive.txt', sep='\t')
TRA_cd8_2=pd.read_csv('S2/TRAcd8_productive.txt', sep='\t')
TRA_cd8_3=pd.read_csv('S3/TRAcd8_productive.txt', sep='\t')
TRA_cd8_4=pd.read_csv('S4/TRAcd8_productive.txt', sep='\t')
TRA_cd8_5=pd.read_csv('S5/TRAcd8_productive.txt', sep='\t')

TRB_cd4_1=pd.read_csv('S1/TRBcd4_productive.txt', sep='\t')
TRB_cd4_2=pd.read_csv('S2/TRBcd4_productive.txt', sep='\t')
TRB_cd4_3=pd.read_csv('S3/TRBcd4_productive.txt', sep='\t')
TRB_cd4_4=pd.read_csv('S4/TRBcd4_productive.txt', sep='\t')
TRB_cd4_5=pd.read_csv('S5/TRBcd4_productive.txt', sep='\t')

TRB_cd8_1=pd.read_csv('S1/TRBcd8_productive.txt', sep='\t')
TRB_cd8_2=pd.read_csv('S2/TRBcd8_productive.txt', sep='\t')
TRB_cd8_3=pd.read_csv('S3/TRBcd8_productive.txt', sep='\t')
TRB_cd8_4=pd.read_csv('S4/TRBcd8_productive.txt', sep='\t')
TRB_cd8_5=pd.read_csv('S5/TRBcd8_productive.txt', sep='\t')

#Get clone sizes and their frequencies for plotting.
CSA_cd4_prod_1, CSA_cd4_freq_prod_1=powerlaw_single_toplot(TRA_cd4_1)
CSA_cd4_prod_2, CSA_cd4_freq_prod_2=powerlaw_single_toplot(TRA_cd4_2)
CSA_cd4_prod_3, CSA_cd4_freq_prod_3=powerlaw_single_toplot(TRA_cd4_3)
CSA_cd4_prod_4, CSA_cd4_freq_prod_4=powerlaw_single_toplot(TRA_cd4_4)
CSA_cd4_prod_5, CSA_cd4_freq_prod_5=powerlaw_single_toplot(TRA_cd4_5)

CSA_cd8_prod_1, CSA_cd8_freq_prod_1=powerlaw_single_toplot(TRA_cd8_1)
CSA_cd8_prod_2, CSA_cd8_freq_prod_2=powerlaw_single_toplot(TRA_cd8_2)
CSA_cd8_prod_3, CSA_cd8_freq_prod_3=powerlaw_single_toplot(TRA_cd8_3)
CSA_cd8_prod_4, CSA_cd8_freq_prod_4=powerlaw_single_toplot(TRA_cd8_4)
CSA_cd8_prod_5, CSA_cd8_freq_prod_5=powerlaw_single_toplot(TRA_cd8_5)

CSB_cd4_prod_1, CSB_cd4_freq_prod_1=powerlaw_single_toplot(TRB_cd4_1)
CSB_cd4_prod_2, CSB_cd4_freq_prod_2=powerlaw_single_toplot(TRB_cd4_2)
CSB_cd4_prod_3, CSB_cd4_freq_prod_3=powerlaw_single_toplot(TRB_cd4_3)
CSB_cd4_prod_4, CSB_cd4_freq_prod_4=powerlaw_single_toplot(TRB_cd4_4)
CSB_cd4_prod_5, CSB_cd4_freq_prod_5=powerlaw_single_toplot(TRB_cd4_5)

CSB_cd8_prod_1, CSB_cd8_freq_prod_1=powerlaw_single_toplot(TRB_cd8_1)
CSB_cd8_prod_2, CSB_cd8_freq_prod_2=powerlaw_single_toplot(TRB_cd8_2)
CSB_cd8_prod_3, CSB_cd8_freq_prod_3=powerlaw_single_toplot(TRB_cd8_3)
CSB_cd8_prod_4, CSB_cd8_freq_prod_4=powerlaw_single_toplot(TRB_cd8_4)
CSB_cd8_prod_5, CSB_cd8_freq_prod_5=powerlaw_single_toplot(TRB_cd8_5)

#SUPPLEMENTARY FIGURE 7A
#CD4 alpha
#Subject 1
plt.plot(CSA_cd4_prod_1, CSA_cd4_freq_prod_1, 'o', color='red', markeredgecolor='none', markersize=5.0,alpha=1,label='S1')
#Subject 2
plt.plot(CSA_cd4_prod_2, CSA_cd4_freq_prod_2, 'o', color='orange', markeredgecolor='none', markersize=5.0,alpha=1,label='S2')
#Subject 3
plt.plot(CSA_cd4_prod_3, CSA_cd4_freq_prod_3, 'o', color='black', markeredgecolor='none', markersize=5.0,alpha=1,label='S3')
#Subject 4
plt.plot(CSA_cd4_prod_4, CSA_cd4_freq_prod_4, 'o', color='green', markeredgecolor='none', markersize=5.0,alpha=1,label='S4')
#Subject 5
plt.plot(CSA_cd4_prod_5, CSA_cd4_freq_prod_5, 'o', color='blue', markeredgecolor='none', markersize=5.0,alpha=1,label='S5')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**4])
plt.ylim([10**-5,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('alpha_CD4_powerlaws_paper.png')

#SUPPLEMENTARY FIGURE 7B
#CD4 beta
#Subject 1
plt.plot(CSB_cd4_prod_1, CSB_cd4_freq_prod_1, 'o', color='red', markeredgecolor='none', markersize=5.0,alpha=1,label='S1')
#Subject 2
plt.plot(CSB_cd4_prod_2, CSB_cd4_freq_prod_2, 'o', color='orange', markeredgecolor='none', markersize=5.0,alpha=1,label='S2')
#Subject 3
plt.plot(CSB_cd4_prod_3, CSB_cd4_freq_prod_3, 'o', color='black', markeredgecolor='none', markersize=5.0,alpha=1,label='S3')
#Subject 4
plt.plot(CSB_cd4_prod_4, CSB_cd4_freq_prod_4, 'o', color='green', markeredgecolor='none', markersize=5.0,alpha=1,label='S4')
#Subject 5
plt.plot(CSB_cd4_prod_5, CSB_cd4_freq_prod_5, 'o', color='blue', markeredgecolor='none', markersize=5.0,alpha=1,label='S5')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**4])
plt.ylim([10**-5,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('beta_CD4_powerlaws_paper.png')

#SUPPLEMENTARY FIGURE 7C
#CD8 alpha
#Subject 1
plt.plot(CSA_cd8_prod_1, CSA_cd8_freq_prod_1, 'o', color='red', markeredgecolor='none', markersize=5.0,alpha=1,label='S1')
#Subject 2
plt.plot(CSA_cd8_prod_2, CSA_cd8_freq_prod_2, 'o', color='orange', markeredgecolor='none', markersize=5.0,alpha=1,label='S2')
#Subject 3
plt.plot(CSA_cd8_prod_3, CSA_cd8_freq_prod_3, 'o', color='black', markeredgecolor='none', markersize=5.0,alpha=1,label='S3')
#Subject 4
plt.plot(CSA_cd8_prod_4, CSA_cd8_freq_prod_4, 'o', color='green', markeredgecolor='none', markersize=5.0,alpha=1,label='S4')
#Subject 5
plt.plot(CSA_cd8_prod_5, CSA_cd8_freq_prod_5, 'o', color='blue', markeredgecolor='none', markersize=5.0,alpha=1,label='S5')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**4])
plt.ylim([10**-5,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('alpha_CD8_powerlaws_paper.png')

#SUPPLEMENTARY FIGURE 7D
#CD8 beta
#Subject 1
plt.plot(CSB_cd8_prod_1, CSB_cd8_freq_prod_1, 'o', color='red', markeredgecolor='none', markersize=5.0,alpha=1,label='S1')
#Subject 2
plt.plot(CSB_cd8_prod_2, CSB_cd8_freq_prod_2, 'o', color='orange', markeredgecolor='none', markersize=5.0,alpha=1,label='S2')
#Subject 3
plt.plot(CSB_cd8_prod_3, CSB_cd8_freq_prod_3, 'o', color='black', markeredgecolor='none', markersize=5.0,alpha=1,label='S3')
#Subject 4
plt.plot(CSB_cd8_prod_4, CSB_cd8_freq_prod_4, 'o', color='green', markeredgecolor='none', markersize=5.0,alpha=1,label='S4')
#Subject 5
plt.plot(CSB_cd8_prod_5, CSB_cd8_freq_prod_5, 'o', color='blue', markeredgecolor='none', markersize=5.0,alpha=1,label='S5')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**4])
plt.ylim([10**-5,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('beta_CD8_powerlaws_paper.png')
