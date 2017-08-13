"""

"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

TRp_1=pd.read_csv('S1/TRpairs_productive.txt', sep='\t')
TRp_2=pd.read_csv('S2/TRpairs_productive.txt', sep='\t')
TRp_3=pd.read_csv('S3/TRpairs_productive.txt', sep='\t')
TRp_4=pd.read_csv('S4/TRpairs_productive.txt', sep='\t')
TRp_5=pd.read_csv('S5/TRpairs_productive.txt', sep='\t')
TRp_all=pd.concat([TRp_1,TRp_2,TRp_3,TRp_4,TRp_5],axis=0)

TRp_list=[TRp_1,TRp_2,TRp_3,TRp_4,TRp_5]

#For pooled data.
df=TRp_all

beta_list=df.cloneB.unique().tolist()
counts_A=pd.DataFrame(index=beta_list,columns=['alpha_count'])
for i in range(len(counts_A)):
    counts_A.ix[counts_A.index[i]]=len(df[df.cloneB==counts_A.index[i]].cloneA.unique())
x_alpha=np.array(counts_A.alpha_count.value_counts().index)
y_alpha=np.array(counts_A.alpha_count.value_counts())

alpha_list=df.cloneA.unique().tolist()
counts_B=pd.DataFrame(index=alpha_list,columns=['beta_count'])
for i in range(len(counts_B)):
    counts_B.ix[counts_B.index[i]]=len(df[df.cloneA==counts_B.index[i]].cloneB.unique())
x_beta=np.array(counts_B.beta_count.value_counts().index)
y_beta=np.array(counts_B.beta_count.value_counts())

y_alpha_freq=np.divide(y_alpha,np.sum(y_alpha)+0.0)
y_beta_freq=np.divide(y_beta,np.sum(y_beta)+0.0)

f, ax = plt.subplots()
ax.plot(x_alpha[0:-1],y_alpha_freq[0:-1],'o',markersize=10,color='yellow',markeredgecolor='black',label='# alphas')
ax.plot(x_beta[0:-1],y_beta_freq[0:-1],'o',markersize=10,color='green',markeredgecolor='black',label='# betas')
ax.legend(loc='upper right',fontsize=12)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('chain count',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.tight_layout()
plt.savefig('unique_chain_count_all.png')


#For individual subjects.
for i in range(5):
    df=TRp_list[i]
    beta_list=df.cloneB.unique().tolist()
    counts_A=pd.DataFrame(index=beta_list,columns=['alpha_count'])
    for i in range(len(counts_A)):
        counts_A.ix[counts_A.index[i]]=len(df[df.cloneB==counts_A.index[i]].cloneA.unique())
    x_alpha=np.array(counts_A.alpha_count.value_counts().index)
    y_alpha=np.array(counts_A.alpha_count.value_counts())

    alpha_list=df.cloneA.unique().tolist()
    counts_B=pd.DataFrame(index=alpha_list,columns=['beta_count'])
    for i in range(len(counts_B)):
        counts_B.ix[counts_B.index[i]]=len(df[df.cloneA==counts_B.index[i]].cloneB.unique())
    x_beta=np.array(counts_B.beta_count.value_counts().index)
    y_beta=np.array(counts_B.beta_count.value_counts())

    y_alpha_freq=np.divide(y_alpha,np.sum(y_alpha)+0.0)
    y_beta_freq=np.divide(y_beta,np.sum(y_beta)+0.0)

    f, ax = plt.subplots()
    ax.plot(x_alpha[0:-1],y_alpha_freq[0:-1],'o',markersize=10,color='yellow',markeredgecolor='black',label='# alphas')
    ax.plot(x_beta[0:-1],y_beta_freq[0:-1],'o',markersize=10,color='green',markeredgecolor='black',label='# betas')
    ax.legend(loc='upper right',fontsize=12)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('chain count',fontsize=15)
    plt.ylabel('frequency',fontsize=15)
    plt.tight_layout()
    plt.savefig('unique_chain_count_'+str(i)+'.png')
