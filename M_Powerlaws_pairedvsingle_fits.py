"""
Plot powerlaws and fits, as well as exponent plots for comparing paired vs single
repertoires.
"""

from __future__ import division
from scipy.optimize import newton
from scipy.special import zeta
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#Defining the functions.
def powerlawpdf(x,alpha):
    """The power law probability function"""
    return (x**(-alpha))/Z(alpha,xmin,xmax)

def Z(alpha,xmin=1,xmax=np.infty):
    """The normalization function Z.  Note that default arguments of xmin and xmax
    make Z equivalent to Riemann zeta function which is already implemented in Scipy as
    zeta(alpha,1)"""
    if np.isfinite(xmax):
        s=0
        for i in xrange(xmin,xmax+1):
            s+=(1/(i**alpha))
    else:
        s=zeta(alpha,xmin)
    return s

def F(alpha):
    """The optimization function F(alpha). C is the second term in the definition
    of F(alpha) and is independent of alpha"""
    h = 1e-8
    Z_prime = (Z(alpha+h,xmin,xmax) - Z(alpha-h,xmin,xmax))/(2*h)
    return (Z_prime/Z(alpha,xmin,xmax))+C

def powerlaw_pvs_tofit(df):
    """This function calculates clone counts.
    Input: pandas dataframe of paired sequences.
    Returns: clone counts."""

    df['joint']=df['NSeqCDR3A']+'_'+df['NSeqCDR3B']

    paired=np.array(pd.DataFrame(df.joint.value_counts())).flatten()
    alpha=np.array(pd.DataFrame(df.NSeqCDR3A.value_counts())).flatten()
    beta=np.array(pd.DataFrame(df.NSeqCDR3B.value_counts())).flatten()

    return paired, alpha, beta

def powerlaw_paired_toplot(df):
    """This function calculates clone sizes and their frequencies.
    Input: pandas dataframe of paired sequences.
    Returns: numpy array of clone sizes, numpy array of clone size frequencies."""

    df['joint']=df['NSeqCDR3A']+'_'+df['NSeqCDR3B']

    clone_size_count=np.array(pd.DataFrame(df.joint.value_counts()).joint.value_counts())
    clone_size=np.array(pd.DataFrame(df.joint.value_counts()).joint.value_counts().index)
    clone_size_frequency=np.divide(clone_size_count, np.sum(clone_size_count).astype(float))

    return clone_size, clone_size_frequency

def powerlaw_pvs_toplot(df,chain):
    """This function calculates clone sizes and their frequencies.
    Input: pandas dataframe of single sequences.
    Returns: numpy array of clone sizes, numpy array of clone size frequencies."""

    clone_size_count=np.array(pd.DataFrame(df[chain].value_counts())[chain].value_counts())
    clone_size=np.array(pd.DataFrame(df[chain].value_counts())[chain].value_counts().index)
    clone_size_frequency=np.divide(clone_size_count, np.sum(clone_size_count).astype(float))

    return clone_size, clone_size_frequency


#Load data
TRp_1=pd.read_csv('S1/TRpairs_productive.txt', sep='\t')
TRp_2=pd.read_csv('S2/TRpairs_productive.txt', sep='\t')
TRp_3=pd.read_csv('S3/TRpairs_productive.txt', sep='\t')
TRp_4=pd.read_csv('S4/TRpairs_productive.txt', sep='\t')
TRp_5=pd.read_csv('S5/TRpairs_productive.txt', sep='\t')

#Find exponents for Xmin = 1, 2, 3, 4, 5
CCp_1,CCa_1,CCb_1=powerlaw_pvs_tofit(TRp_1)
CCp_2,CCa_2,CCb_2=powerlaw_pvs_tofit(TRp_2)
CCp_3,CCa_3,CCb_3=powerlaw_pvs_tofit(TRp_3)
CCp_4,CCa_4,CCb_4=powerlaw_pvs_tofit(TRp_4)
CCp_5,CCa_5,CCb_5=powerlaw_pvs_tofit(TRp_5)

sample_list_paired=[CCp_1,CCp_2,CCp_3,CCp_4,CCp_5]
sample_list_alpha=[CCa_1,CCa_2,CCa_3,CCa_4,CCa_5]
sample_list_beta=[CCb_1,CCb_2,CCb_3,CCb_4,CCb_5]

#Create zeros array where each row is a different subject, and each column is a different Xmin.

exponent_p=np.zeros((5,5))
exponent_a=np.zeros((5,5))
exponent_b=np.zeros((5,5))

#Find exponent values.

#Paired
for s in range(5):
    for i in range(5):
        xmin=i+1 # change this if you want. Best results at xmin=1 for our simulations.
        xmax=np.infty # change this if you want.
        sample=sample_list_paired[s]

        data=sample[(sample>=xmin) & (sample<=xmax)] # filters the sample
        N=int(len(data)) # restricted sample size
        C=np.sum(np.log(data))/N # the constant to be used in powerlawpdf function above
        initial_guess=2.0 # method requires an initial guess for the parameter alpha
        alpha_ML=newton(F,initial_guess)
        exponent_p[s,i]=alpha_ML

#Alpha chain
for s in range(5):
    for i in range(5):
        xmin=i+1 # change this if you want. Best results at xmin=1 for our simulations.
        xmax=np.infty # change this if you want.
        sample=sample_list_alpha[s]

        data=sample[(sample>=xmin) & (sample<=xmax)] # filters the sample
        N=int(len(data)) # restricted sample size
        C=np.sum(np.log(data))/N # the constant to be used in powerlawpdf function above
        initial_guess=2.0 # method requires an initial guess for the parameter alpha
        alpha_ML=newton(F,initial_guess)
        exponent_a[s,i]=alpha_ML

#Beta chain
for s in range(5):
    for i in range(5):
        xmin=i+1 # change this if you want. Best results at xmin=1 for our simulations.
        xmax=np.infty # change this if you want.
        sample=sample_list_beta[s]

        data=sample[(sample>=xmin) & (sample<=xmax)] # filters the sample
        N=int(len(data)) # restricted sample size
        C=np.sum(np.log(data))/N # the constant to be used in powerlawpdf function above
        initial_guess=2.0 # method requires an initial guess for the parameter alpha
        alpha_ML=newton(F,initial_guess)
        exponent_b[s,i]=alpha_ML


#Get mean exponent values and standard errors.
mean1_p=np.mean(exponent_p[:,0])
mean2_p=np.mean(exponent_p[:,1])
mean3_p=np.mean(exponent_p[:,2])
mean4_p=np.mean(exponent_p[:,3])
mean5_p=np.mean(exponent_p[:,4])

mean1_a=np.mean(exponent_a[:,0])
mean2_a=np.mean(exponent_a[:,1])
mean3_a=np.mean(exponent_a[:,2])
mean4_a=np.mean(exponent_a[:,3])
mean5_a=np.mean(exponent_a[:,4])

mean1_b=np.mean(exponent_b[:,0])
mean2_b=np.mean(exponent_b[:,1])
mean3_b=np.mean(exponent_b[:,2])
mean4_b=np.mean(exponent_b[:,3])
mean5_b=np.mean(exponent_b[:,4])

std1_p=np.std(exponent_p[:,0])
std2_p=np.std(exponent_p[:,1])
std3_p=np.std(exponent_p[:,2])
std4_p=np.std(exponent_p[:,3])
std5_p=np.std(exponent_p[:,4])

std1_a=np.std(exponent_a[:,0])
std2_a=np.std(exponent_a[:,1])
std3_a=np.std(exponent_a[:,2])
std4_a=np.std(exponent_a[:,3])
std5_a=np.std(exponent_a[:,4])

std1_b=np.std(exponent_b[:,0])
std2_b=np.std(exponent_b[:,1])
std3_b=np.std(exponent_b[:,2])
std4_b=np.std(exponent_b[:,3])
std5_b=np.std(exponent_b[:,4])

sem1_p=std1_p/(5**(1/2))
sem2_p=std2_p/(5**(1/2))
sem3_p=std3_p/(5**(1/2))
sem4_p=std4_p/(5**(1/2))
sem5_p=std5_p/(5**(1/2))

sem1_a=std1_a/(5**(1/2))
sem2_a=std2_a/(5**(1/2))
sem3_a=std3_a/(5**(1/2))
sem4_a=std4_a/(5**(1/2))
sem5_a=std5_a/(5**(1/2))

sem1_b=std1_b/(5**(1/2))
sem2_b=std2_b/(5**(1/2))
sem3_b=std3_b/(5**(1/2))
sem4_b=std4_b/(5**(1/2))
sem5_b=std5_b/(5**(1/2))

#Count clone sizes and frequences for plotting
CS_prod_1, CS_freq_prod_1=powerlaw_paired_toplot(TRp_1)
CS_prod_2, CS_freq_prod_2=powerlaw_paired_toplot(TRp_2)
CS_prod_3, CS_freq_prod_3=powerlaw_paired_toplot(TRp_3)
CS_prod_4, CS_freq_prod_4=powerlaw_paired_toplot(TRp_4)
CS_prod_5, CS_freq_prod_5=powerlaw_paired_toplot(TRp_5)

CSA_prod_1, CSA_freq_prod_1=powerlaw_pvs_toplot(TRp_1, 'NSeqCDR3A')
CSA_prod_2, CSA_freq_prod_2=powerlaw_pvs_toplot(TRp_2, 'NSeqCDR3A')
CSA_prod_3, CSA_freq_prod_3=powerlaw_pvs_toplot(TRp_3, 'NSeqCDR3A')
CSA_prod_4, CSA_freq_prod_4=powerlaw_pvs_toplot(TRp_4, 'NSeqCDR3A')
CSA_prod_5, CSA_freq_prod_5=powerlaw_pvs_toplot(TRp_5, 'NSeqCDR3A')

CSB_prod_1, CSB_freq_prod_1=powerlaw_pvs_toplot(TRp_1, 'NSeqCDR3B')
CSB_prod_2, CSB_freq_prod_2=powerlaw_pvs_toplot(TRp_2, 'NSeqCDR3B')
CSB_prod_3, CSB_freq_prod_3=powerlaw_pvs_toplot(TRp_3, 'NSeqCDR3B')
CSB_prod_4, CSB_freq_prod_4=powerlaw_pvs_toplot(TRp_4, 'NSeqCDR3B')
CSB_prod_5, CSB_freq_prod_5=powerlaw_pvs_toplot(TRp_5, 'NSeqCDR3B')

#Plot mean exponents for each Xmin to compare paired and single.
plt.plot(range(1,6),([mean1_p,mean2_p,mean3_p,mean4_p,mean5_p]),color='red',  marker='o',label='paired')
plt.errorbar(range(1,6), ([mean1_p,mean2_p,mean3_p,mean4_p,mean5_p]), yerr=([sem1_p,sem2_p,sem3_p,sem4_p,sem5_p]), ls='none', color='red', elinewidth=1, capsize=4)
plt.plot(range(1,6),([mean1_a,mean2_a,mean3_a,mean4_a,mean5_a]),color='blue',  marker='o',label='alpha')
plt.errorbar(range(1,6), ([mean1_a,mean2_a,mean3_a,mean4_a,mean5_a]), yerr=([sem1_a,sem2_a,sem3_a,sem4_a,sem5_a]), ls='none', color='blue', elinewidth=1, capsize=4)
plt.plot(range(1,6),([mean1_b,mean2_b,mean3_b,mean4_b,mean5_b]),color='green',  marker='o',label='beta')
plt.errorbar(range(1,6), ([mean1_b,mean2_b,mean3_b,mean4_b,mean5_b]), yerr=([sem1_b,sem2_b,sem3_b,sem4_b,sem5_b]), ls='none', color='green', elinewidth=1, capsize=4)
plt.xlabel('Xmin',fontsize=15)
plt.ylabel('exponent',fontsize=15)
plt.ylim(2,5)
plt.xlim(0.5,5.5)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('mean_powerlaw_pvs_exponents.png')

#Plot exponents for each individual subject for each Xmin to compare CD4 and CD8.
for i in range(5):
    plt.figure()
    plt.plot(np.array(range(1,6)),exponent_p[i,:],'red')
    plt.scatter(np.array(range(1,6)),exponent_p[i,:],marker='o',color='red', edgecolor='black',s=100, label='paired')
    plt.plot(np.array(range(1,6)),exponent_a[i,:],'blue')
    plt.scatter(np.array(range(1,6)),exponent_a[i,:],marker='o',color='blue', edgecolor='black',s=100, label='alpha')
    plt.plot(np.array(range(1,6)),exponent_b[i,:],'green')
    plt.scatter(np.array(range(1,6)),exponent_b[i,:],marker='o',color='green', edgecolor='black',s=100,label='beta')

    plt.xlabel('Xmin',fontsize=15)
    plt.ylabel('exponent',fontsize=15)
    plt.ylim(1,6)
    plt.xlim(0.5,5.5)
    plt.legend(loc='upper left',fontsize=12)
    plt.tight_layout()

    plt.savefig(('S'+str(i+1)+'_powerlaw_pvs_exponents.png'))

CS_prod_list=[CS_prod_1,CS_prod_2,CS_prod_3,CS_prod_4,CS_prod_5]
CS_freq_prod_list=[CS_freq_prod_1,CS_freq_prod_2,CS_freq_prod_3,CS_freq_prod_4,CS_freq_prod_5]
CSA_prod_list=[CSA_prod_1,CSA_prod_2,CSA_prod_3,CSA_prod_4,CSA_prod_5]
CSA_freq_prod_list=[CSA_freq_prod_1,CSA_freq_prod_2,CSA_freq_prod_3,CSA_freq_prod_4,CSA_freq_prod_5]
CSB_prod_list=[CSB_prod_1,CSB_prod_2,CSB_prod_3,CSB_prod_4,CSB_prod_5]
CSB_freq_prod_list=[CSB_freq_prod_1,CSB_freq_prod_2,CSB_freq_prod_3,CSB_freq_prod_4,CSB_freq_prod_5]

#Clone size distribution
xmin=1
for i in range(5):
    plt.figure()
    #paired
    xpmaxP=((Z(exponent_p[i,0],xmin,xmax)/len(TRp_1))**(-1/exponent_p[i,0]))*2
    xpP=np.arange(1,xpmaxP)
    plt.plot(CS_prod_list[i], CS_freq_prod_list[i], 'o', color='red', markeredgecolor='none', markersize=5.0,alpha=1,label='paired')
    plt.plot(xpP,powerlawpdf(xpP,exponent_p[i,0]),'red',linewidth=1)
    #alpha
    xpmaxA=((Z(exponent_a[i,0],xmin,xmax)/len(TRp_1))**(-1/exponent_a[i,0]))*2
    xpA=np.arange(1,xpmaxA)
    plt.plot(CSA_prod_list[i], CSA_freq_prod_list[i], 'o', color='blue', markeredgecolor='none', markersize=5.0,alpha=1,label='alpha')
    plt.plot(xpA,powerlawpdf(xpA,exponent_a[i,0]),'blue',linewidth=1)
    #beta
    xpmaxB=((Z(exponent_b[i,0],xmin,xmax)/len(TRp_1))**(-1/exponent_b[i,0]))*2
    xpB=np.arange(1,xpmaxB)
    plt.plot(CSB_prod_list[i], CSB_freq_prod_list[i], 'o', color='green', markeredgecolor='none', markersize=5.0,alpha=1,label='paired')
    plt.plot(xpB,powerlawpdf(xpB,exponent_b[i,0]),'green',linewidth=1)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('clone size',fontsize=15)
    plt.ylabel('frequency',fontsize=15)
    plt.xlim([10**0-0.2,10**3])
    plt.ylim([10**-6,10**0+1])
    plt.legend()
    plt.tight_layout()
    plt.savefig(('S'+str(i+1)+'_powerlaw_pvs.png'))
