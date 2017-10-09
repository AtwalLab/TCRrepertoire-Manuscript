"""
This script inputs paired alpha/beta dataframes divided by CD4 and CD8.
The clone size distribution of CD4 and CD8 repertoires are compared by fitting
powerlaw exponents.
The script plots powerlaw distributions and acquires the exponents 'gamma' by using a
maximum likelihood approach. Exponents are acquired for Xmin=[1,2,3,4,5].
"""

from __future__ import division
from scipy.optimize import newton
from scipy.special import zeta
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#Defining the functions.
def powerlawpdf(x,gamma):
    """The power law probability function
    Input: x - array of clone sizes; gamma - desired exponent.
    Output: array of probabilities of each clone size x."""
    return (x**(-gamma))/Z(gamma,xmin,xmax)

def Z(gamma,xmin=1,xmax=np.infty):
    """The normalization function Z.
    Note that default arguments of xmin and xmax make Z equivalent to Riemann zeta function which is already
    implemented in Scipy as zeta(gamma,1)"""
    if np.isfinite(xmax):
        s=0
        for i in xrange(xmin,xmax+1):
            s+=(1/(i**gamma))
    else:
        s=zeta(gamma,xmin)
    return s

def F(gamma):
    """The optimization function F(gamma). C is the second term in the definition of F(gamma) and is independent of
    gamma. C need to be defined before running the function. Function Z must be run beforehand as well."""
    h = 1e-8
    Z_prime = (Z(gamma+h,xmin,xmax) - Z(gamma-h,xmin,xmax))/(2*h)
    return (Z_prime/Z(gamma,xmin,xmax))+C

def powerlaw_paired_tofit(df):
    """This function calculates clone sizes.
    Input: pandas dataframe of paired sequences.
    Returns: clone sizes of paired alpha beta sequences."""

    df['joint']=df['NSeqCDR3A']+'_'+df['NSeqCDR3B']

    return np.array(pd.DataFrame(df.joint.value_counts())).flatten()

def powerlaw_paired_toplot(df):
    """This function calculates clone sizes and their frequencies.
    Input: pandas dataframe of paired sequences.
    Returns: numpy array of clone sizes, numpy array of clone size frequencies of paired dataframes."""

    df['joint']=df['NSeqCDR3A']+'_'+df['NSeqCDR3B']

    clone_size_count=np.array(pd.DataFrame(df.joint.value_counts()).joint.value_counts())
    clone_size=np.array(pd.DataFrame(df.joint.value_counts()).joint.value_counts().index)
    clone_size_frequency=np.divide(clone_size_count, np.sum(clone_size_count).astype(float))

    return clone_size, clone_size_frequency

#Load data paired alpha/beta CD4 and CD8 sequencing data.
TRpCD4_prod_1=pd.read_csv('S1/TRcd4_productive.txt', sep='\t')
TRpCD4_prod_2=pd.read_csv('S2/TRcd4_productive.txt', sep='\t')
TRpCD4_prod_3=pd.read_csv('S3/TRcd4_productive.txt', sep='\t')
TRpCD4_prod_4=pd.read_csv('S4/TRcd4_productive.txt', sep='\t')
TRpCD4_prod_5=pd.read_csv('S5/TRcd4_productive.txt', sep='\t')

TRpCD8_prod_1=pd.read_csv('S1/TRcd8_productive.txt', sep='\t')
TRpCD8_prod_2=pd.read_csv('S2/TRcd8_productive.txt', sep='\t')
TRpCD8_prod_3=pd.read_csv('S3/TRcd8_productive.txt', sep='\t')
TRpCD8_prod_4=pd.read_csv('S4/TRcd8_productive.txt', sep='\t')
TRpCD8_prod_5=pd.read_csv('S5/TRcd8_productive.txt', sep='\t')


#Find exponents for Xmin = 1, 2, 3, 4, 5
#Find clone sizes of paired alpha beta sequences for CD4 and CD8 types of each subject.
CC4_1=powerlaw_paired_tofit(TRpCD4_prod_1)
CC4_2=powerlaw_paired_tofit(TRpCD4_prod_2)
CC4_3=powerlaw_paired_tofit(TRpCD4_prod_3)
CC4_4=powerlaw_paired_tofit(TRpCD4_prod_4)
CC4_5=powerlaw_paired_tofit(TRpCD4_prod_5)

CC8_1=powerlaw_paired_tofit(TRpCD8_prod_1)
CC8_2=powerlaw_paired_tofit(TRpCD8_prod_2)
CC8_3=powerlaw_paired_tofit(TRpCD8_prod_3)
CC8_4=powerlaw_paired_tofit(TRpCD8_prod_4)
CC8_5=powerlaw_paired_tofit(TRpCD8_prod_5)

#Put each array into a list in preparation for 'for-loops' when acquiring exponents.
sample_list_4=[CC4_1,CC4_2,CC4_3,CC4_4,CC4_5]
sample_list_8=[CC8_1,CC8_2,CC8_3,CC8_4,CC8_5]

#Create zeros array where each row is a different subject, and each column is a different Xmin.
exponent4=np.zeros((5,5))
exponent8=np.zeros((5,5))

#Find exponent values.
#CD4
for s in range(5):
    for i in range(5):
        xmin=i+1 #define Xmin
        xmax=np.infty #define Xmax (in this case always infinity)
        sample=sample_list_4[s] #Choose a sample from previously prepared sample lists.
        data=sample[(sample>=xmin) & (sample<=xmax)] # filters the sample by Xmin and Xmax
        N=int(len(data)) # Find sample size post filtering
        C=np.sum(np.log(data))/N # Find the constant to be used in powerlawpdf function above
        initial_guess=2.0 #initial guess for the parameter gamma
        gamma_ML=newton(F,initial_guess)
        exponent4[s,i]=gamma_ML

#CD8
for s in range(5):
    for i in range(5):
        xmin=i+1
        xmax=np.infty
        sample=sample_list_8[s]
        data=sample[(sample>=xmin) & (sample<=xmax)]
        N=int(len(data))
        C=np.sum(np.log(data))/N
        initial_guess=2.0
        gamma_ML=newton(F,initial_guess)
        exponent8[s,i]=gamma_ML

#Get mean exponent values and standard errors for each Xmin across the 5 subjects.
mean1_4=np.mean(exponent4[:,0])
mean2_4=np.mean(exponent4[:,1])
mean3_4=np.mean(exponent4[:,2])
mean4_4=np.mean(exponent4[:,3])
mean5_4=np.mean(exponent4[:,4])

mean1_8=np.mean(exponent8[:,0])
mean2_8=np.mean(exponent8[:,1])
mean3_8=np.mean(exponent8[:,2])
mean4_8=np.mean(exponent8[:,3])
mean5_8=np.mean(exponent8[:,4])

std1_4=np.std(exponent4[:,0])
std2_4=np.std(exponent4[:,1])
std3_4=np.std(exponent4[:,2])
std4_4=np.std(exponent4[:,3])
std5_4=np.std(exponent4[:,4])

std1_8=np.std(exponent8[:,0])
std2_8=np.std(exponent8[:,1])
std3_8=np.std(exponent8[:,2])
std4_8=np.std(exponent8[:,3])
std5_8=np.std(exponent8[:,4])

sem1_4=std1_4/(5**(1/2))
sem2_4=std2_4/(5**(1/2))
sem3_4=std3_4/(5**(1/2))
sem4_4=std4_4/(5**(1/2))
sem5_4=std5_4/(5**(1/2))

sem1_8=std1_8/(5**(1/2))
sem2_8=std2_8/(5**(1/2))
sem3_8=std3_8/(5**(1/2))
sem4_8=std4_8/(5**(1/2))
sem5_8=std5_8/(5**(1/2))

#Now get clone sizes and their frequencies for plotting.
CS_cd4_prod_1, CS_cd4_freq_prod_1=powerlaw_paired_toplot(TRpCD4_prod_1)
CS_cd4_prod_2, CS_cd4_freq_prod_2=powerlaw_paired_toplot(TRpCD4_prod_2)
CS_cd4_prod_3, CS_cd4_freq_prod_3=powerlaw_paired_toplot(TRpCD4_prod_3)
CS_cd4_prod_4, CS_cd4_freq_prod_4=powerlaw_paired_toplot(TRpCD4_prod_4)
CS_cd4_prod_5, CS_cd4_freq_prod_5=powerlaw_paired_toplot(TRpCD4_prod_5)

CS_cd8_prod_1, CS_cd8_freq_prod_1=powerlaw_paired_toplot(TRpCD8_prod_1)
CS_cd8_prod_2, CS_cd8_freq_prod_2=powerlaw_paired_toplot(TRpCD8_prod_2)
CS_cd8_prod_3, CS_cd8_freq_prod_3=powerlaw_paired_toplot(TRpCD8_prod_3)
CS_cd8_prod_4, CS_cd8_freq_prod_4=powerlaw_paired_toplot(TRpCD8_prod_4)
CS_cd8_prod_5, CS_cd8_freq_prod_5=powerlaw_paired_toplot(TRpCD8_prod_5)

#Plot FIGURE 4A
#Clone size distribution
#CD4
xmin=1
#Subject 1
xpmax1=((Z(exponent4[0,0],xmin,xmax)/len(CC4_1))**(-1/exponent4[0,0]))*2
xp1=np.arange(1,xpmax1)
plt.plot(CS_cd4_prod_1, CS_cd4_freq_prod_1, 'o', color='red', markeredgecolor='none', markersize=5.0,alpha=1,label='S1')
plt.plot(xp1,powerlawpdf(xp1,exponent4[0,0]),'red',linewidth=1)
#Subject 2
xpmax2=((Z(exponent4[1,0],xmin,xmax)/len(CC4_2))**(-1/exponent4[1,0]))*2
xp2=np.arange(1,xpmax2)
plt.plot(CS_cd4_prod_2, CS_cd4_freq_prod_2, 'o', color='orange', markeredgecolor='none', markersize=5.0,alpha=1,label='S2')
plt.plot(xp2,powerlawpdf(xp2,exponent4[1,0]),'orange',linewidth=1)
#Subject 3
xpmax3=((Z(exponent4[2,0],xmin,xmax)/len(CC4_3))**(-1/exponent4[2,0]))*2
xp3=np.arange(1,xpmax3)
plt.plot(CS_cd4_prod_3, CS_cd4_freq_prod_3, 'o', color='black', markeredgecolor='none', markersize=5.0,alpha=1,label='S3')
plt.plot(xp3,powerlawpdf(xp3,exponent4[2,0]),'black',linewidth=1)
#Subject 4
xpmax4=((Z(exponent4[3,0],xmin,xmax)/len(CC4_4))**(-1/exponent4[3,0]))*2
xp4=np.arange(1,xpmax4)
plt.plot(CS_cd4_prod_4, CS_cd4_freq_prod_4, 'o', color='green', markeredgecolor='none', markersize=5.0,alpha=1,label='S4')
plt.plot(xp4,powerlawpdf(xp4,exponent4[3,0]),'green',linewidth=1)
#Subject 5
xpmax5=((Z(exponent4[4,0],xmin,xmax)/len(CC4_5))**(-1/exponent4[4,0]))*2
xp5=np.arange(1,xpmax5)
plt.plot(CS_cd4_prod_5, CS_cd4_freq_prod_5, 'o', color='blue', markeredgecolor='none', markersize=5.0,alpha=1,label='S5')
plt.plot(xp5,powerlawpdf(xp5,exponent4[4,0]),'blue',linewidth=1)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**3])
plt.ylim([10**-6,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('CD4_fitted_powerlaws_paper.png')


#Plot FIGURE 4B
#CD8
xmin=1
#Subject 1
#calculate max clone size used to generate a straight line fit.
xpmax1=((Z(exponent8[0,0],xmin,xmax)/len(CC8_1))**(-1/exponent8[0,0]))*2
xp1=np.arange(1,xpmax1) #generate a straight line fit.
#Plot a scatter plot of clone size distribution.
plt.plot(CS_cd8_prod_1, CS_cd8_freq_prod_1, 'o', color='red', markeredgecolor='none', markersize=5.0,alpha=1,label='S1')
#Plot the linear fit of the clone size distribution.
plt.plot(xp1,powerlawpdf(xp1,exponent8[0,0]),'red',linewidth=1)
#Subject 2
xpmax2=((Z(exponent8[1,0],xmin,xmax)/len(CC8_2))**(-1/exponent8[1,0]))*2
xp2=np.arange(1,xpmax2)
plt.plot(CS_cd8_prod_2, CS_cd8_freq_prod_2, 'o', color='orange', markeredgecolor='none', markersize=5.0,alpha=1,label='S2')
plt.plot(xp2,powerlawpdf(xp2,exponent8[1,0]),'orange',linewidth=1)
#Subject 3
xpmax3=((Z(exponent8[2,0],xmin,xmax)/len(CC8_3))**(-1/exponent8[2,0]))*2
xp3=np.arange(1,xpmax3)
plt.plot(CS_cd8_prod_3, CS_cd8_freq_prod_3, 'o', color='black', markeredgecolor='none', markersize=5.0,alpha=1,label='S3')
plt.plot(xp3,powerlawpdf(xp3,exponent8[2,0]),'black',linewidth=1)
#Subject 4
xpmax4=((Z(exponent8[3,0],xmin,xmax)/len(CC8_4))**(-1/exponent8[3,0]))*2
xp4=np.arange(1,xpmax4)
plt.plot(CS_cd8_prod_4, CS_cd8_freq_prod_4, 'o', color='green', markeredgecolor='none', markersize=5.0,alpha=1,label='S4')
plt.plot(xp4,powerlawpdf(xp4,exponent8[3,0]),'green',linewidth=1)
#Subject 5
xpmax5=((Z(exponent8[4,0],xmin,xmax)/len(CC8_5))**(-1/exponent4[4,0]))*2
xp5=np.arange(1,xpmax5)
plt.plot(CS_cd8_prod_5, CS_cd8_freq_prod_5, 'o', color='blue', markeredgecolor='none', markersize=5.0,alpha=1,label='S5')
plt.plot(xp5,powerlawpdf(xp5,exponent8[4,0]),'blue',linewidth=1)

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**3])
plt.ylim([10**-6,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('CD8_fitted_powerlaws_paper.png')

#Plot FIGURE 4C
#Mean exponents for each Xmin to compare CD4 and CD8 clone size distributions.
plt.plot(range(1,6),([mean1_4,mean2_4,mean3_4,mean4_4,mean5_4]),color='red',  marker='o',label='CD4')
plt.errorbar(range(1,6), ([mean1_4,mean2_4,mean3_4,mean4_4,mean5_4]), yerr=([sem1_4,sem2_4,sem3_4,sem4_4,sem5_4]), ls='none', color='red', elinewidth=1, capsize=4)
plt.plot(range(1,6),([mean1_8,mean2_8,mean3_8,mean4_8,mean5_8]),color='purple',  marker='o',label='CD8')
plt.errorbar(range(1,6), ([mean1_8,mean2_8,mean3_8,mean4_8,mean5_8]), yerr=([sem1_8,sem2_8,sem3_8,sem4_8,sem5_8]), ls='none', color='purple', elinewidth=1, capsize=4)
plt.xlabel('Xmin',fontsize=15)
plt.ylabel('exponent',fontsize=15)
plt.ylim(1,6)
plt.xlim(0.5,5.5)
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('mean_powerlaw_exponents.png')

#Plot SUPPLEMENTARY FIGURE 9
#Exponents for each individual subject for each Xmin to compare CD4 and CD8 clone size distributions.
for i in range(5):
    plt.figure()
    plt.plot(np.array(range(1,6)),exponent4[i,:],'red')
    plt.scatter(np.array(range(1,6)),exponent4[i,:],marker='o',color='red', edgecolor='black',s=100, label='CD4')
    plt.plot(np.array(range(1,6)),exponent8[i,:],'purple')
    plt.scatter(np.array(range(1,6)),exponent8[i,:],marker='o',color='purple', edgecolor='black',s=100,label='CD8')

    plt.xlabel('Xmin',fontsize=15)
    plt.ylabel('exponent',fontsize=15)
    plt.ylim(1,9)
    plt.xlim(0.5,5.5)
    plt.legend(loc='upper left',fontsize=12)
    plt.tight_layout()

    plt.savefig(('S'+str(i+1)+'_powerlaw_exponents.png'))

#SUPPLEMENTARY FIGURE 8
#CD4 vs CD8
#Subject 1
plt.figure()
plt.plot(CS_cd4_prod_1, CS_cd4_freq_prod_1, 'o', color='red', markeredgecolor='none', markersize=8.0,alpha=1,label='CD4')
plt.plot(CS_cd8_prod_1, CS_cd8_freq_prod_1, 'o', color='purple', markeredgecolor='none', markersize=8.0,alpha=1,label='CD8')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**3])
plt.ylim([10**-6,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('S1_CD4vCD8_powerlaw_paper.png')

#Subject 2
plt.figure()
plt.plot(CS_cd4_prod_2, CS_cd4_freq_prod_2, 'o', color='red', markeredgecolor='none', markersize=8.0,alpha=1,label='CD4')
plt.plot(CS_cd8_prod_2, CS_cd8_freq_prod_2, 'o', color='purple', markeredgecolor='none', markersize=8.0,alpha=1,label='CD8')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**3])
plt.ylim([10**-6,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('S2_CD4vCD8_powerlaw_paper.png')

#Subject 3
plt.figure()
plt.plot(CS_cd4_prod_3, CS_cd4_freq_prod_3, 'o', color='red', markeredgecolor='none', markersize=8.0,alpha=1,label='CD4')
plt.plot(CS_cd8_prod_3, CS_cd8_freq_prod_3, 'o', color='purple', markeredgecolor='none', markersize=8.0,alpha=1,label='CD8')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**3])
plt.ylim([10**-6,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('S3_CD4vCD8_powerlaw_paper.png')

#Subject 4
plt.figure()
plt.plot(CS_cd4_prod_4, CS_cd4_freq_prod_4, 'o', color='red', markeredgecolor='none', markersize=8.0,alpha=1,label='CD4')
plt.plot(CS_cd8_prod_4, CS_cd8_freq_prod_4, 'o', color='purple', markeredgecolor='none', markersize=8.0,alpha=1,label='CD8')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**3])
plt.ylim([10**-6,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('S4_CD4vCD8_powerlaw_paper.png')

#Subject 5
plt.figure()
plt.plot(CS_cd4_prod_5, CS_cd4_freq_prod_5, 'o', color='red', markeredgecolor='none', markersize=8.0,alpha=1,label='CD4')
plt.plot(CS_cd8_prod_5, CS_cd8_freq_prod_5, 'o', color='purple', markeredgecolor='none', markersize=8.0,alpha=1,label='CD8')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('clone size',fontsize=15)
plt.ylabel('frequency',fontsize=15)
plt.xlim([10**0-0.2,10**3])
plt.ylim([10**-6,10**0+1])
plt.legend(fontsize=12)
plt.tight_layout()
plt.savefig('S5_CD4vCD8_powerlaw_paper.png')
