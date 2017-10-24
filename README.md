# TCRrepertoire-Manuscript
Jupyter notebooks for generating every data figure in the manuscript using Juno sequencing data of single T cell receptor alpha and beta chains. 

For quality control and filtering see Filtering scripts.

Prerequisites to run all codes:
Python 2.7
numpy 1.11.3
matplotlib 2.0.0
pandas 0.10.1
scipy 0.18.1
seaborn 0.7.1


## Data Analysis

### Figure 2A and B. Supplementary Figure 1 and 2.
M_Powerlaws_pairedvsingle_fits - inputs paired alpha/beta dataframes.

The dataframes are split by alpha and beta chain and clone size distributions of paired, single alpha and single beta datasets are compared by fitting power law exponents.

Exponent values are acquired by solving a maximum likelihood equation using Newton-Raphson method.

Figure 2A - power law distributions with a linear fit given the exponent when Xmin=1 for all subject pooled together.

Figure 2B - mean exponent values across all 5 subjects for Xmin=[1,2,3,4,5] for paired and single alpha and beta repertoires.

Supplementary Figure 1 - clone size distributions with a linear fit given the exponent when Xmin=1 for each individual subject.

Supplemetary Figure 2 - exponent values for Xmin=[1,2,3,4,5] of each individual subject for paired and single alpha and beta repertoires.

### Figure 2C. Supplementary Figure 3.
M_unique_chain_count - inputs paired alpha/beta dataframes.

Calculates the number of unique alpha sequences each beta sequence pairs with and other way around.

Figure 2C - plot of unique pairing counts for each alpha and beta chain against the frequencies of those counts for all subjects pooled.

Supplementary Figure 3 - unique pairing counts for each alpha and beta chain against the frequencies of those counts of each individual subject.

### Figure 3. Supplementary Figure 4 and 5.
M_heatmaps - inputs paired alpha/beta dataframes.

Calculates gene usage and CDR3 length usage across the alpha and beta chains and outputs results in a form of heatmaps.

Frequency = p(state1,state2)/(sample size)

Figure 3 - heatmaps of each gene pair (aVbV, aVbJ, aJbV, aJbJ, aVbD, aJbD) and CDR3 length usage for all subjects pooled.

Supplementary Figure 4 - heatmaps of each gene pair usage of each individual subject.

Supplementary Figure 5 - heatmaps of CDR3 length usage of each individual subject.

### Table 1.
M_mutual_information - inputs paired alpha/beta and single alpha and beta dataframes.

Performs mutual information calculations and applies correction: MI(estimated)=MI(true)+a/n+b/n^2+...

Corrections is performed by subsamplling from the dataframes and estimating mutual information for each sample fraction.

Table 1 - MI(true) values for each gene pair within and across the chains as well as CDR3 lengths across the alpha and beta chains.

Additional plots (not present in the manuscript) - plots if 1/(sample size) against the MI estimate. The curve is then fitted to some order polynomial and plotted along with the scatter plot. The order of polynomial can be adjusted given what the plots looks like.

### Supplementary Figure 6.
M_bargraphs - inputs single alpha and beta chain dataframes split by CD4 and CD8 type.

Calculates frequences of gene usage (V and J) for each dataframe.

Supplementary Figure 6 - bargraph plots for each gene for each dataframe comparing gene usage across all 5 subjects.

### Figure 4A, B and C. Supplementary Figure 8 and 9.
M_Powerlaws_alldata_fits - inputs paired alpha/beta dataframes divided by CD4 and CD8.

The clone size distributions of CD4 and CD8 repertoires are compared by fitting power law exponents. Exponent values are acquired by solving a maximum likelihood equation using Newton-Raphson method.

Figure 4A,B - power law distributions with a linear fit given the exponent when Xmin=1 for CD4 and CD8 for each individual subject.

Figure 4C - mean exponent values across all 5 subjects for Xmin=[1,2,3,4,5] for CD4 and CD8.

Supplementary Figure 8 - power law distribution of CD4 and CD8 repertoires for each individual subject.

Supplementary Figure 9 - exponent values for Xmin=[1,2,3,4,5] for each individual subject for CD4 and CD8 repertoires.

### Supplementary Figure 7.
M_powerlaw_singlechain - inputs single alpha and beta chain dataframes divided by CD4 and CD8.

Supplementary Figure 7 - plots of clone size distributions of CD4 and CD8 single chain alpha and beta repertoires of each subject.

### Figure 4D.
M_Clonality - inputs paired alpha/beta dataframes divided by CD4 and CD8.

The dataframes are split by alpha and beta chain and clonality scores of paired, single alpha and single beta datasets are compared.

Clonality score is calculated by subtracting Shannon's Entropy from 1.

Figure 4D - clonality scores comparing paired, and single alpha and beta datasets for CD4 and CD8 for each individual subject.

### Figure 5A, B and C.
M_shared_sequences - inputs paired alpha/beta and single alpha and beta dataframes.

Figure 5A, B and C - number of shaired paired, and single alpha and beta CDR3 sequences among any 2, 3, 4, and all 5 individuals.

### Figure 5D and E.
M_Shared_sequences_simulation - simulation during which an artificial dataset with a clone size distribution following a power law is generated.

The generated dataset is then subsampled twice given the specific fractions. The idea is to see how many shared sequences we can observe when sampling from the same distribution twice.

Figure 5D - the clone size distribution of the generated dataset (n=1 million)

Figure 5E - plot of a fraction of shared sequences between the two subsampled datasets against the fraction to which the dataset was subsampled to.
