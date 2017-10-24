# TCR repertoire - Manuscript
Scripts used to filter and quality control raw JUNO sequencing data of single T cell receptor alpha and beta chains.
Scripts for generating every data figure in the manuscript.

Prerequisites to run all codes: Python 2.7, numpy 1.11.3, matplotlib 2.0.0, pandas 0.10.1, scipy 0.18.1, seaborn 0.7.1, biopython 1.68.

## Filtering and quality control

### FILTERING STEP 1. 
M_Filter_prealign.py - takes in the raw fasta files containing TCR alpha and beta sequences and puts sequence information (droplet barcode, molecular barcode, read count, and cell type) into separate columns followed by the full sequence column.

The created dataframe is separated into two - TRA and TRB for alpha and beta chains.

The reads that are only sequenced once are removed.

The new dataframes are then converted into fasta files for mixcr alignment.

### ALIGNMENT STEP.
M_Alignment_mixcr.sh - takes fasta files created in FILTERING STEP 1 and performs sequence alignment. 

mixcr_2.1.1 software is used: Bolotin, Dmitriy A., Stanislav Poslavsky, Igor Mitrophanov, Mikhail Shugay, Ilgar Z. Mamedov, Ekaterina V. Putintseva, and Dmitriy M. Chudakov. "MiXCR: software for comprehensive adaptive immunity profiling." Nature methods 12, no. 5 (2015): 380-381.

The output is text files containing all alignment information - best V(D)J hits, CDR and FR nucleotide and amino acid sequences and lenghs.

### FILTERING STEP 2.
M_Filter_postalign.py - takes in the text files generated in the ALIGNMENT STEP and removes all sequences that has missing information.

Amino acid sequences that have a stop codon identified by '_' in any region except for FR4 are considered unproductive and are removed.

The new dataframe is then divided into two: 1) is all droplet barcodes that are only associated with one unique sequence, and 2) all the droplet barcodes that are associated with more than one sequence (doublets, triplets, etc.)

Multiple droplet dataframes are then converted into fasta files for clonality analysis.

### CLONALITY STEP.
M_Clonality_mixcr.sh - takes fasta files created in FILTERING STEP 2 and performs clonality analysis.

Finds sequences that are differrent by only one nucleotide and consider them identical by assigning the same clone ID.

The output is a text file containing droplet barcode and clone ID associated with it.

### FILTERING STEP 3.
M_Filter_postclonality.py - takes text files generated in CLONALITY STEP containing clone IDs of each sequence.

Finds sequences with identical clone ID and collapses them into one, keeping the one with a highest read count.

If after collapsing, the droplet is now associated with only one unique sequences, it gets added to the main unique sequence dataframe generated in FILTERING STEP 2.

The new unique sequence dataframe is converted into fasta file to the second clonality analysis.

### CLONALITY STEP 2.
M_Clonality_mixcr2.sh - takes fasta files created in FILTERING STEP 3 and performs clonality analysis.

Finds sequences that are different by only one nucleotide and consider them identical by assigning the same clone ID.

The output is a text file containing droplet barcode and clone ID associated with it.

### ADDING CLONE ID.
M_Add_cloneID.py - takes text files containing clone IDs generated in CLONALITY STEP 2.

Finds sequences that got dropped by mixcr in clonality analysis and removes them.

Adds the corresponding clone IDs as a first column to the main unique sequence datafrane.

### GET ALPHA BETA PAIRS.
M_Pair_recovery.py - takes single chain files and finds matching droplet barcodes for alpha and beta chains.

Generates a dataframe of paired alpha beta chains put in one row. The column names end with A or B for alpha and beta chains in the resulting dataframe.

Some droplets will not have a pair and will be dropped from the paired chain dataframe.

### STRATIFY BY CD4/CD8
M_CD4_CD8_type.py - takes all tab files containing CD4/CD8 information for each droplet barcode, and takes both paired and single chain dataframes.

CD4 and CD8 type will be assigned by a number, and the type will be considered true if that number is >0 for only one type.

Some barcodes have no type assigned with number 0 for both, and some have both types assigned to them - these were dropped.

The output are separate CD4 and CD8 dataframes for each input dataframe.

## Data Analysis

### Figure 2A and B. Supplementary Figure 1 and 2.
M_Powerlaws_pairedvsingle_fits.py - inputs paired alpha/beta dataframes.

The dataframes are split by alpha and beta chain and clone size distributions of paired, single alpha and single beta datasets are compared by fitting power law exponents.

Exponent values are acquired by solving a maximum likelihood equation using Newton-Raphson method.

Figure 2A - power law distributions with a linear fit given the exponent when Xmin=1 for all subject pooled together.

Figure 2B - mean exponent values across all 5 subjects for Xmin=[1,2,3,4,5] for paired and single alpha and beta repertoires.

Supplementary Figure 1 - clone size distributions with a linear fit given the exponent when Xmin=1 for each individual subject.

Supplemetary Figure 2 - exponent values for Xmin=[1,2,3,4,5] of each individual subject for paired and single alpha and beta repertoires.

### Figure 2C. Supplementary Figure 3.
M_unique_chain_count.py - inputs paired alpha/beta dataframes.

Calculates the number of unique alpha sequences each beta sequence pairs with and other way around.

Figure 2C - plot of unique pairing counts for each alpha and beta chain against the frequencies of those counts for all subjects pooled.

Supplementary Figure 3 - unique pairing counts for each alpha and beta chain against the frequencies of those counts of each individual subject.

### Figure 3. Supplementary Figure 4 and 5.
M_heatmaps.py - inputs paired alpha/beta dataframes.

Calculates gene usage and CDR3 length usage across the alpha and beta chains and outputs results in a form of heatmaps.

Frequency = p(state1,state2)/(sample size)

Figure 3 - heatmaps of each gene pair (aVbV, aVbJ, aJbV, aJbJ, aVbD, aJbD) and CDR3 length usage for all subjects pooled.

Supplementary Figure 4 - heatmaps of each gene pair usage of each individual subject.

Supplementary Figure 5 - heatmaps of CDR3 length usage of each individual subject.

### Table 1.
M_mutual_information.py - inputs paired alpha/beta and single alpha and beta dataframes.

Performs mutual information calculations and applies correction: MI(estimated)=MI(true)+a/n+b/n^2+...

Corrections is performed by subsamplling from the dataframes and estimating mutual information for each sample fraction.

Table 1 - MI(true) values for each gene pair within and across the chains as well as CDR3 lengths across the alpha and beta chains.

Additional plots (not present in the manuscript) - plots if 1/(sample size) against the MI estimate. The curve is then fitted to some order polynomial and plotted along with the scatter plot. The order of polynomial can be adjusted given what the plots looks like.

### Supplementary Figure 6.
M_bargraphs.py - inputs single alpha and beta chain dataframes split by CD4 and CD8 type.

Calculates frequences of gene usage (V and J) for each dataframe.

Supplementary Figure 6 - bargraph plots for each gene for each dataframe comparing gene usage across all 5 subjects.

### Figure 4A, B and C. Supplementary Figure 8 and 9.
M_Powerlaws_alldata_fits.py - inputs paired alpha/beta dataframes divided by CD4 and CD8.

The clone size distributions of CD4 and CD8 repertoires are compared by fitting power law exponents. Exponent values are acquired by solving a maximum likelihood equation using Newton-Raphson method.

Figure 4A,B - power law distributions with a linear fit given the exponent when Xmin=1 for CD4 and CD8 for each individual subject.

Figure 4C - mean exponent values across all 5 subjects for Xmin=[1,2,3,4,5] for CD4 and CD8.

Supplementary Figure 8 - power law distribution of CD4 and CD8 repertoires for each individual subject.

Supplementary Figure 9 - exponent values for Xmin=[1,2,3,4,5] for each individual subject for CD4 and CD8 repertoires.

### Supplementary Figure 7.
M_powerlaw_singlechain.py - inputs single alpha and beta chain dataframes divided by CD4 and CD8.

Supplementary Figure 7 - plots of clone size distributions of CD4 and CD8 single chain alpha and beta repertoires of each subject.

### Figure 4D.
M_Clonality.py - inputs paired alpha/beta dataframes divided by CD4 and CD8.

The dataframes are split by alpha and beta chain and clonality scores of paired, single alpha and single beta datasets are compared.

Clonality score is calculated by subtracting Shannon's Entropy from 1.

Figure 4D - clonality scores comparing paired, and single alpha and beta datasets for CD4 and CD8 for each individual subject.

### Figure 5A, B and C.
M_shared_sequences.py - inputs paired alpha/beta and single alpha and beta dataframes.

Figure 5A, B and C - number of shaired paired, and single alpha and beta CDR3 sequences among any 2, 3, 4, and all 5 individuals.

### Figure 5D and E.
M_Shared_sequences_simulation.py - simulation during which an artificial dataset with a clone size distribution following a power law is generated.

The generated dataset is then subsampled twice given the specific fractions. The idea is to see how many shared sequences we can observe when sampling from the same distribution twice.

Figure 5D - the clone size distribution of the generated dataset (n=1 million)

Figure 5E - plot of a fraction of shared sequences between the two subsampled datasets against the fraction to which the dataset was subsampled to.





