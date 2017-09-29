# TCR repertoire - Manuscript
Scripts used to filter and quality control raw JUNO sequencing data of single T cell receptor alpha and beta chains.
Scripts for generating every data figure in the manuscript post

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




