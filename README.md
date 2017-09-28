# TCR repertoire - Manuscript
Scripts used to filter and quality control raw JUNO sequencing data of single T cell receptor alpha and beta chains.
Scripts for generating every data figure in the manuscript post

## Filtering and quality control

### FILTERING STEP 1. 
M_filter_prealign.py - takes in the raw fasta files containing TCR alpha and beta sequences and puts sequence information (droplet barcode, molecular barcode, read count, and cell type) into separate columns followed by the full sequence column.
The created dataframe is separated into two - TRA and TRB for alpha and beta chains.
The reads that are only sequenced once are removed.
The new dataframes are then converted into fasta files for mixcr alignment.

### ALIGNMENT STEP.
M_Alignment_mixcr.sh - takes fasta files created in FILTERING STEP 1 and performs sequence alignment. 
mixcr_2.1.1 software is used: Bolotin, Dmitriy A., Stanislav Poslavsky, Igor Mitrophanov, Mikhail Shugay, Ilgar Z. Mamedov, Ekaterina V. Putintseva, and Dmitriy M. Chudakov. "MiXCR: software for comprehensive adaptive immunity profiling." Nature methods 12, no. 5 (2015): 380-381.

The output is text files containing all alignment information - best V(D)J hits, CDR and FR nucleotide and amino acid sequences and lenghs.


