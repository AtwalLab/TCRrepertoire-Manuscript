#!/bin/bash

"""
This bash script takes files created during Filter_postclonality:
TRAsFinal_productive.fasta and TRBsFinal_productive.fasta and performes clonality analysis with mixcr-2.1.1
The output of this step is text files with droplet barcode and clone ID.
Identical clone ID means, identical clone even if the sequence is different.

BEFORE RUNNING THE CODE:
1.Go to the locations where the fasta files are saved.
2.Input a full address for mixcr location for each mixcr process.
e.g. /Users/your_domain_name/Desktop/Mixcr-2.1.1/mixcr
"""

clear

#go to the file location from Filter_prealign step:
#---Insert command here---

#Make sure to enter full mixcr location for each mixcr commant in this script.

#Align TRA productive
/sonas-hs/atwal/hpc/home/kgrigait/new_data/Scripts/mixcr align --save-description -c TRA TRAsFinal_productive.fasta alignments.vdjca
#Assemble
/sonas-hs/atwal/hpc/home/kgrigait/new_data/Scripts/mixcr assemble --index indexFile alignments.vdjca clones.clns
#Export
/sonas-hs/atwal/hpc/home/kgrigait/new_data/Scripts/mixcr exportAlignments -cloneId indexFile -descrR1 alignments.vdjca TRAs_final_prod_clones.txt
#remove intermediate files
rm alignments.vdjca
rm clones.clns
rm indexFile

#Align TRB productive
/sonas-hs/atwal/hpc/home/kgrigait/new_data/Scripts/mixcr align --save-description -c TRB TRBsFinal_productive.fasta alignments.vdjca
#Assemble
/sonas-hs/atwal/hpc/home/kgrigait/new_data/Scripts/mixcr assemble --index indexFile alignments.vdjca clones.clns
#Export
/sonas-hs/atwal/hpc/home/kgrigait/new_data/Scripts/mixcr exportAlignments -cloneId indexFile -descrR1 alignments.vdjca TRBs_final_prod_clones.txt
#remove intermediate files
rm alignments.vdjca
rm clones.clns
rm indexFile
