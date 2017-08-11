#!/bin/bash

"""
This bash script takes files created during Filter_prealign: TRAtoalign.fasta and TRBtoalign.fasta
and performes alingmnet with mixcr-2.1.1
The output of this step is text files with post alignment information:
V(D)J hits, CDR and FR nucleotide and amino acid sequences and lengths.

BEFORE RUNNING THE CODE:
1.Go to the locations where the fasta files are saved.
2.Input a full address for mixcr location for each mixcr process.
e.g. /Users/your_domain_name/Desktop/Mixcr-2.1.1/mixcr
"""

clear

#go to the file location from Filter_prealign step:
#---Insert command here---

#Make sure to enter full mixcr location for each mixcr commant in this script.

#Align TRA
/mixcr-2.1.1/mixcr align --save-description -c TRA \
TRAtoalign.fasta TRAalignments.vdjca

#Export TRA
/mixcr-2.1.1/mixcr exportAlignments -descrR1 -sequence \
-vHit -dHit -jHit -nFeature CDR1 -lengthof CDR1 -nFeature CDR2 -lengthof CDR2 \
-nFeature CDR3 -lengthof CDR3 -nFeature FR1 -lengthof FR1 -nFeature FR2 -lengthof FR2 \
-nFeature FR3 -lengthof FR3 -nFeature FR4 -lengthof FR4 -aaFeature CDR1 -lengthof CDR1 \
-aaFeature CDR2 -lengthof CDR2 -aaFeature CDR3 -lengthof CDR3 -aaFeature FR1 -lengthof FR1 \
-aaFeature FR2 -lengthof FR2 -aaFeature FR3 -lengthof FR3 -aaFeature FR4 \
-lengthof FR4 TRAalignments.vdjca TRAaligned.txt

echo "Exported TRA file"

#Align TRB files
/mixcr-2.1.1/mixcr align --save-description -c TRB \
TRBtoalign.fasta TRBalignments.vdjca

#Export TRB files
/mixcr-2.1.1/mixcr exportAlignments -descrR1 -sequence \
-vHit -dHit -jHit -nFeature CDR1 -lengthof CDR1 -nFeature CDR2 -lengthof CDR2 \
-nFeature CDR3 -lengthof CDR3 -nFeature FR1 -lengthof FR1 -nFeature FR2 -lengthof FR2 \
-nFeature FR3 -lengthof FR3 -nFeature FR4 -lengthof FR4 -aaFeature CDR1 -lengthof CDR1 \
-aaFeature CDR2 -lengthof CDR2 -aaFeature CDR3 -lengthof CDR3 -aaFeature FR1 -lengthof FR1 \
-aaFeature FR2 -lengthof FR2 -aaFeature FR3 -lengthof FR3 -aaFeature FR4 \
-lengthof FR4 TRBalignments.vdjca TRBaligned.txt

echo "Exported TRB file"

#Remove vdjca files
rm TRAalignments.vdjca
rm TRBalignments.vdjca

echo "COMPLETE"
