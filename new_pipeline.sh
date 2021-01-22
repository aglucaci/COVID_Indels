#!/bin/bash

## Step 0 
## Get sequences from gisaid

mkdir -p for_observable

## Step 1
## Clean and create geneome length and locale csv
sed 's/ /_/g' gisaid_hcov-19_2021_01_22_16.fasta > underscore.fa
python3 WG_Variation_DayByDay_ColoredByCountry.py underscore.fa /home/aglucaci/SouthAfrica_Indels/for_observable


# Step 2
makeblastdb -in underscore.fa -dbtype nucl -out SA_SARS2_db
blastn -query S.fasta -db SA_SARS2_db -max_hsps 1 -max_target_seqs 4000 -outfmt '6 sseqid sseq' | awk 'BEGIN{FS="\t"; OFS="\n"}{gsub(/-/, "", $2); print ">"$1,$2}' > SA_SARS2_Spike.fasta


# Step 2a, filter the blast results for 90% and above of expected sequence length
# This is important to remove very bad results from blast.
python3 filter_nuc.py SA_SARS2_Spike.fasta 0.9 Run > SA_SARS2_Spike_filtered.fas

# ###
# Step 3

BASEDIR="/home/aglucaci/SouthAfrica_Indels"
HYPHY=$BASEDIR"/hyphy-develop/HYPHYMP"
RES=$BASEDIR"/hyphy-develop/res"
PREMSA=$BASEDIR"/hyphy-analyses/codon-msa/pre-msa.bf"
POSTMSA=$BASEDIR"/hyphy-analyses/codon-msa/post-msa.bf"

FASTA=$BASEDIR"/SA_SARS2_Spike_filtered.fas"
REF=$BASEDIR"/S.fasta"

$HYPHY LIBPATH=$RES $PREMSA --input $FASTA --reference $REF --keep-reference No

mafft --auto $FASTA"_protein.fas" > $FASTA"_protein.msa"

$HYPHY LIBPATH=$RES $POSTMSA --protein-msa $FASTA"_protein.msa" --nucleotide-sequences $FASTA"_nuc.fas" --output $FASTA"_aligned.fasta" --duplicates $FASTA"_duplicates.json"


# Step 3 continued, now that we have a genes msa, process it
# This creates the csv, with genomic positions.
mkdir -p for_observable
python3 gaps_as_indels_withGenomicPosition.py $FASTA"_protein.msa" /home/aglucaci/SouthAfrica_Indels/for_observable S


# This creates the nucleotide chart 
python3 nucl_distribution.py $FASTA"_nuc.fas" > $FASTA"_nuc.fas.table"
mv *.table for_observable



