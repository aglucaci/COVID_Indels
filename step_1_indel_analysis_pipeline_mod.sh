# bash

# ###################################################################### #
# Information
# ###################################################################### #
# @Author: Alexander G Lucaci
# More to SARS-CoV-2 evolution than recombination and purifying selection
# Deletion as an adaptive evolutionary mechanism.
# ###################################################################### #

# ###################################################################### #
# Goals
# ###################################################################### #
# Part 1. Patterns, what can we detect
# Part 2. Processes, what is causing viral deletions?
# Part 3. Prediction, what can we learn from previous deletions to inform us, city by city/
# Part 4. Product, whats a tool or resource for this? Can we create a profile of deletions to place incoming/new sequences into.
# ###################################################################### #


# ###################################################################### #
# Pre-pipeline
# ###################################################################### #
# Pre-pipeline-steps.
# Filter for complete genomes (> 29k bp)
# Filter for Human hosts
# Download from GISAID all of the complete genomes (>29,000bp)
# Download the metadata .xls format
# Convert metadata xls format to csv.
# ###################################################################### #

# ###################################################################### #
# Pipeline supporting software
# ###################################################################### #
# Requirements:
# Python >3.6
# Mafft
# FastTree
# blast+
# parallel
# PRANK
# HyPhy
# ###################################################################### #

# ###################################################################### #
# Pipeline supporting software to install.
# ###################################################################### #
# conda install -c bioconda prank
# conda install -c bioconda hyphy
# ###################################################################### #


# ###################################################################### #
# @Usage: parallel ::: "bash step_1_indel_analysis_pipeline.sh"
# ###################################################################### #


#Need a filtering step
#Human, genome size >27000, steven doesnt do this.


# ###################################################################### #
# Declaress
# ###################################################################### #
ANALYSISDATE="05112020-genomicpos"

#Lets setup the latest HyPhy
#bash install_hyphy.sh

HYPHY="hyphy-develop/HYPHYMP"
RES="hyphy-develop/res"
MAFFT="mafft"

#GISAID="../data/GISAID/"$ANALYSISDATE"/gisaid_cov2020_sequences.fasta"
GISAID="../Extra_data/2020-05-11_gisaid/2020-05-11.fasta"

#confirm it exists, or exit
if [ -f "$GISAID" ]
then
    echo "# $GISAID found."
else
    echo "# $GISAID not found."
    exit
fi

#Correct for spaces, turn them into underscores
echo "# Running sed to correct for spaces"
sed -i '' 's/ /_/g' $GISAID

REFCDSDIR="../data/reference_genes"
REFGENOME=""

# ###################################################################### #
# Pipeline, only use numbered steps (X)
# ###################################################################### #
OUTPUTDIR="../analysis/"$ANALYSISDATE
mkdir $OUTPUTDIR

#echo "Generating violin plot"
# (1) Generate a violin plot of whole genome
#python WG_Variation.py $GISAID $OUTPUTDIR

echo "Generating Country by Country variation plot"
# (2) Colored locale by locale
python WG_Variation_DayByDay_ColoredByCountry.py $GISAID $OUTPUTDIR


# ###################################################################### #
# End of file.
# ###################################################################### #


