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
ANALYSISDATE="05112020"

#Lets setup the latest HyPhy
#bash install_hyphy.sh

HYPHY="hyphy-develop/HYPHYMP"
RES="/Users/user/Documents/hyphy-develop/res"
MAFFT="/usr/local/bin/mafft"

#GISAID="../data/GISAID/"$ANALYSISDATE"/gisaid_cov2020_sequences.fasta"
GISAID="../Extra_data/2020-05-11_gisaid/2020-05-11.fasta"

#confirm it exists, or exit
if [ -f "$GISAID" ]
then
    echo "$GISAID found."
else
    echo "$GISAID not found."
    exit
fi

#Correct for spaces, turn them into underscores
echo "Running sed to correct for spaces"
sed -i '' 's/ /_/g' $GISAID

REFCDSDIR="../data/reference_genes"
REFGENOME=""

# ###################################################################### #
# Pipeline, only use numbered steps (X)
# ###################################################################### #
OUTPUTDIR="../analysis/"$ANALYSISDATE
mkdir $OUTPUTDIR

HYPHYALIGNMENTS=$OUTPUTDIR"/hyphy_alignments"
mkdir $HYPHYALIGNMENTS

#hyphy standalone analyses
echo "Checking for hyphy-analyses"
git clone https://github.com/veg/hyphy-analyses.git

PREMSA="hyphy-analyses/codon-msa/pre-msa.bf"
POSTMSA="hyphy-analyses/codon-msa/post-msa.bf"

echo "Generating violin plot"
# (1) Generate a violin plot of whole genome
python WG_Variation.py $GISAID $OUTPUTDIR

echo "Generating Country by Country variation plot"
# (2) Colored locale by locale
python WG_Variation_DayByDay_ColoredByCountry.py $GISAID $OUTPUTDIR




echo Done


exit



HYPHYALIGNMENTS=$OUTPUTDIR"/hyphy_alignments"

echo "Creating output directory for alignments: "$HYPHYALIGNMENTS
mkdir $HYPHYALIGNMENTS

function run_a_gene {
    GENE=$1
    REFERENCE_SEQUENCE=$2
    TRIM_FROM=$3
    TRIM_TO=$4
    N_FRAC=$5
    REFERENCE_PROTEIN=$6
    
    #If check.
    
    #will output to data directory, lets rename and move it to alignments directory for this analysis date.
    # uses the Gene reference sequence
    $HYPHY LIBPATH=$RES $PREMSA --input $GISAID --reference $REFERENCE_SEQUENCE --trim-from $TRIM_FROM --trim-to $TRIM_TO --N-fraction $N_FRAC --E 0.01
    f="$(basename -- $GISAID)"
    
    mv ${GISAID}_protein.fas $HYPHYALIGNMENTS/${f}.${GENE}_protein.fas
    mv ${GISAID}_nuc.fas $HYPHYALIGNMENTS/${f}.${GENE}_nuc.fas
    
    #Protein - msa
    $MAFFT $HYPHYALIGNMENTS/${f}.${GENE}_protein.fas > $HYPHYALIGNMENTS/${f}.${GENE}.msa
    
    #Run twice, compressed and uncompressed.
    $HYPHY LIBPATH=$RES $POSTMSA --protein-msa $HYPHYALIGNMENTS/${f}.${GENE}.msa --nucleotide-sequences $HYPHYALIGNMENTS/${f}.${GENE}_nuc.fas --output $HYPHYALIGNMENTS/${f}.${GENE}.compressed.fas --duplicates $HYPHYALIGNMENTS/${f}.${GENE}.duplicates.json
    
    $HYPHY LIBPATH=$RES $POSTMSA --protein-msa $HYPHYALIGNMENTS/${f}.${GENE}.msa --nucleotide-sequences $HYPHYALIGNMENTS/${f}.${GENE}_nuc.fas --output $HYPHYALIGNMENTS/${f}.${GENE}.all.fas --compress No
    
    # Add the reference back to the alignment so we can make calls.
    # This is for codon aware.
    $MAFFT --add $REFERENCE_SEQUENCE --reorder $HYPHYALIGNMENTS/${f}.${GENE}.all.fas > $HYPHYALIGNMENTS/${f}.${GENE}.codon_all_withref.fas
    
    # This is for the protein msa.
    # NEED THE PROTEIN Refernece seq.
    $MAFFT --add $REFERENCE_PROTEIN --reorder $HYPHYALIGNMENTS/${f}.${GENE}.msa > $HYPHYALIGNMENTS/${f}.${GENE}.protein_all_withref.fas
}
#end method

# Main subroutine. ----------------------------------------------------------------------------------------------------------------------------------------
echo "Running Codon aware alignment"

run_a_gene "S" "../data/reference_genes/S.fasta" "20000" "27000" 0.005 "../data/reference_proteins/S.fasta"
run_a_gene "M" "../data/reference_genes/M.fasta" "25000" "30000" 0.01 "../data/reference_proteins/M.fasta"
run_a_gene "N" "../data/reference_genes/N.fasta" "26000" "35000" 0.01 "../data/reference_proteins/N.fasta"
run_a_gene "ORF3a" "../data/reference_genes/ORF3a.fasta" "24000" "27000" 0.05 "../data/reference_proteins/ORF3a.fasta"
run_a_gene "ORF6" "../data/reference_genes/ORF6.fasta" "26000" "30000" 0.05 "../data/reference_proteins/ORF6.fasta"
run_a_gene "ORF7a" "../data/reference_genes/ORF7a.fasta" "26000" "35000" 0.05 "../data/reference_proteins/ORF7a.fasta"
run_a_gene "ORF8" "../data/reference_genes/ORF8.fasta" "26000" "35000" 0.05 "../data/reference_proteins/ORF8.fasta"
run_a_gene "ORF1a" "../data/reference_genes/ORF1a.fasta" "1" "15000" 0.001 "../data/reference_proteins/ORF1a.fasta"
run_a_gene "ORF1b" "../data/reference_genes/ORF1b.fasta" "12000" "24000" 0.001 "../data/reference_proteins/ORF1b.fasta"
run_a_gene "E" "../data/reference_genes/E.fasta" "26000" "27000" 0.01 "../data/reference_proteins/E.fasta"
run_a_gene "ORF7b" "../data/reference_genes/ORF7b.fasta" "27000" "28000" 0.01 "../data/reference_proteins/ORF7b.fasta"
run_a_gene "ORF10" "../data/reference_genes/ORF10.fasta" "29000" "29800" 0.01 "../data/reference_proteins/ORF10.fasta"
run_a_gene "ORF1ab" "../data/reference_genes/ORF1ab.fasta" "1" "24000" 0.001 "../data/reference_proteins/ORF1ab.fasta"

#Move/Organize the outputs ---------------------------------------------------------------------------------------------------------------------------------
mkdir $HYPHYALIGNMENTS/nucl
mkdir $HYPHYALIGNMENTS/protein
mkdir $HYPHYALIGNMENTS/all
mkdir $HYPHYALIGNMENTS/compressed
mkdir $HYPHYALIGNMENTS/duplicates
mkdir $HYPHYALIGNMENTS/codon_all_withref
mkdir $HYPHYALIGNMENTS/protein_all_withref

mv $HYPHYALIGNMENTS/*_nuc.fas $HYPHYALIGNMENTS/nucl
mv $HYPHYALIGNMENTS/*_protein.fas $HYPHYALIGNMENTS/protein
mv $HYPHYALIGNMENTS/*.all.fas $HYPHYALIGNMENTS/all
mv $HYPHYALIGNMENTS/*.compressed.fas $HYPHYALIGNMENTS/compressed
mv $HYPHYALIGNMENTS/*.duplicates.json $HYPHYALIGNMENTS/duplicates
mv $HYPHYALIGNMENTS/*.msa $HYPHYALIGNMENTS/protein
mv $HYPHYALIGNMENTS/*.codon_all_withref.fas $HYPHYALIGNMENTS/codon_all_withref
mv $HYPHYALIGNMENTS/*.protein_all_withref.fas $HYPHYALIGNMENTS/protein_all_withref

# ###################################################################### #
# End of file.
# ###################################################################### #


