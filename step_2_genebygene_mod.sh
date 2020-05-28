# bash

# ###################################################################### #
# Declares
# ###################################################################### #
ANALYSISDATE="05112020"


# ###################################################################### #
# (1) This step will create a figure "gap_analysis.html"
# ###################################################################### #
#echo "Step 1"
#@Usage: python gene by gene variationanalysis.py {Input dir} {Output dir}

#input is codon aware alignments
#Codons
#python variation_analysis.py "../2020-05-11" "../analysis/"$ANALYSISDATE

#exit

#Proteins
#python variation_analysis.py "../analysis/"$ANALYSISDATE"/hyphy_alignments/protein_all_withref" "../analysis/"$ANALYSISDATE


# ###################################################################### #
# This step will
# output_csv = fasta_data.split("/")[-1] + "_table.csv"
# ###################################################################### #
#(2)
echo "Step 2"

# for file in protein msa directry
FILES=../Extra_data/2020-05-11/*.msa
for msa in $FILES; do
    #continue
    echo "Processing: "$msa
    # python gaps_as_indels.py {input fasta} {output directory}
    # outputs
    python gaps_as_indels.py $msa "../analysis/"$ANALYSISDATE
done

exit

# ###################################################################### #
# This step will

# Cause I dont have unaligned proteins,
# Lets look at the nucleotide distibution for S.
# ###################################################################### #
#(3)
echo "Step 3"

INPUT="../analysis/04142020/hyphy_alignments/nucl/gisaid_cov2020_sequences.fasta.S_nuc.fas"
python nucl_distribution.py $INPUT > $INPUT".csv"


# ###################################################################### #
# End of file.
# ###################################################################### #
