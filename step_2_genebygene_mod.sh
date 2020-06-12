# bash

# ###################################################################### #
# Declares
# ###################################################################### #
ANALYSISDATE="05112020-genomicpos"

# ###################################################################### #
# This step will
# output_csv = fasta_data.split("/")[-1] + "_table.csv"
# ###################################################################### #
#(2)
echo "# (Model) In Step 2"

# for file in protein msa directry
FILES=../Extra_data/2020-05-11/*.msa

for msa in $FILES; do
    #continue
    echo "# Processing: "$msa
    # python gaps_as_indels.py {input fasta} {output directory}
    # outputs
    #python gaps_as_indels.py $msa "../analysis/"$ANALYSISDATE
    python gaps_as_indels_withGenomicPosition.py $msa "../analysis/"$ANALYSISDATE
done


# ###################################################################### #
# End of file.
# ###################################################################### #
