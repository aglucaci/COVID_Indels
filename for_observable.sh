#bash

#updated for each date.
analysisdate="05112020"

#make the directory
mkdir -p ../analysis/$analysisdate/observable

echo "Copying variation csv to observable directory"

#this needs to be uploaded to observablehq
cp ../analysis/$analysisdate/wg_variation_daybyday_coloredbycountry.csv "../analysis/$analysisdate/observable/"$analysisdate"_wg_variation_daybyday_coloredbycountry.csv"
#Need to filter by bad dates.

echo "Creating nucleotide distribution tables"
#python nucl_distribution.py ../analysis/$analysisdate/hyphy_alignments/nucl/gisaid_cov2020_sequences.fasta.S_nuc.fas > ../analysis/$analysisdate/observable/gisaid_cov2020_sequences.fasta.S_nuc.fas.table
#../Extra_data/2020-05-11/*.msa
for nuclseq in ../Extra_data/2020-05-11/*_nuc.fas; do
    f="$(basename -- $nuclseq)"
    echo "Processing: "$nuclseq
    echo "    Output: "../analysis/$analysisdate/observable/$f".table"
    python nucl_distribution.py $nuclseq > ../analysis/$analysisdate/observable/$f".table"
    #python nucl_distribution.py $nuclseq
done

#for nuclseq in ../analysis/$analysisdate/hyphy_alignments/nucl/*.fas; do
#    f="$(basename -- $nuclseq)"
#    python nucl_distribution.py $nuclseq > ../analysis/$analysisdate/observable/$f".table"
#done


echo "Creating output for Spike"

#python gaps_as_indels.py ../analysis/$analysisdate/hyphy_alignments/protein_all_withref/gisaid_cov2020_sequences.fasta.S.protein_all_withref.fas ../analysis/$analysisdate/observable/

#python gaps_as_indels.py {Protein alignment} {output_dir}
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.S.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.M.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.N.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.ORF1a.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.ORF3a.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.ORF6.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.ORF7a.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels.py ../Extra_data/2020-05-11/sequences.ORF8.compressed.fas ../analysis/$analysisdate/observable/

echo .
