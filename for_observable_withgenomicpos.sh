#bash

#updated for each date.
analysisdate="05112020"

#make the directory
mkdir -p ../analysis/$analysisdate/observable

#exit

python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.S.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.M.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.N.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.ORF1a.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.ORF3a.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.ORF6.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.ORF7a.compressed.fas ../analysis/$analysisdate/observable/
python gaps_as_indels_withGenomicPosition.py ../Extra_data/2020-05-11/sequences.ORF8.compressed.fas ../analysis/$analysisdate/observable/

