
#Imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

# Usage@: python3 filter_nuc.py {FASTA} {THRESHOLD} {MODE: Test or Run} > output.fasta
FASTA = sys.argv[1]
THRESHOLD = float(sys.argv[2])
MODE = sys.argv[3]
# Declares

seq_length_list = []

# Main
if MODE == "Test":
    print("# Processing:", FASTA)

with open(FASTA, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"): 
        #ID = record.id
        SEQ = record.seq
        #DESC = record.description

        seq_length_list.append(len(SEQ))
    #end for
#end with

average = sum(seq_length_list)/len(seq_length_list)
seq_threshold = int(max(seq_length_list) * THRESHOLD)
#How many would pass?
passed_threshold_count = len([x for x in seq_length_list if x > seq_threshold])
capture_percentage = passed_threshold_count / len(seq_length_list)

if MODE == "Test": 
    print("# Average sequence length:", int(average))
    print("# Maximum sequence length:", int(max(seq_length_list)))                          
    print("# Minimum sequence length:", int(min(seq_length_list)))
    print("# NT Threshold (bp):", seq_threshold)
    print("# Total number of sequences:", len(seq_length_list))
    print("# Number of sequences which would pass this filter:", passed_threshold_count, "or", float(capture_percentage))
#end if




with open(FASTA, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        ID = record.id
        SEQ = record.seq
        DESC = record.description

        #seq_length_list.append(len(SEQ))

        if len(SEQ) > seq_threshold:
            pass
            if MODE == "Run":
                print(">" + ID + "\n" + SEQ)
            #end if
    #end for
#end with



# End of file


