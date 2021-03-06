#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 10:28:15 2020

@author: alexander lucaci

conda install -c conda-forge biopython
conda install -c plotly plotly
conda install -c synthicity prettytable
"""

# =============================================================================
# Imports
# =============================================================================

import os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import AlignIO
from Bio.Align import AlignInfo
import plotly.graph_objects as go
import plotly
import numpy as np
import pandas as pd
#from prettytable import PrettyTable

# =============================================================================
# Declares
# =============================================================================

# Gene file
fasta_data = sys.argv[1]
reference_gene = ""

print("# (Model) Input file is:", fasta_data)
print("# Creating consensus sequence.")

alignment = AlignIO.read(fasta_data, 'fasta')
summary_align = AlignInfo.SummaryInfo(alignment)
Consensus_Sequence = summary_align.dumb_consensus(float(0.5))

print("# Done.")

#Site: [Ydel, IDs]
data_dict = {}

#data_dict["ID"] = []
# Keep count of total number of gap characters
count = 1

#output_directory = "../analysis/04052020"
output_directory = sys.argv[2]

# =============================================================================
# Main
# =============================================================================
print("# Looping over fasta")
counter = 1
with open(fasta_data, "r") as handle:
    for i, record in enumerate(SeqIO.parse(handle, "fasta")): 
        #continue
        ID = record.id
        SEQ = record.seq
        #num_gaps = str(SEQ).count("-") / 3 #codon aware alignment
        seq_length = len(str(SEQ))
        #print(ID, seq_length)

        #Loop over sequence.
        for n, char in enumerate(str(SEQ)):
            if char == "-":
                data_dict[counter] = [ID, int(n+1), Consensus_Sequence[n]]
                counter += 1
                count += 1    
            #end if
        #end for
    #end for
#end with

#print(data_dict)
print("# Done..")
print("# Character count '-':", count, "\n")

#print(data_dict)
#sys.exit(10)
print("# Creating dataframe")
df = pd.DataFrame.from_dict(data_dict, orient='index', columns=['Sequence_ID', 'Site_Indel', 'Consensus(AA)'])
#df = pd.DataFrame.from_dict(data_dict)
print("# Done...")

print("# Sorting dataframe by Site number")
df = df.sort_values(by=['Site_Indel'])
print("# Done. ")

output_csv = fasta_data.split("/")[-1] + "_table.csv"

print("# Saving dataframe to file:", output_csv)
df.to_csv(os.path.join(output_directory, output_csv), index=False)
print("# Done....")

print(df)
    
# =============================================================================
# End of file
# =============================================================================

