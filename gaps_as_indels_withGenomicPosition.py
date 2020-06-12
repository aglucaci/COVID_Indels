#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 10:28:15 2020

@author: Alexander G. Lucaci

conda install -c conda-forge biopython
conda install -c plotly plotly
conda install -c synthicity prettytable

SARS2 Reference Genome
https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
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

# =============================================================================
# Declares
# =============================================================================
# Gene file
fasta_data = sys.argv[1]
reference_gene = ""

print("# Input file is:", fasta_data)
print("# Creating consensus sequence.")

alignment = AlignIO.read(fasta_data, 'fasta')
summary_align = AlignInfo.SummaryInfo(alignment)
Consensus_Sequence = summary_align.dumb_consensus(float(0.5))

print("# Done creating consensus sequence.")
data_dict = {}
count = 1
output_directory = sys.argv[2]

Gene_map_to_Genomic = {}

gene = fasta_data.split("/")[-1].replace(".msa", "").replace("sequences.", "")
print("# Gene we are analyzing:", gene)

Gene_map_to_Genomic["S"] = 21563
Gene_map_to_Genomic["M"] = 26523
Gene_map_to_Genomic["N"] = 28274
Gene_map_to_Genomic["ORF1a"] = 266
Gene_map_to_Genomic["ORF3a"] = 25393
Gene_map_to_Genomic["ORF6"] = 27202
Gene_map_to_Genomic["ORF7a"] = 27394
Gene_map_to_Genomic["ORF8"] = 27894

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
            if char == "-": #Will look at the nucleotide position.
                site_position = int(n+1)
                genomic_position = Gene_map_to_Genomic[gene] + (site_position*3) 
                
                data_dict[counter] = [ID, str(genomic_position) + "-" + str(genomic_position + 2), site_position, Consensus_Sequence[n]]
                counter += 1
                count += 1    
            #end if
        #end for
    #end for
#end with

                    
print("# Done..")
print("# Character count '-':", count, "\n")


print("# Creating dataframe")
df = pd.DataFrame.from_dict(data_dict, orient='index', columns=['Sequence_ID', 'Genomic_Position(NT)', 'Site_Indel(AA)', 'Consensus(AA)'])
print("# Done...")

print("# Sorting dataframe by Site number")
df = df.sort_values(by=['Site_Indel(AA)'])
print("# Done. ")


output_csv = "genomic_" + fasta_data.split("/")[-1] + "_table.csv"

print("# Saving dataframe to file:", output_csv)
df.to_csv(os.path.join(output_directory, output_csv), index=False)
print("# Done....")

print(df)

    
# =============================================================================
# End of file
# =============================================================================

