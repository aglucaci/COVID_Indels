#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:57:14 2020

@author: alexander g. lucaci

https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/

https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
Reference genome size 29,903
"""

# Python Imports
import os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import plotly.graph_objects as go
import plotly
import numpy as np

fasta_file = sys.argv[1]
output_dir = sys.argv[2]

sequence_lengths = {}

with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"): 
        ID = record.id
        SEQ = record.seq
        if len(SEQ) < 27000: continue
        num_gaps = str(SEQ).count("-") / 3 #codon aware alignment
        seq_length = len(str(SEQ))

        #print(ID, seq_length)
        sequence_lengths[ID] = seq_length
    #end for
#end with


#init figure
fig = go.Figure()

array = []
ID_array = []

for item in sequence_lengths:
    array += [sequence_lengths[item]]
    ID_array += [item]

fig.add_trace(go.Violin(y=array, name="WG Size Variation (SARS-CoV-2)", hovertext=ID_array, points="all"))

#overlay the box plots
fig.update_layout(barmode='overlay', title="Whole Genome Size Variation (SARS-CoV-2) ", xaxis_title="WG SARS-CoV-2", yaxis_title="WG Size (bp)",)

#plot
#fig.update_traces(orientation='h') # horizontal box plots#
#fig.show()
plotly.offline.plot(fig, filename=output_dir+'/WG_indel_analysis.html')
