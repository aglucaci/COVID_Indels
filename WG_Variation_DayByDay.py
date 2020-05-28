#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:57:14 2020

@author: alexander g. lucaci

https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/

https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
Reference genome size 29,903
"""

# Imports
import os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import plotly.graph_objects as go
import plotly
import numpy as np


# Main subroutine
sequence_lengths = {}
count = 0

fasta_file = sys.argv[1]
output_dir = sys.argv[2]

with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"): 
        ID = record.id
        SEQ = record.seq
        num_gaps = str(SEQ).count("-") / 3 #codon aware alignment
        seq_length = len(str(SEQ))

        #print(ID, seq_length)
        #sequence_lengths[ID] = seq_length
        date = str(ID.split("|")[-1])
        if date[:4] == "2020" or date[:4] == "2019" and len(date) > 4:
            sequence_lengths[count] = {}
            sequence_lengths[count]["ID"] = ID
            sequence_lengths[count]["Date"] = date
            sequence_lengths[count]["WG_Length"] = seq_length
            count += 1
    #end for
#end with

#print(sequence_lengths)


import pandas as pd

df = pd.DataFrame.from_dict(sequence_lengths, orient='index')

#df.rename( columns={ 0 :'ID'}, inplace=True )

#print(df.iloc[:, 0].tolist())
#df.index.name = 'new_name'

print(df)

#dfObj = pd.DataFrame(sequence_lengths)

#print(dfObj)

#sys.exit(1)

#init figure
fig = go.Figure()

fig = go.Figure(data=go.Scatter(
                x=df["Date"],
                y=df["WG_Length"],
                mode='markers',
                hovertext=df["ID"],
                marker=dict(
                     #color='rgb(255, 178, 102)',
                     color=df["WG_Length"],
                     size=9,
                     line=dict(
                        color='DarkSlateGrey',
                        width=1
                      )
               )
))

fig.update_layout(
    title='Price vs. Average Playtime',
    xaxis_title='Price (GBP)',
    yaxis_title='Average Playtime (Minutes)'
)

fig.show()

plotly.offline.plot(fig, filename=output_dir+'/WG_indel_analysis_daybyday.html')

print("done")

















#for item in sequence_lengths:
#    #print(item)
#    for item2 in sequence_lengths[item2]
#        print("\t", "Date:", sequence_lengths[item][0], "\n\t", "WG length (nt):",sequence_lengths[item][1])




























sys.exit(1)


#lets transform the above dict.
dated_dict = {}

for ID in sequence_lengths:
    #e.g. hCoV-19/USA/WA-UW56/2020|EPI_ISL_415621|2020-03-09
    date = str(ID.split("|")[-1])
    if date[:4] == "2020" or date[:4] == "2019":
        #print(date)
        dated_dict[date] = {}
        
        try:
            dated_dict[date][ID] += sequence_lengths[ID]
        except: 
            dated_dict[date][ID] = sequence_lengths[ID]

        #end try
    #end if
#end for

print(dated_dict)
#sys.exit(1)

# PLOT

#init figure
fig = go.Figure()


for the_date in dated_dict:
    # '2020-03-09': {'hCoV-19/USA/WA-UW42/2020|EPI_ISL_415607|2020-03-09': 29857}, '2020-03-08': {'hCoV-19/USA/WA-UW43/2020|EPI_ISL_415608|2020-03-08': 29882}, '2020-03-06': {'hCoV-19/Slovakia/SK-BMC1/2020|EPI_ISL_417877|2020-03-06': 29903}, '2020-03-15': {'hCoV-19/Iceland/8/2020|EPI_ISL_417864|2020-03-15': 29903}, '2020-03-16': {'hCoV-19/Australia/VIC12/2020|EPI_ISL_416518|2020-03-16': 29872}
    data_array = []
    ID_array = []
    dates = []

    for ID in dated_dict[the_date]:
        data_array += [dated_dict[the_date][ID]]
        ID_array += [ID]
    #end inner for

    #print(ID)
    #print(data_array)
    fig.add_trace(go.Box(y=data_array, name=the_date, hovertext=ID_array))
#end for
    
    
#overlay the box plots
fig.update_layout(barmode='overlay', title="Whole Genome Size Variation (SARS-CoV-2) [03-29-2020]", xaxis_title="WG SARS-CoV-2", yaxis_title="WG Size (bp)",)


fig.show()

#plotly.offline.plot(fig, filename='../WG_indel_analysis_daybyday.html')




sys.exit(1)


data_array = []
ID_array = []
dates = []

for the_date in dated_dict:
    data_array = []
    data_array.append(dated_dict[the_date][1])
    ID_array.append(dated_dict[the_date][0])
    dates += [the_date]
    
    print("\t", data_array)
    
    fig.add_trace(go.Box(y=data_array, x=[the_date], name="WG Size Variation (SARS-CoV-2) Day by Day", hovertext=ID_array))
    
#overlay the box plots
fig.update_layout(barmode='overlay', title="Whole Genome Size Variation (SARS-CoV-2) [03-29-2020]", xaxis_title="WG SARS-CoV-2", yaxis_title="WG Size (bp)",)


fig.show()
plotly.offline.plot(fig, filename='../WG_indel_analysis_daybyday.html')


sys.exit(1)







#init figure
fig = go.Figure()

array = []
ID_array = []

for item in sequence_lengths:
    array += [sequence_lengths[item]]
    ID_array += [item]

fig.add_trace(go.Violin(y=array, name="WG Size Variation (SARS-CoV-2)", hovertext=ID_array, points="all"))

#overlay the box plots
fig.update_layout(barmode='overlay', title="Whole Genome Size Variation (SARS-CoV-2) [03-29-2020]", xaxis_title="WG SARS-CoV-2", yaxis_title="WG Size (bp)",)

#plot
#fig.update_traces(orientation='h') # horizontal box plots#
fig.show()
plotly.offline.plot(fig, filename='../WG_indel_analysis.html')


# =============================================================================
# End of file.
# =============================================================================
