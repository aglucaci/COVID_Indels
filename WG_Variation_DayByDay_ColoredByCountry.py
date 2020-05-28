#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 29 13:57:14 2020

@author: alexander g. lucaci

https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/

https://www.ncbi.nlm.nih.gov/nuccore/NC_045512
Reference genome size 29,903
"""

# =============================================================================
# Imports
# =============================================================================
import os, sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import plotly.graph_objects as go
import plotly
import numpy as np
import pandas as pd
import plotly.express as px

# =============================================================================
# Declares
# =============================================================================
sequence_lengths = {}
count = 0
fasta_file = sys.argv[1]
output_dir = sys.argv[2]

# =============================================================================
# Helper functions
# =============================================================================
def lookup_collection_date():
    pass

# =============================================================================
# Main, lets gather the data into a dict.
# =============================================================================
with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"): 
        ID = record.id
        SEQ = record.seq
        #num_gaps = str(SEQ).count("-") / 3 #codon aware alignment
        seq_length = len(str(SEQ))
        if seq_length < 27000: continue
        
        #print(ID, seq_length)
        #sequence_lengths[ID] = seq_length
        date = str(ID.split("|")[-1])
        
        
        locale = ID.split("/")[1]
        #print(ID, date)
        if " " in ID: print(ID, date)
        
        if date[:4] == "2020" or date[:4] == "2019" and len(date) >= 7:
            sequence_lengths[count] = {}
            sequence_lengths[count]["ID"] = ID
            sequence_lengths[count]["Date"] = date
            sequence_lengths[count]["WG_Length"] = seq_length
            sequence_lengths[count]["Locale"] = locale
            #sequence_lengths[count]["Color"] = "red"
            count += 1
        else:
            print("# ERR:", ID, seq_length)
    #end for
#end with

#print(sequence_lengths)


# =============================================================================
# Lets create a df from dict.
# =============================================================================

df = pd.DataFrame.from_dict(sequence_lengths, orient='index')

#df.rename( columns={ 0 :'ID'}, inplace=True )
#print(df.iloc[:, 0].tolist())
#df.index.name = 'new_name'

# DEBUG
#print(set([sequence_lengths[x]["Locale"] for x in sequence_lengths]), len(set([sequence_lengths[x]["Locale"] for x in sequence_lengths])))


# =============================================================================
# Lets create a discrete color set
# =============================================================================
color_set = {}
color_count = 1

for place in set([sequence_lengths[x]["Locale"] for x in sequence_lengths]):
    color_set[place] = color_count
    color_count += 1*100
#end for

# DEBUG
#print(color_set)

# =============================================================================
# add color set to dataframe
# =============================================================================

# create column
df['Color'] = 0

#col x row
for index, row in df.iterrows(): 
    # print(df["Locale"][i]) 
    #df["Color"][i] = color_set[df["Locale"][i]]
    #df["Color"][j] = 0
    #print(row["Locale"], color_set[row["Locale"]])
    
    df["Color"][index] = color_set[row["Locale"]]
#end for
    
print(df)

#dfObj = pd.DataFrame(sequence_lengths)
#print(dfObj)

# =============================================================================
# Save df to file.
# =============================================================================

df.to_csv(output_dir+'/wg_variation_daybyday_coloredbycountry.csv')

# =============================================================================
# Lets plot.
# =============================================================================


#fig = px.scatter(df, x="Date", y="WG_Length", color="Locale",
#                 size=[1]*len([x for x in sequence_lengths]), hover_data=["ID"], marginal_y="histogram", marginal_x="violin")

fig = px.scatter(df, x="Date", y="WG_Length", color="Locale",
                 size=[1]*len([x for x in sequence_lengths]), hover_data=["ID"], marginal_y="histogram", marginal_x="histogram")
fig.update_layout(
    title='[SARS-CoV-2] Whole genome size variation n=' + str(len([x for x in sequence_lengths])),
    xaxis={'title': 'Collection Date'},
    yaxis={'title': 'Whole Genome size (nt)'})

#fig.show()

plotly.offline.plot(fig, filename=output_dir+'/WG_indel_analysis_daybyday_coloredbylocale.html')







sys.exit(1)


# =============================================================================
# Does not reach this.
# =============================================================================

#init figure
fig = go.Figure()

fig = go.Figure(data=go.Scatter(
                x=df["Date"],
                y=df["WG_Length"],
                mode='markers',
                hovertext=df["ID"],
                marker=dict(
                     #color='rgb(255, 178, 102)',
                     color=df["Locale"],
                     size=9,
                    
               )
))

fig.update_layout(
    title='[SARS-CoV-2] Whole genome size variation n=' + str(len([x for x in sequence_lengths])),
    xaxis_title='Date',
    yaxis_title='Whole Genome size (nt)',
    showlegend=True,
    legend=color_set
)

#fig.show()

plotly.offline.plot(fig, filename='../WG_indel_analysis_daybyday_coloredbylocale.html')





# =============================================================================
# End of file.
# =============================================================================