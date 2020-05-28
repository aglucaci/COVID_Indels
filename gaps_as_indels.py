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
os.system("clear")

# Gene file
#fasta_data = "data/gisaid_cov2020_sequences_03272020.fasta"
#fasta_data = "../analysis/04052020/hyphy_alignments/withref/gisaid_cov2020_sequences.fasta.S.withref.fas"
#fasta_data = "../analysis/04052020/hyphy_alignments/protein/gisaid_cov2020_sequences.fasta.S.msa"
fasta_data = sys.argv[1]
reference_gene = ""

print("# Input file is:", fasta_data)
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
                #print(i, ID, Consensus_Sequence[n], n+1, "del") # I want to add what the reference char would be at this position.
                #if n + 1 == 265: continue
                #try:
                    #data_dict[ID] += [n+1, Consensus_Sequence[n]]
                    #data_dict["ID"][ID] += [n+1, Consensus_Sequence[n]]
                    #data_dict["ID"] += [ID]
                    #data_dict["ID"] += {ID : {"SiteDeletion" : [n+1]}, {"Consensus(AA)": Consensus_Sequence[n]}}
                    #pass
                    
                #    data_dict["ID"] += [ID]
                #    data_dict["ID"][ID] = "1"
                #except:
                    #data_dict[ID] = [n+1, Consensus_Sequence[n]]
                    #pass
                #    data_dict["ID"] = [ID]
                #    data_dict["ID"][ID] = "0"
                
                #data_dict[n+1] = {}
                #data_dict[counter] = [{"ID": ID}, {"SiteDeletion" : [n+1]}, {"Consensus(AA)": Consensus_Sequence[n]}]
                # columns=['first_col', 'second_col', 'third_col']
                
                data_dict[counter] = [ID, int(n+1), Consensus_Sequence[n]]
                counter += 1
                #end try
                count += 1    
            #end if
        #end for
        
    #end for
#end with



# pretty table.
                    
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


sys.exit(10)
# End here for now.

"""
What do I want to know?
The site that has experienced the most deletions (frequency)
The sequence with the most deletions (len of the dict[id])
"""


#print({k: v for k, v in sorted(data_dict.items(), key=lambda item: item[1])})
"""
most_deletions_dict = {}

for k in data_dict:
    #key is the id
    num = len(data_dict[k]) / 2 #number of deletions
    most_deletions_dict[k] = num
#end for

sorted_most_deletions_dict = {k: v for k, v in sorted(most_deletions_dict.items(), key=lambda item: item[1], reverse=True)}
    
for item in sorted_most_deletions_dict:
    #print(item, sorted_most_deletions_dict[item])
    print(item, data_dict[item])
    
"""    



# =============================================================================
# Want to get a better handle on sites
# =============================================================================
# Site, Del, IDs
sites_del_ids = {}

for key in data_dict:
    ID = key
    number_of_exchanges = len(data_dict[key]) / 2
    
    #sites are the evens
    #amino acid is odds
    sites = data_dict[key][0::2]
    
    amino_acids = data_dict[key][1::2]
    
    for n, item in enumerate(sites):
        #item is the site number
        # the site then contains the AA and the IDs
        try:
            sites_del_ids[item][amino_acids[n]] += [ID]
        except:
            sites_del_ids[item] = {amino_acids[n] : [ID]}
        #end try
    #end for
#end for
            
print(sites_del_ids)
# =============================================================================
# Output table.
# =============================================================================

#sys.exit(100)

from prettytable import PrettyTable
    
pt = PrettyTable()

pt.field_names = ["#", "Site", "Potential Deletion", "IDs where this is present"]
count = 1

for site in sites_del_ids:

    for item in sites_del_ids[site]:
        
        pt.add_row([count, site, item, sites_del_ids[site][item]])

        count += 1
        

print(pt)
    
# =============================================================================
# End of file
# =============================================================================

