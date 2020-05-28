#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 14 11:03:23 2020

@author: alexander g lucaci
"""

# Imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import numpy as np
from itertools import islice
import sys
import plotly.express as px


#input_fasta = "../analysis/04142020/hyphy_alignments/nucl/gisaid_cov2020_sequences.fasta.S_nuc.fas"
input_fasta = sys.argv[1]

print("#,ID,Nucleotides,Date")

with open(input_fasta, "r") as handle:
    for n, record in enumerate(SeqIO.parse(handle, "fasta")): 
        #if n == 1: break
        ID = record.id
        SEQ = record.seq
        #date = ID.split("_")[-3] + "-" + ID.split("_")[-2] + "-" + ID.split("_")[-1]
        #date = ID.split("_")[-3] + "-" + ID.split("_")[-2] + "-" + ID.split("_")[-1]
        date = ID.split("_")[-2]
        date = date[:4] + "-" + date[4:6] + "-" + date[6:8]
        if date[:4] == "2022": continue
        #print(date)
        
        #Debug.
        #print(ID, len(SEQ), date)
        
        print(",".join([str(n), str(ID), str(len(SEQ)), date]))
    #end for
#end with


# END OF FILE
