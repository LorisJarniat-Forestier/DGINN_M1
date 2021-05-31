# -*- coding: utf-8 -*-
"""
Created on Tue May 17 10:26:49 2021

@author: Loris

"""

import pandas as pd
import csv
import matplotlib.pyplot as plt
import numpy as np
import PyQt5
from ete3 import PhyloTree, Tree
from svgutils.compose import *

#============================================================================================

from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='''Tool to plot positive sites from two result table, a tree and a MSA.''')

parser.add_argument('-mod1', metavar="table.log", type=FileType('r'),
    help='the result table for M1', required=True)

parser.add_argument('-mod2', metavar="table.log", type=FileType('r'),
    help='the result table for M2', required=True)

parser.add_argument('-msa', metavar="msa.fa", type=FileType('r'),
    help='the msa file (fasta format)', required=True)

parser.add_argument('-tree', metavar="tree.nhx", type=FileType('r'),
    help='the tree file (NHX format with the "Condition" tag)', required=True)

parser.add_argument('-t', dest="threshold", type=str,
    metavar="\"meth1:0.6\"", help="threshold to filter site by method in the table file",
    default=None)

#parser.add_argument('-out', metavar="output.svg", type=str,
#                    help='the output file (svg format)', required=True)


# names of files to read from
#M1 = PSMD2_CCDS_results_M1.log'
#M2 = PSMD2_CCDS_results_M2.log'
#g_tree=PSMD2_CCDS_M0.nwk.annotated'
#ali_nf=PSMD2_sequences_longestORFs_mafft_mincov_prank.best.fas'
#parser

args = parser.parse_args()

M1 = args.mod1
M2 = args.mod2
ali_nf = args.msa
g_tree = args.tree
out_file = args.out
lim = args.t

# read the data
df_M1 = pd.read_csv(M1, sep='\t')
df_M2 = pd.read_csv(M2, sep='\t')

#site
sites =df_M1.iloc[:,0]

#valeur ll
meanM1 = df_M1.iloc[:,6]
meanM2 = df_M2.iloc[:,6]
diff = meanM2-meanM1

#Mean
mean=[]
for i,element in enumerate(diff):
        if element<-lim :
            mean.append(df_M2.at[i,"mean"])
        else:
            mean.append(0)
    
#============================================================================================
#Création du barplot 

x = sites
y = mean


f = plt.figure()
#f.set_figwidth(50)

plt.bar(x,y)
plt.xlabel("Sites")
plt.ylabel("Omega mean for M2")
plt.title("Omega mean for M2 per sites")

image_format ='svg'
image_name = 'PSMD2_CCDS_results'

f.savefig('PSMD2_CCDS_results_pos.svg')

#=============================================================================================
#Alignement de l'arbre et des séquences

t = PhyloTree(g_tree, ali_nf)
t.render('PSMD2_CCDS_results_tree_ali.svg')


#=============================================================================================
#merging the svg files

#Figure("16cm", "6.5cm",
#       Panel(
#          Text("A", 25, 20),
#          SVG("PSMD2_CCDS_results.svg")
#          ),
#       Panel(
#          Text("B", 25, 20),
#          SVG("C:/Users/Loris/OneDrive/Documents/Stage/DGINN/Example/PSMD2_CCDS_M0.svg").scale(0.5)
#          ).tile(1, 2)
#       )





