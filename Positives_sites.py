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

"""from argparse import ArgumentParser, FileType
parser = ArgumentParser(description='''Tool to plot convergent sites from a result table, a tree and a MSA.
ex:\n\n
python script/plot_convergent_sites.py  -tsv example/tree4plot.tsv -msa example/tree4plot.fa -tree example/tree4plot_annotated.nw -out example/tree4plot.svg -meth Meth3,Meth1 -t Meth3:0.8,Meth1:70
''')
parser.add_argument('-csvM1', metavar="table.csv", type=FileType('r'),
    help='the result table for M1', required=True)
parser.add_argument('-csvM2', metavar="table.csv", type=FileType('r'),
    help='the result table for M2', required=True)"""



# names of files to read from
M1 = 'C:/Users/Loris/OneDrive/Documents/Stage/DGINN/Example/PSMD2_CCDS_results_M1.log'
M2 = 'C:/Users/Loris/OneDrive/Documents/Stage/DGINN/Example/PSMD2_CCDS_results_M2.log'
g_tree='C:/Users/Loris/OneDrive/Documents/Stage/DGINN/Example/PSMD2_CCDS_M0.nwk.annotated'
ali_nf='C:/Users/Loris/OneDrive/Documents/Stage/DGINN/Example/PSMD2_sequences_longestORFs_mafft_mincov_prank.best.fas'
#parser



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
        if element<-0.6 :
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

f.savefig('PSMD2_CCDS_results_pos.png')

#=============================================================================================
#Alignement de l'arbre et des séquences

t = PhyloTree(g_tree, ali_nf)
t.render('PSMD2_CCDS_results_tree_ali.png')


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







