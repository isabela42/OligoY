#!/usr/bin/python
# coding=utf-8

#################################################################
#         This script was developed by Isabela Almeida          #
#                    mb.isabela42@gmail.com                     #
#                                                               #
# For ploting YGS2 (Y Genome Scan v02) results according to the #
# stablished vsck and pvscuk thresholds.                        #
#                                                               #
# This was created while I was a MSc Bioinformatics student at  #
# the Universidade de São Paulo (USP/Brazil). For this project  #
# I had the guidence of Eduardo Dupim, who did the first and    #
# more relevant changes in YGS2.                                #
# This was done in the year of 2021  ~~~~~(pandemic vibes)~~~~~ #
#                                                               #
# <Envinroment: Sublime Text>                                   #
# <OS: Linux (Ubuntu)>                                          #
#################################################################

import numpy as np
import pandas as pd
from pandas import DataFrame
import seaborn as sns
sns.set()
from matplotlib import cm
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('classic')
from matplotlib.gridspec import GridSpec
import argparse
import sys

## Get command-line args
userInput = argparse.ArgumentParser(description=\
	"Takes 3 .histo files produced by jellyfish histo and plots the data "
    "into output .png file.")
requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-hom', '--HomogameticFile', action='store', 
                          required=True, type=argparse.FileType('r'),
                          help='Homogametic reads .histo file produced by jellyfish histo')
requiredNamed.add_argument('-het', '--HeterogameticFile', action='store', 
                          required=True, type=argparse.FileType('r'),
                          help='Heterogametic reads .histo file produced by jellyfish histo')
requiredNamed.add_argument('-gen', '--GenomeFile', action='store', 
                          required=True, type=argparse.FileType('r'),
                          help='genome .histo file produced by jellyfish histo')

requiredNamed.add_argument('-homT', '--HomogameticTitle', action='store', 
                          required=True, type=str,
                          help='Homogametic title as should appear on plot')
requiredNamed.add_argument('-hetT', '--HeterogameticTitle', action='store', 
                          required=True, type=str,
                          help='Heterogametic title as should appear on plot')
requiredNamed.add_argument('-genT', '--GenomeTitle', action='store', 
                          required=True, type=str,
                          help='genome title as should appear on plot')

requiredNamed.add_argument('-o', '--output', action='store', 
                           required=True,
                           type=str, help='Specify the output filename')

## Import user-specified command line values.
args = userInput.parse_args()
homFile = args.HomogameticFile
hetFile = args.HeterogameticFile
genFile = args.GenomeFile
hom_title = args.HomogameticTitle
het_title = args.HeterogameticTitle
gen_title = args.GenomeTitle
outputfile = args.output

## Create and treat dataframes for each input file
files = [homFile, hetFile, genFile]
headers =[hom_title, het_title, gen_title]
dfs_names = ('female', 'male', 'genome')
dfs ={}
for dfn,file in zip(dfs_names, files):
    ## Open data file
    dfs[dfn] = pd.read_csv(file, sep=' ', usecols = [0,1], names=["k_mer", "count"])

## Plot Histo
mpl.style.use('seaborn')
mpl.rc('ytick', labelsize = 10)
mpl.rc('xtick', labelsize = 10)
mpl.rc('font',family='serif')#, serif='Tlwg Typist')

# Create plot
fig = plt.figure(figsize=(12,4))
gs = GridSpec(nrows=1, ncols=3)

plot1 = fig.add_subplot(gs[0, 0])
sns.scatterplot(ax=plot1, x="k_mer", y="count", s=20, data=np.log10(dfs[dfs_names[0]]), marker='p', palette="icefire", linewidth=0.1, alpha=0.6,legend=False)
plot1.set(xlim=(0,None))
plot1.set(ylim=(0,None))
plot1.set(xlabel='')
plot1.set(ylabel='')
plot1.set_title(headers[0], fontsize=8, pad=20, fontweight="bold")

plot2 = fig.add_subplot(gs[0, 1])
sns.scatterplot(ax=plot2, x="k_mer", y="count", s=20, data=np.log10(dfs[dfs_names[1]]), marker='p', palette="icefire", linewidth=0.1, alpha=0.6,legend=False)
plot2.set(xlim=(0,None))
plot2.set(ylim=(0,None))
plot2.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plot2.set(xlabel='')
plot2.set_ylabel('')
plot2.set_title(headers[1], fontsize=8, pad=20, fontweight="bold")

plot3 = fig.add_subplot(gs[0, 2])
sns.scatterplot(ax=plot3, x="k_mer", y="count", s=20, data=np.log10(dfs[dfs_names[2]]), marker='p', palette="icefire", linewidth=0.1, alpha=0.6,legend=False)
plot3.set(xlim=(0,None))
plot3.set(ylim=(0,None))
plot3.set(xlabel='')
plot3.set_ylabel('')
plot3.set_title(headers[2], fontsize=8, pad=20, fontweight="bold")

# Plot info and save
fig.text(0.5, -0.01, 'k-mer frequency', ha='center', fontsize=10,fontweight="bold")
fig.text(0.08, 0.5, 'Count of distinct k-mers', ha='center', va='center', fontsize=10, fontweight="bold", rotation='vertical')

plt.savefig(outputfile, dpi=300, bbox_inches='tight')
# plt.show()
plt.close()