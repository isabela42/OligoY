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

import pandas as pd
from pandas import DataFrame
import seaborn as sns
sns.set()
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.use('classic')
from matplotlib.gridspec import GridSpec
import argparse
import sys

## Get command-line args
userInput = argparse.ArgumentParser(description=\
	"Takes a .ygs file produced by YGS2 and plots the data into output "
    ".png file according VSCK and PVSCUK thresholds provided.")
requiredNamed = userInput.add_argument_group('required arguments')
requiredNamed.add_argument('-f', '--file', action='store', 
                          required=True, type=argparse.FileType('r'),
                          help='.ygs file produced by YGS2')
requiredNamed.add_argument('-vsck', '--vsck_threshold', action='store',        
                          required=True, type=int,
                          help='Valid Single Copy k-mers threshold')
requiredNamed.add_argument('-pvscuk', '--pvscuk_threshold', action='store',        
                          required=True, type=int,
                          help='Percentage of Valid Single Copy Unmatched k-mers threshold')
requiredNamed.add_argument('-o', '--output', action='store', 
                           required=True,
                           type=str, help='Specify the output filename')

## Import user-specified command line values.
args = userInput.parse_args()
ygsFile = args.file
no_vsck = args.vsck_threshold
per_vscuk = args.pvscuk_threshold
outputfile = args.output

## Open data file
df = pd.read_csv(ygsFile, sep='\t', usecols = [8,10])

## Add Classfication info
classification = []
df = df.values.tolist()
for gi in range(0,len(df)):
    if df[gi][1] == 'NA': df[gi][1] = 0.0
    else: df[gi][1] = pd.to_numeric(df[gi][1])

    if df[gi][0] >= no_vsck and df[gi][1] >= per_vscuk:
        classification.append('Inferred Chr Y')
    elif df[gi][0] >= no_vsck and df[gi][1] < per_vscuk: classification.append('Inferred X/A')
    elif df[gi][0] < no_vsck: classification.append('Lower than VSCK Threshold')
df = DataFrame(df,columns=['VSC_K','P_VSC_UK'])
df['Classification'] = classification
df = df.sort_values('Classification')

## Plot Style
# Colors: BuPu_r, CMRmap, colorblind,GnBu_r,PRGn_r,YlGnBu_r,afmhot,cubehelix,gist_earth,gist_heat,gist_stern,icefire
# https://medium.com/@morganjonesartist/color-guide-to-seaborn-palettes-da849406d44f

mpl.style.use('seaborn')
mpl.rc('ytick', labelsize = 7)
mpl.rc('xtick', labelsize = 7)
mpl.rc('font',family='serif')

## Plot the data
fig = plt.figure(figsize=(14,3.5))
gs = GridSpec(nrows=1, ncols=4)

firstPlot = fig.add_subplot(gs[0, 0])
sns.scatterplot(ax=firstPlot, x="P_VSC_UK", y="VSC_K", s=20, data=df, marker='p', hue="Classification", palette="icefire", alpha=0.3)#,legend=False)
firstPlot.set(xlim=(-2,102))
firstPlot.set(ylim=(-2,None))
firstPlot.set(xlabel='')
firstPlot.set_ylabel('Valid Single Copy Kmers (no. VSCK)', fontsize=8)
handles, labels = firstPlot.get_legend_handles_labels()
plt.legend(bbox_to_anchor=(0.5,-0.04),handles=handles[1:], labels=labels[1:], loc="upper center", ncol=3, scatterpoints=1, markerscale=1, fontsize=6)

secondPlot = fig.add_subplot(gs[0, 1])
sns.scatterplot(ax=secondPlot, x="P_VSC_UK", y="VSC_K", s=20, data=df, marker='p', hue="Classification", palette="icefire", alpha=0.3,legend=False)
secondPlot.set(xlim=(-2,102))
secondPlot.set(ylim=(-2,1000))
secondPlot.set(xlabel='')
secondPlot.set_ylabel('Valid Single Copy Kmers (no. VSCK)', fontsize=7,labelpad=0.05)
secondPlot.tick_params(labelsize=4)

thridPlot = fig.add_subplot(gs[0, 2])
sns.scatterplot(ax=thridPlot, x="P_VSC_UK", y="VSC_K", s=20, data=df, marker='p', hue="Classification", palette="icefire", alpha=0.3,legend=False)
thridPlot.set(xlim=(-2,102))
thridPlot.set(ylim=(-0.5,50))
thridPlot.set(xlabel='')
thridPlot.set_ylabel('Valid Single Copy Kmers (no. VSCK)', fontsize=7, labelpad=0.05)
thridPlot.axhline(20, ls='--',linewidth=3, color='brown')
thridPlot.tick_params(labelsize=4)

## Remove information with no validation
# df = pd.read_csv(ygsFile, sep='\t', usecols = [8,10])

# df = df[df.P_VSC_UK != 'NA']
# df = df.values.tolist()
# for gi in range(0,len(df)):
#     df[gi][1] = pd.to_numeric(df[gi][1])
    
#     df[gi][0] = pd.to_numeric(df[gi][0])
# df = DataFrame(df,columns=['VSC_K','P_VSC_UK'])

fourthPlot = fig.add_subplot(gs[0, 3])
sns.distplot(df["P_VSC_UK"],kde=False)
sns.rugplot(df["P_VSC_UK"], color='black', ax=fourthPlot)
fourthPlot.set(xlim=(0,100))
fourthPlot.set(xlabel='')
fourthPlot.set_ylabel('Number of valid scaffolds', fontsize=7,labelpad=0.05)
fourthPlot.tick_params(labelsize=4)

## Plot info and save
fig.text(0.5, 0.01, '% kmers unmatched by female reads', ha='center', fontsize=10)
plt.savefig(outputfile, dpi=300)
plt.close()
# plt.show()