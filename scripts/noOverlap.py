#!/usr/bin/python
# coding=utf-8

'''
         This script was developed by Isabela Almeida          
 in order to remove the probes that overlap for up to n
 positions to the previous one, given a .bed file input.       
 
 This was created while I was a MSc Bioinformatics student at  
 the Universidade de São Paulo (USP/Brazil) under supervision  
 of USP PhD Prof. Maria Vibranovski and under co-supervion of  
 UFRJ PhD Prof. Bernardo de Carvalho in the year of 2020       
 ~~~~~~~~~~~~~~~~~~~~~~(pandemic vibes)~~~~~~~~~~~~~~~~~~~~~~~ 

 <Envinroment: Sublime Text>                                   
 <OS: Linux (Ubuntu)>                                          
'''

import argparse
import sys

# Get command-line args
parser = argparse.ArgumentParser(description="Takes a probes .bed file containing "
								"starting position in the second column and "
								"removes the probes that overlap for up to "
								"n positions to the previous one.")
parser.add_argument("input", help="a probes.bed file",
					default=sys.stdin,
                    type=argparse.FileType('r'),
                    nargs='?')
parser.add_argument("-n", "--noOverlap", help="number of positions that must "
					"not have probe overlap", type=int,default=20)
args = parser.parse_args()

# Open data and create output file
data = args.input.readlines()
# print "argument:", sys.argv
inputfile = sys.argv[3]
outputfile = inputfile.split(".")[0] + "F4_nO.bed"
output = open(outputfile, "w+" )
args.input.close()

# Get probes that do not overlap with the following ones for up to n positions 
output.write(data[0])

first_probe=data[0].split("\t")
previous=int(first_probe[1])

nO_probes = 1

for probe in range (1, len(data)):
	line = data[probe].split("\t")
	starts = int(line[1])

	for i in range (0,args.noOverlap+1):
		if (previous+i) == (starts): break

	if (previous+i) != (starts):
		previous = starts
		output.write(data[probe])
		nO_probes += 1

# Print info about the results to terminal.
print('noOverlap identified {0} of {1} / {2:.4f}% probes that do not overlap for up to {3} positions'.format(nO_probes, len(data), float (nO_probes) / float(len(data)) * 100, args.noOverlap))