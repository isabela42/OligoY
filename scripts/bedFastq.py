#!/usr/bin/python

'''
         This script was developed by Isabela Almeida          
 in order to get a .fastq file out of a .bed one. It was       
 originally designed for transforming probes.bed files into    
 probes.fastq from bash (bellow), and it should work normally 
 if the .bed is in the following format:
        	 "Chr"	"start"	"end"	"seq"	"etc"

 Bash version:

 cat $input | while read line; do chr=`echo $line | \
 cut -d " " -f1`; start=`echo $line | cut -d " " -f2`; \
 end=`echo $line | cut -d " " -f3`; seq=`echo $line | \
 cut -d " " -f4`; seqSize=`echo ${#seq}` ; \
 seqSize=$((seqSize+1)); echo @${chr}:${start}-${end} \
 >> fastq ; echo ${seq} >> fastq; echo + >> fastq; \
 seq -s~ $seqSize| tr -d '[:digit:]' >> fastq ; done

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
parser = argparse.ArgumentParser(description="Takes a .bed file structured as: "
								"'Chr'	'start'	'end'	'seq'	'etc'"
								" and returns the equivalent .fastq file."
								" (Note: Designed for probes.bed files)")
parser.add_argument("input", help="a .bed file (Specifically, designed for probes.bed files",
					default=sys.stdin,
                    type=argparse.FileType('r'),
                    nargs='?')

args = parser.parse_args()

# Open data and create output file
data = args.input.readlines()
# print "argument:", sys.argv
inputfile = sys.argv[1]
outputfile = inputfile.split(".")[0] + ".fastq"
output = open(outputfile, "w+" )
args.input.close()

# Generate fastq info for each sequence and print to output file
for seq in range (0, len(data)):
	thisSeq = data[seq].split("\t")
	chromosome = thisSeq[0]
	start = thisSeq[1]
	end = thisSeq[2]
	sequence = thisSeq[3]
	seqSize = len(sequence)
	tio="~"
	tio=tio.replace('~', '~'*seqSize, 1)

	output.write('@' + chromosome + ':' + start + '-' + end + '\n')
	output.write(sequence + '\n')
	output.write("+\n")
	output.write(tio + '\n')

# Print info about the results to terminal.
print('BedFastq.py generated a .fastq file using the inicial {0} sequences.bed'.format(len(data)))