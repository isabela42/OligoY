#!/usr/bin/python
# coding=utf-8

'''
         This script was developed by Isabela Almeida          
 in order to get the first report of each probe found in the
 given .bed file input.       
 
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
								"probes in the fourth column and "
								"returns the first report of each probe in that file.")
parser.add_argument("input", help="a probes.bed file",
					default=sys.stdin,
                    type=argparse.FileType('r'),
                    nargs='?')

args = parser.parse_args()

# Open data and create output file
data = args.input.readlines()
# print "argument:", sys.argv
inputfile = sys.argv[1]
outputfile = inputfile.split(".")[0] + "F5_nR.bed"
output = open(outputfile, "w+" )
args.input.close()

# Get the first report of each probe
out_probes=[]
times_probes=[]

#first_probeLine=data[0].split("\t")
#out_probes.append(data[0])

nR_probes = 0

for probe in range (0, len(data)):
	thisLine = data[probe].split("\t")
	thisProbe = thisLine[3]
	unique=True

	for eachProbe in range (0,len(out_probes)):
		outLine = out_probes[eachProbe].split("\t")
		outProbe = outLine[3]

		if thisProbe == outProbe:
			times_probes[eachProbe]  += 1
			unique=False
			break

	if unique == True:
		out_probes.append(data[probe])
		times_probes.append(1)
		output.write(data[probe])
		nR_probes += 1

oneTime = TwotoTen = ElevenToTwenty = TwentyoneToThirty = ThirtyoneToFourty = FourtyPlus = 0 
for number in range (1, max(times_probes)+1): #se der erro, tira esse +1
	times = times_probes.count(number)
	if number == 1: oneTime = times
	elif 1 < number <=10: TwotoTen += times
	elif 10 < number <=20: ElevenToTwenty += times
	elif 20 < number <=30: TwentyoneToThirty += times
	elif 30 < number <=40: ThirtyoneToFourty += times
	elif number > 40: FourtyPlus += times

print('nonRedudant identified {0} first reports out of the {1} probes / {2:.4f}%'.format(nR_probes, len(data), float (nR_probes) / float(len(data)) * 100))
print('{0} probes occured only once / {1:.4f}% of non redundant probes and {2:.4f}% of input probes'.format(oneTime, float (oneTime) / float(len(out_probes)) * 100, float (oneTime) / float(len(data)) * 100))
print('{0} probes occured from 2 to 10 times / {1:.4f}% of non redundant probes and {2:.4f}% of input probes'.format(TwotoTen, float (TwotoTen) / float(len(out_probes)) * 100, float (TwotoTen) / float(len(data)) * 100))
print('{0} probes occured from 11 to 20 times / {1:.4f}% of non redundant probes and {2:.4f}% of input probes'.format(ElevenToTwenty, float (ElevenToTwenty) / float(len(out_probes)) * 100, float (ElevenToTwenty) / float(len(data)) * 100))
print('{0} probes occured from 21 to 30 times / {1:.4f}% of non redundant probes and {2:.4f}% of input probes'.format(TwentyoneToThirty, float (TwentyoneToThirty) / float(len(out_probes)) * 100, float (TwentyoneToThirty) / float(len(data)) * 100))
print('{0} probes occured from 31 to 40 times / {1:.4f}% of non redundant probes and {2:.4f}% of input probes'.format(ThirtyoneToFourty, float (ThirtyoneToFourty) / float(len(out_probes)) * 100, float (ThirtyoneToFourty) / float(len(data)) * 100))
print('{0} probes occured more than 40 times  / {1:.4f}% of non redundant probes and {2:.4f}% of input probes'.format(FourtyPlus, float (FourtyPlus) / float(len(out_probes)) * 100, float (FourtyPlus) / float(len(data)) * 100))
print('This values add up to the original number of probes found in the input .bed file ({0})'.format(sum(times_probes)))