#!/usr/bin/env python

import argparse
import os

#####################################
# Description of script functionality   #
#               #
#####################################

__author__ = "Patrick Barth" 
__email__ = "patrick.barth@computational.bio.uni-giessen.de"

parser = argparse.ArgumentParser()
parser.add_argument('--sequences','-s', 
					help='Sequences to count nucleotides')
parser.add_argument('--output','-o', 
					help='File to write output to')
parser.add_argument('--output_percent','-p', 
					help='File to write percentage output to')
parser.add_argument('--length','-l', 
					help='Length of sequences', default='21')
parser.add_argument('what_shall_i_write_here', nargs=argparse.REMAINDER)
args = parser.parse_args()

#######################
#######################
###    variables    ###
#######################
#######################



########################
########################
###    parameters    ###
########################
########################
length = int(args.length)


##########################################################################################################################################

###################
###################
### main script ###
###################
###################


def main():

	nucleotides_alignments = count_nucleotides(args.sequences)

	if args.output_percent:
		nucleotides_percent = [None] * length
		total_count = sum(nucleotides_alignments[0].values())
		for pos in list(range(0,length)):
			nucleotides_percent[pos] = {"A":0,"C":0,"G":0,"T":0,"N":0}
			for nuc,count in nucleotides_alignments[pos].items():
				nucleotides_percent[pos][nuc] = count/total_count * 100 if not total_count == 0 else 0
			

	out_file = open(os.path.realpath(args.output), 'w')

	out_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
		"Alignments_" + str(length),
		"A",
		"C",
		"G",
		"T",
		"N"
	))
	counter = 1
	for pos in nucleotides_alignments:
		out_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
			str(counter),
			pos["A"],
			pos["C"],
			pos["G"],
			pos["T"],
			pos["N"],
		))
		counter += 1

	if args.output_percent:
		out_file = open(os.path.realpath(args.output_percent), 'w')

		out_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
			"percentage_" + str(length),
			"A",
			"C",
			"G",
			"T",
			"N"
		))
		counter = 1
		for pos in nucleotides_percent:
			out_file.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (
				str(counter),
				pos["A"],
				pos["C"],
				pos["G"],
				pos["T"],
				pos["N"],
			))
			counter += 1

        

#######################
#######################
###    functions    ###
#######################
#######################

def count_nucleotides( file ):
	nucleotides = [None] * length
	for i in list(range(1,length+1)):
		nucleotides[i-1] = {"A":0,"C":0,"G":0,"T":0,"N":0}

	with open ( file, 'rt' ) as parse_file:
		for line in parse_file:
			line = line.strip()
			split_line = list(line)

			counter = 0
			for nuc in split_line:
				nucleotides[counter][nuc] += 1
				counter += 1
	return nucleotides



##########################
### starts main script ###
##########################
main()
