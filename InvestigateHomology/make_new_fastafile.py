#!/usr/bin/env/python3

import Bio
from Bio import SeqIO
from Bio.Seq import Seq
import sys
import re

def search(arr, x):
	for i in range(len(arr)):
		if arr[i] == x:
			return True;

	return  False



faa_name = sys.argv[1] #Put the name of the fasta file in the python call
ogp = open(faa_name, "r") #Open and read the protein file

tx_name = sys.argv[2]
longlist = open(tx_name, "r")

outname = sys.argv[3]
outfile = open(outname, "w+")

longesttx = []

#Append the names of the longest transcripts to a list
for lines in longlist.readlines():
	lines = lines.rstrip()
	longesttx.append(lines)

#for lines in ogp.readlines():
#	if re.search(">", lines):
#		firstpos = re.search(">", lines).end()
#		lastpos = re.search(" ", lines).end()
#		print(lines[firstpos:lastpos])

#Read the protein fasta file
for record in SeqIO.parse(ogp, "fasta"):
	name = record.id
#	print(name)
	if search(longesttx, name):
		outfile.write(">")
		outfile.write(str(record.description))
		outfile.write("\n")
		outfile.write(str(record.seq))
		outfile.write("\n")

ogp.close()
longlist.close()
outfile.close()
