import sys
import re
filename = sys.argv[1]
outputname = sys.argv[2]
infile = open(filename,'r')
outfile = open(outputname,'w')

outfile.write("#"+"filename"+"-ids\n")

for line in infile.readlines():
	if re.search("#",line):
		continue
	else:
		line = line.split("\t")
		outfile.write(line[1]+"\n")
infile.close()
outfile.close()
