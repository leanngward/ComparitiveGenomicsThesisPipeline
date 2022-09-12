import sys
import re
from Bio import SeqIO
mapname = sys.argv[1]
genomename = sys.argv[2]
outname = sys.argv[3]


gfile = open(genomename,'r')
chrom = {}
for record in SeqIO.parse(gfile,"fasta"):
	recid = str(record.id)
	recid = recid.strip()
	recdes = str(record.description)
	chromid = ''
	if re.search("BANY",recdes):
		start = re.search("BANY",recdes).start()
		chromid = recdes[start:start+9]
	chrom[recid] = chromid

#print(chrom)

mapfile = open(mapname,'r')
outfile = open(outname,'w')
for lines in mapfile.readlines():
	if lines[0] != "#":
		linelist = lines.split("\t")
		nwchrom = linelist[4]
		nwchrom = nwchrom.strip()
		bchrom = chrom[nwchrom]
		lines = linelist[0]+"\t"+linelist[1]+"\t"+linelist[2]+"\t"+linelist[3]+"\t"+linelist[4]+"/"+bchrom+"\t"+linelist[5]+"\t"+linelist[6]
		outfile.write(lines)
	else:
		outfile.write(lines)

gfile.close()
mapfile.close()
outfile.close()
