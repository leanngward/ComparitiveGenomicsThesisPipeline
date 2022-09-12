##Using this for when the CDS list is longer than the updated CDS file. I wanna
##see which CDS ID's are missing
##Update 6/21/2022 I think the CDS missing were just duplicates. so this script finds duplicates


import re
import sys
from Bio import SeqIO
cds_input_name = sys.argv[1]
faname = sys.argv[2]
ogp = open(faname,'r')
def duplicates(seq):
	seen = set()
	seen_twice = [x for x in seq if x in seen or (seen.add(x) or False)]
	return list( seen_twice )

cds_input = open(cds_input_name, 'r')
cdslist = []
count = 0
for items in cds_input.readlines():
	items = items.strip()
	cdslist.append(items)

seenlist = duplicates(cdslist)
print(len(seenlist))


## for checking if there are no duplicates
#reclist = []
#for record in SeqIO.parse(ogp, "fasta"):
#	name = record.id
#	reclist.append(name)

#for item in cdslist:
#	if item not in reclist:
#		print(item)

cds_input.close()

