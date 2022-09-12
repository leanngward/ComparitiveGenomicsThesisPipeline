import sys
import re


inName = sys.argv[1]
inFile = open(inName,'r')

outName = sys.argv[2]
outFile = open(outName,'w')
sitelist = []

outFile.write("GroupName\tTotal.Sig.Sites\n")
for line in inFile.readlines():
	count = 0
	line = line.strip()
	sitelist = line.split("\t")
	groupname = sitelist[0]
	for item in sitelist:
		item = item.strip(" ")
		item = item.strip(",")
		if re.search("\*",item):
#			AApos = re.search("\D",item).start()
			#AA = item[AApos+1])
#			sitecount = item[0:AApos]
			count += 1
	outFile.write(groupname+"\t"+str(count)+"\n")
