import os
import sys
import re

#Provide directory that contains the files of interest
dirName = sys.argv[1]

#Provide the string pattern for the file names that you want to record
stringFlag = sys.argv[2]

#Provide the output file name for the list of files
outName = sys.argv[3]
outFile = open(outName,'w')

for filename in os.listdir(dirName):
	f = os.path.join(dirName,filename)
	if os.path.isfile(f):
		if re.search(stringFlag,f):
			f = f.strip()
			outFile.write(f+"\n")


outFile.close()
