import sys
import subprocess
#Provide the name of list of the control files
inName = sys.argv[1]

out = sys.argv[2]

outfile = open(out, 'w')
outfile.write("#!/bin/bash \n")
outfile.write("#$ -pe threads 16 \n")
outfile.write("module load PAML/4.9 \n")
#activate conda source so you can run GNU parallel
outfile.write("source /mnt/home/software/anaconda/anaconda2/bin/activate ete3 \n")
outfile.write("filename=\""+inName+"\" \n")
outfile.write("parallel --will-cite yn00 {1}")

outfile.close()
