import os
import sys
def codeml_run(ctrl_file,scriptfile_b,gene_name):
	import subprocess

	with open(scriptfile_b, 'w') as fh:
		fh.write(ctrl_file)
		fh.write(' \n')


#Provide directory of the gene name subdirectories
dirName = sys.argv[1]

#provide the model name for unique files name
flag = sys.argv[2]


dirList = os.listdir(dirName)
geneDirName = ''
geneName = ''


for items in dirList:
	geneDirName = dirName + items
	if(os.path.isdir(geneDirName)):
		geneName = items
		controlFileName = dirName+geneName+"_"+flag+".ctl"
		nullScriptName = dirName+geneName+"_"+flag+".txt"
		codeml_run(controlFileName,nullScriptName,geneName)
