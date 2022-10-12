import os
import sys
def fitmodel_run(prefix,seq_file,treefile,outfile,scriptfile_b,gene_name):
        import subprocess

        with open(scriptfile_b, 'w') as fh:
                fh.write(seq_file)
                fh.write(",")
                fh.write(treefile)
                fh.write(",")
                fh.write(outfile)
                fh.write("\n")



#Provide directory of the gene name subdirectories
dirName = sys.argv[1]

#Provide location of edited tree files
treeDir = sys.argv[2]


dirList = os.listdir(dirName)
geneDirName = ''
geneName = ''


for items in dirList:
        geneDirName = dirName + items
        if(os.path.isdir(geneDirName)):
                geneName = items
                seqFileName = dirName+geneName+"/"+"nuc_"+geneName+"_aligned.fasta"
                treeFileName = treeDir+geneName+"_edited_tree.nwk" #LOCATION OF EDITED TREEFILES
                nulloutFileName = dirName+geneName+"/"+geneName+"_fitmodel.out"
                nullScriptName = dirName+geneName+"_fitmodel.txt"
                fitmodel_run(dirName,seqFileName,treeFileName,nulloutFileName,nullScriptName,geneName)
