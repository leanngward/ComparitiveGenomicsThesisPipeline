import os
import sys


def codeml_branches_runner(codeml_prefix,seq_file,treefile,outfile,gene_name,groupName):
    import subprocess
    #create the control file for codeml to run
    control_file_name = codeml_prefix  + gene_name + '_'+groupName+'_branchsite_null.ctl'
    with open(control_file_name, 'w') as fh: 
        fh.write('seqfile = ' + seq_file +'\n')
        fh.write('treefile = '+ treefile +'\n')
        fh.write('outfile = ' + outfile +'\n')
        fh.write('noisy = 9 \n')
        fh.write('verbose = 0 \n')
        fh.write('runmode = 0 \n')
        fh.write('seqtype = 1 \n')
        fh.write('CodonFreq = 2 \n')
        fh.write('clock = 0 \n')
        fh.write('aaDist = 0 \n')
        fh.write('model = 2 \n')
        fh.write('NSsites = 2 \n')
        fh.write('Mgene = 0 \n')
        fh.write('fix_kappa = 0 \n')
        fh.write('kappa = 2 \n')
        fh.write('fix_omega = 1 \n')
        fh.write('omega = 1 \n')
        fh.write('fix_alpha = 1 \n')
        fh.write('alpha = 0 \n')
        fh.write('Malpha = 0 \n')
        fh.write('ncatG = 8  \n')
        fh.write('getSE = 0 \n')
        fh.write('RateAncestor = 1 \n')
        fh.write('Small_Diff = .5e-6 \n')
        fh.write('method = 0 \n')


#Provide directory of the gene name subdirectories
dirName = sys.argv[1]

#Provide location of edited tree files
treeDir = sys.argv[2]

#Provide group name flag
groupName = sys.argv[3]

flag = ''
mnmFlag = False
#Provide True or False for whether you are tested files with removed sites affected by multi-nucleotide mutations
try:
	flag = sys.argv[4]
	if flag == "T":
		mnmFlag = True
	else:
		mnmFlag = False
except IndexError:
	pass

dirList = os.listdir(dirName)
geneDirName = ''
geneName = ''

for items in dirList:
	geneDirName = dirName + items
	if(os.path.isdir(geneDirName)):
		geneName = items
		treeFileName = treeDir+geneName+"_edited_tree.nwk" #LOCATION OF EDITED TREEFILES
		if mnmFlag == False:
			seqFileName = dirName+geneName+"/"+"nuc_"+geneName+"_aligned.phy"
			altoutFileName = dirName+geneName+"/"+geneName+"_"+groupName+"_null_branchsite.out"
		else:
			seqFileName = dirName+geneName+"/edited_sites_nuc_"+geneName+"_aligned.phy"
			altoutFileName = dirName+geneName+"/"+geneName+"_"+groupName+"_mnm_nullbranchsite.out"

		codeml_branches_runner(dirName,seqFileName,treeFileName,altoutFileName,geneName,groupName)
