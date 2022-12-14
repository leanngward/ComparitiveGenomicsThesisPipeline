import os
import sys
def null_codeml_branches_runner(codeml_prefix,seq_file,treefile,outfile,scriptfile_b,gene_name):
    import subprocess
    #create the control file for codeml to run
    control_file_name = codeml_prefix  + gene_name + '_branch_null.ctl'
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
        fh.write('model = 0 \n')
        fh.write('NSsites = 0 \n')
        fh.write('Mgene = 0 \n')
        fh.write('fix_kappa = 0 \n')
        fh.write('kappa = 2 \n')
        fh.write('fix_omega = 0 \n')
        fh.write('omega = 0.4 \n')
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

dirList = os.listdir(dirName)
geneDirName = ''
geneName = ''


for items in dirList:
	geneDirName = dirName + items
	if(os.path.isdir(geneDirName)):
		geneName = items
		seqFileName = dirName+geneName+"/"+"nuc_"+geneName+"_aligned.phy"
		treeFileName = treeDir+geneName+"_edited_tree.nwk" #LOCATION OF EDITED TREEFILES
		nulloutFileName = dirName+geneName+"/"+geneName+"_branch_null.out"
		nullScriptName = dirName+geneName+"_branch_null_runcod.sh"
		##note: nullScriptName is not currently being used
		null_codeml_branches_runner(dirName,seqFileName,treeFileName,nulloutFileName,nullScriptName,geneName)
