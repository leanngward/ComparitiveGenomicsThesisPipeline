## Pipe lines to run PAML's codeml Branch-Model

### Note: These scripts were created to run on Orthogroups output from OrthoFinder. The OrthoGroup Parsing Scripts should be used to create directories for each Orthogroup and generate the correct additional files.

#### 1. check_trees1.py
	Required Package: ete3
	Need: You will need a newick tree file with a master tree of all the species in your experiment.
		***Note: Master tree will need to have specific edits. Species in the file must be writen with an underscore following the name. [speciesname]_
                        The branches/clades of interest must have a ":2" the species name. All other branch lengths (:###) must be removed.
	      You will need ete3 installed. (As of 1/27/2022 for LeAnn: source ~/.bashrc and conda activate ete3 will allow you to run it.)
	Run Command: python3 check_trees1.py [master tree] [tree directory]
	How it Works: Code double checks for spelling errors, etc., by seeing if all the species in the treefiles
        exist in the mastertree by making use of ete3. Helps make sure that both are correct.

#### 2. prune_trees2.py
	Need: The directory filled with tree files from the OrthoGroup Parsing Scripts.
	      The same master tree file used with check_trees1.py
	Run Command: python3 prune_tree2.py [tree file directory] [desired output directory name] [mastertree]
		***Notes: Do not create your output directory before running. It does this for you, and will fail if you have already made one.
			  Make sure to include a slash at the end of your paths
	How it Works: Loops through the tree files in the input directory. Creates a pruned tree file with the topology of the given supertree.
		      For use with codeml!

#### 3. parse_trees3.py
	Need: The directory filled with pruned tree from prune_trees2.py
	Run Command: python3 parse_trees3.py [pruned tree directory name] [output directory name]
			***Note: This creates the directory for you, so if it is already created then delete it or pick new output name.
	How it Works: Final curation step for codeml trees. Adds back the group/gene names to each node in the tree. Changes
                      The 2 to a $1 on the species of interest for the codeml run. See PAML manual for details about choosing branches of interest.

	UPDATE 3/02/2022: Scripts now removes ":" and "1" for a cleaner output treefiles.
	UPDATE NEEDED: If there is '3' for a second foreground species, change the 3 to a $2


#### 4. create_branchmodel.py
	Need: You will need a directory full of gene group subdirectories and a directory of you final tree files. These should be created by the previous scripts.
	Run Command: python3 create_codeml_files.py [directory of group directories] [tree file directory] [Group Name Flag]
	How it Works: This file uses the subdirectory groupname to create codeml control files for both an alternative and null branch model. It also creates .sh files to submit all the control files to the job queue using the "qsub" command. All the files it creates will be in the main directory, outside the group directories. If all the files don't completely run, the created job files can still be run individually.
			***Note: re-running script will create duplicate files. Delete these before running again.
	The group name flag allows the addition of more strings to the output files (and paths in the control file) to allow for multiple models to be run.
	
	Important Parameters:
		Alt:
			model = 2
			Nssites = 0
		Null:
			model = 0
			Nssites = 0
	Note: These are only suggestions, parameters should be determined from the PAML manual.

#### 5. See "Submitting Jobs" Section under Thesis Project Pipeline/RunEvolutionModels /

#### 6. read_codeml_output.py
	Need: Directory of gene group subdirectories and the NAME of the output file you want. It is tab delimited
              so a .tsv file is recommended but not required.
        Run Command: python3 read_codeml_output.py [runfile directory] [results output file name]
	How it Works: This code systematically checks in the subdirectories for the codeml outputs. If it finds one, it stores the appropraite lnl and omega
			values depending on if it is a null or alternative output file. If it does not find one, it raises a flag to the output screen.
			It also creates a list called "undonelist.txt" (it is the same everytime and will replace itself when you run it again). The list consists
			of the appropriate run statement for the file you want, along with the qsub, so that you can directory copy it and run it again (once errors are solved).
                       You can continue running this output until you have your complete results.
			***Note: This is designed for how mine are named by the previous scripts: ending in _null.out or _alt.out
			for the null and alternative models
	

#### Additional: delete_extras.py
	Run Command: python3 delete_extras.py [out directory filled with subdirectories] [string flag to be deleted]
	How it Works: This python file is for my specific file structure. It's helpful for automaticaly deleting
			files created in my gene subdirectories that I don't want anymore. For example,
			if I wanted to delete all the codeml results that I created, this would do it automatically.



