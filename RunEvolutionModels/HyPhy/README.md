## Scripts for running different HyPhy models

### Note: These scripts were created to run on Orthogroups output from OrthoFinder. The OrthoGroup Parsing Scripts should be used to create directories for each Orthogroup and generate the correct additional files.

#### 1. check_trees1.py
	Required Package: ete3
	Need: You will need a nwk file with a master tree of all the species in your experiment.
		***Note: Master tree will need to have specific edites. Species in the file must be writen with an underscore following the name. [speciesname]_
                        The branches/clades of interest must have a ":2" the species name. All other branch lengths (:###) must be removed.
	Run Command: python3 check_trees1.py [master tree] [tree directory]
	How it Works: Code double checks for spelling errors, etc., by seeing if all the species in the treefiles
        exist in the mastertree by making use of ete3. Helps make sure that both are correct.

#### 2. prune_trees2.py
	Need: The directory filled with tree files from the curate_pipeline.sh that has been checked.
	      The same master tree file used with check_trees1.py
	Run Command: python3 prune_tree2.py [tree file directory] [desired output directory name] [mastertree]
		***Notes: Do not create your output directory before running. It does this for you, and will fail if you have already made one.
			  Make sure to include a slash at the end of your paths
	How it Works: Loops through the tree files in the input directory. Creates a pruned tree file with the topology of the given supertree.
		      For use with codeml!

#### 3. parse_tree_forhyphy.py
	Need: A directory of pruned trees with a '2' labeled on the branches of interest.
	Run Command: python3 parse_tree_forhyphy.py [pruned tree directory name] [output directory name]
	How it Works: It edits the nodes of interest with the label {Foreground}. "Foreground" must be
	explicity described when running a Hyphy model using --branches to define the hypothesis being tested.

#### 4. Running Different HyPhy Models
	Run Command: python3 [file name from below] [directory of group directories] [tree file directory]
	Options:
		run_absrel_files.py
		Run Command: Same as above
		How it Works: calls the Hyphy absrel model and outputs a .json file with naming
				scheme ending with "_absrel_model.out"
		run_busted_files.py
		Run Command: Same as above
		How it Works: calls the Hyphy BUSTED model and outputs a .json file with naming
				scheme ending with ".fasta.BUSTED.json"

		run_fitmulti_model.py
		Run Command: Same as above

#### 5. Reading HyPhy output files

		read_absrel_output.py
		Run Command: python3 read_absrel_output.py [outer directory] [output file name]
		How it Works: gather information from .json file and puts into a tab-separated output file

		read_busted_output.py
		Run Command: python3 read_busted_output.py [outer directory] [output filename]
		How it Works: read about from .json files with naming scheme ".fasta.BUSTED.json"
		6/24/2022: noted that some outputs did not have the "constrained model" ... may need to figure out why

		read_fitmulti_output.py
		Run Command: same as above
		How it works: reads the output from .json files with namign scheme ".fasta.FITTER.json"

#### 6. For use with FitMultiModel:

##### edited_mnm_site_removal.sh
	Run Command: run the .sh file in the outer directory for the gene groups subdirectories
	How it Works: Creates a list of sites that should be removed from the alignment because
			they are impacted by mnm

##### changes_sites_to_gaps.py
	Run Command: python3 changes_sites_to_gaps.py [outer directory]
	How it Works:Finds the list output by edited_mnm_site_removal.sh. Takes the aligned phylip
			file in the subdirectory and then adds gaps to the sites that should be removed
			then outputs a new phylip files starting with "edited_sites_"

		How it Works: calls the Hyphy FitMulti model and outputs a .json file with naming
				scheme ending with ".fasta.FITTER.json"
