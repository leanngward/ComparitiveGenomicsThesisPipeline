## Scripts for running different HyPhy models

### Note: These scripts were created to run on Orthogroups output from OrthoFinder. The OrthoGroup Parsing Scripts should be used to create directories for each Orthogroup and generate the correct additional files.

#### Running FitMulti HyPhy Model:
	Step 1 (Create File List for GNU parallel) Run Command: python3 create_fitmulti_inputfiles.py [directory of group directories] [tree file directory]
	
	Step 2 (run GNU parallel) Run Command: ./run_fitmodel_inparallel.sh [input file name]

#### 5. Reading HyPhy output files

		read_fitmulti_output.py
		Run Command: same as above
		How it works: reads the output from .json files with namign scheme ".fasta.FITTER.json"

#### 6. Using the FitMulti output

	Step 1: edited_mnm_site_removal.sh
	Run Command: run the .sh file in the outer directory for the gene groups subdirectories
	How it Works: Creates a list of sites that should be removed from the alignment because
			they are impacted by mnm

	Step 2: changes_sites_to_gaps.py
	Run Command: python3 changes_sites_to_gaps.py [outer directory]
	How it Works:Finds the list output by edited_mnm_site_removal.sh. Takes the aligned phylip
			file in the subdirectory and then adds gaps to the sites that should be removed
			then outputs a new phylip files starting with "edited_sites_"

		How it Works: calls the Hyphy FitMulti model and outputs a .json file with naming
				scheme ending with ".fasta.FITTER.json"
