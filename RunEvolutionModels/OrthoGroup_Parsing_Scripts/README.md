## LeAnn's Evolution Model pipeline using the output of Orthofinder, codeML, HyPhy.

#### 1. Copy results directory from OrthoFinder filled with single copy orthogroups (FASTA files). 

##### Note: If the fasta file were translated into proteins and THEN run with Orthofinder,you will have to process the files with Fasta2Phylip to backtranslate them...
#####	1.2 Use back_trans_orthogroups_gene.py
		Run Command: python3 back_trans_orthogroups_gene.py [NUC genomes dir] [PROT orthogroups genomes dir]
				[output dir name of choice]
		How it Works: It grabs all the nucleotide sequences from the original files based off the record ID's
				This allows it to fit within the rest of the pipeline as is.

##### 2. Use rename_orthogroup_genes.py
	Need: Directory filled with single copy orthogroup fasta files, directory filled with genome files used to run orthofinder
	Run Command: python3 rename_orthogroup_genes.py [genomes dir] [ortho/groups genomes dir] [output dir name of choice]
		**reminder: make sure to include the last path slash for your directories!
	How it Works: For every file in the single copy orthogroups directory, it searches the record ID's within the file against a dictionary
	created with the genomes. Within the output directory, a new dir for every orthogroup (named with the orthogroup name) is created and
	an edited sequence file is written there. The edited sequence file has the renamed record.id's for each sequence in the format:
	>speciesname_orthogroup.
		**Note: If the original files are not named as [ORTHOGROUP.fa], like 0G0001832.fa, the new directories, files, and record id edits
			will instead by named with what your file names are. However, the last three characters will be cut off (with the orthogroups,
			".fa" is cut off).
	
##### 3. translate_sequences.py
	Need: Directory filled with your gene group (or orthogroups) directories. There must be one nucleotide fasta  sequence file in each directory.
	Run Command: python3 translate_sequences.py [outer directory]
	How it Works: Removes stop codons and translates the nucleotide sequence files within each subdirectory. Creates new files named:
	"protein_[filename]"  

##### 4. curate_pipeline.sh [name for directory of gene trees]
	MUST BE IN THE DIRECTORY TO RUN
	Need: Directory filled with subdirectories, each subdirectory should be named as the ORTHOGROUP/genegroup (the result directory from rename_orthogroup_genes.py). Provide name of the directory where gene trees should be copied. 
	      Each subdirectory should only contain the nucleotide fasta file and the protein fasta file. Named like they from rename_orthogroup_genes.py and
	      translate_sequences.py (for example, OG0001832.fa and protein_OG0001832.fa)
	NOTES: Remember to make the following edits per run: change the out directory location for tree files to be copied to
	       
	Run Command: ./curatepipeline.sh
		***Note: Must be run in the Directory (still outside the subdirectories)
	How It works: Within each subdirectory it run the following commands:

			Will only start running on files with "protein_". Be careful not to create duplicates.
			Runs mafft on the protein files to align them  and output is named as "aligned_protein_[filename]" 
                        Calls Pal2Nal to back-translate the aligned protein files. Output is named "nuc_[filename]_aligned.fasta"
			Calls iqtree to create trees for "nuc_[filename]_aligned.fasta" and creates "nuc_[filename]_aligned.fasta.nwk"
			Calls Fasta2Phy to convert the "nuc_[filename]_aligned.fasta" from a fasta to a phylip of the same name.
  			Copies the "nuc_[filename]_aligned.fasta.nwk" to a new subdir (this should be edited by user) with files called "[filename]_tree.nwk"
