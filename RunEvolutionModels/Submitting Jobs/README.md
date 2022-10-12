## LeAnn's Evolution Model pipeline using the output of Orthofinde. For running codeML and HyPhy.

#### 1. To be used after creating control files for codeml or run files for HyPhy
	
##### 2. get_control_files_for_input.py [NOTE: This is unecessary if you use script #3. However, it can make files with the locations of the control files if you want it]
	Need: Directory filled with your gene group (or orthogroups) directories and created control files.
	Run Command: python3 get_control_files_for_input.py [outer directory] [flag name]
	Makes a list of PATHS for each control file in the directory provided that contains the flag name of interest.
	
##### 3. makelist_of_control_files.py
	Need: Directory with control .txt files created in the previous step.
	Run Command: python3 make_submission_list.py [directory name] [flag name] [output name]
	Concatenates the paths for all control files with the given string flag. This provides a list for GNU parallel to read in the next script.

##### 4. run_codeml_controlfiles_inparallel.sh 
	Need: List of control file paths created in the previous step. You can concatenate different models as needed. This simply runs codeml and a provided control file.
	       
	Run Command: ./run_codeml_controlfiles_inparallel.sh [list file name]
		
##### 5. submit_yn00_joblist_parallel.py 
	Similar to the prevous script, but instead of codeml it runs the yn00 model.
