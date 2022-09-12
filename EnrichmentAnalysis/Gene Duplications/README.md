## Performing Functional Enrichment Analysis for Paralogs 
### Paralogs examined with these script came from OrthoFinder software (LINK)

#### Note: Functional enrichment analyses with topGO require certain types of IDs OR custom annotations. I retrieved custom GO annotations for all of my species' genes using HMMER2GO. (LINK)

#### 1. grab_geneids_fromduplicationslist.R
	Description: This R script allowed me to get a list of genes that resulted from	duplication events found by Orthofinder.
	You must provided the two nodes of interest where you want to return genes from. It will grab gene id's that were
	duplication at these nodes that have greater than 50% support values. 
	
	Once you have all the genes, there are two sections for parsing out genes from specific species. All of my genes id's contained
	the species name as a result of Orthofinder. The first section will exclude 1 specific species. The second section will include 1 specific species.
	The "non-species" string will include/exclude gene id's with that string present.

	The last line of code will write an output file with only the genes of interest.

#### 2. BP/CC/MFduplication_groups_analysis.R
	Description: These threes script are the same with the exception of which gene ontology terms it focuses on (Biological Process (BP), Cellular Components (CC),
	or Molecular Functions(MF). You can provide the names of files (created by the previous script). In my case, I included every species of interest so that I could
	compare the results at the end.

	If you want to add or subtract a species group, simply copy a line and change the name or comment out one of the repeated lines.

	This script determines which genes, among all species genes (in my case, there were 10 total species but 4 species of interest) were amongst genes duplicated
	in each group. Then, I used the topGO object function on each group to determine with genes functions had significant enrichment among all background genes.
	

	The script also run topGO statistical analysis using runTest. Using this information, I formatted the data to find significantly upregulated and downregulated
	functions among the gene duplications against the background of all species genes. I used the compareCluster (LINK) to visualize the differences in functions
	represented by all the groups.