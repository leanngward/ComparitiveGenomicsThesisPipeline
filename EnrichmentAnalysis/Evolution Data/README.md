## These scripts were written for a functional enrichment analysis of CODEML model results.

#### Note: Functional enrichment analysis requires either specific gene ID's or custom annotations. I acquire functional annotations using HMMER2GO 
on all of single-copy orthologous groups returned by OrthoFinder.

#### go_correlations_with_omega.R
	Description: The goal of this script was to analyze functional enrichment of genes with the highest and lower dN/dS values.

	For this script, you must provide the table of Orthogroups with their dN/dS results from both the branch model and the null model.
	This script will analyze which gene groups fit significantly better under the branch model than the null model.

	From there, it uses the quantile function to determine a cut-off value for the top 10% or bottom 10% of dN/dS values. Following the topGO manual,
	this script does a functional enrichment analysis of the top/bottom 10% of genes against all significant genes.

	Then, I used the enrichment barplot function to look at the top functions represented by both groups. 
	
	 


#### branchsite_go_enrichment.R
	Description: The goal of this script was to analyze functional enrichment of genes with the most and least sites under positive selection.
	Using scripts scripts under the branch-site analysis section, I had the count of sites under positive selection for each gene. 

	For this script, you must provide the counts and the table of Orthogroups with their dN/dS results from both the branch-site model and the null model.
	This script will analyze which gene groups fit significantly better under the branch-site model than the null model.

	From there, it uses the quantile function to determine a cut-off value for the top 10% or bottom 10% of sites under positive selection. Following the topGO manual,
	this script does a functional enrichment analysis of the top/bottom 10% of genes against all significant genes.

	Then, I used the enrichment barplot function to look at the top functions represented by both groups. 
	
	 