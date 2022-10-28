## Examining dN/dS distribution between Evolution Models

### Note: The data files used in this analysis were produced by scripts under the Evolution Models section.

#### 3ratio_branch_model_analysis.Rmd
	Description: Compare and analyze omega results for three sets of branches.
	Uses output from PAML's codeml.

#### branch_model.Rmd
	Description: Compare and analyze omega results for two sets of branches.
	Uses output from PAML's codeml.

#### enrichment_analysis_branch_sitemodel.Rmd
	Description: Run a topGO functional enrichment analysis on the results of PAML's branch-site modified model A.

#### enrichment_analysis_with_omegas.Rmd
	Description: Run a topGO functional enrichment analysis on the results of PAML's branch-model. Enrichment analysis is performed on the top and bottom 10% of omega rates.

#### exploring_aginggenes_ssp.Rmd
	Description: Run a topGO functional enrichment analysis on the results of PAML's branch-model. Enrichment analysis is performed on the omega results where the omega's for the 

#### exploring_blastedlongevitygenes_fgandbg.Rmd
	Description: Check if a pre-defined list of genes are present in results of PAML's branch-site model for genes containing significant site-specific positive selection.
	
#### exploring_bastedlongevitygenes_topandbottompercent.Rmd
	Description: Check if a pre-defined list of genes are present in results of PAML's branch-site model for genes containing significant site-specific positive selection.

#### single_copy_make_dnds_dist.R
	Description: The goal here is the see if the distribution of dN/dS values for my foreground branches 
	were significantly different than the distribution of my brackground branches.
	
	The foreground branch is determined when CODEML is run using the species tree.
	
	This will produced visualizations of the distributions. It will also return the p-values from the Komolgorov-Smirnov test.
	
