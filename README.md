## These are the most up-to-date versions of the pipeline for the comparitive genomics analysis for my master's thesis.

### Note: Currently, this pipeline is not stand-alone. It is formatted to run with specific software such as OrthoFinder, CODEML,
HyPhy, HMMER2GO and topGO more. Many PATHS are also specific to my computer set-up. If you would like to make use of these scripts and need help adjusting them
for your own use, please contact leann.g.ward@gmail.com.

#### InvestigateHomology

	Useful for parsing CDS/Proteome files for use with OrthoFinder.

#### RunEvolutionModels
	
	Useful for parsing OrthoFinder output to run with CODEML and HYPHY. I utilized scripts for the branch-model, branch-site model, and the FitMultiModel.
	Scritps for BUSTED and Absrel were not used for my thesis project, so they have not been tested as thoroughly.	

#### EnrichmentAnalysis

	Useful for examining CODEML results using GO annotations and Enrichment tests.

#### ModelOutputAnalysis

	Useful for examining dN/dS distributions for genes. Uses results returned by CODEML.