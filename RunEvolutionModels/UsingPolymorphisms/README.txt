## These scripts were used to run the McDonald-Kreitman test. It requires a full genome file, a gff annotation file, and sites/frame information for genes of interest.

### 2v_cntsnps.py
Version Finished: 8/29/2022
Run Command: python3 2v_cntsnps.py [Population File .VCF] [Gene Mapping File] [Full Genome FASTA File] [Output File Name]
Description: This script counts non-synonymous and synonymous SNP's for genes of interest.
	1. VCF FILE:
		Needs VCF Format Type: 4. Can be compresssed or uncompressed.
	2. Gene mapping File:
		For my research, this script is formatted to read a gene mapping file with the following tab-separated columns:
		GroupID \t GeneID \t TranscriptID \t CdsID \t Chromosome/Scaffold Name \t CDS Position \t CDS Reading Frame

		GroupID: I retrieved genes from Orthogroups created by Orthofinder, so I provided Orthogroup ID's in the mapping file.
			 There are no specific requirements for this column. The information is written into the output file. <br>
		GeneID:  The GeneID was defined as the Gene ID from the annotation file before it was provided to the script. The script
			 does not use the information for anything except the output file. <br>			
		TranscriptID: The TranscriptID was defined by the annotation file before it was provded to the script. The script does not use
			      the information for anything except the output file. For gff files, these TranscriptID and GeneID are often the same. <br>
		CdsID: This ID must matched the Record ID's in the genome FASTA file.
		Chromosome/Scaffold Name: This ID must match the chromosome ID's the first column of the VCF File.
		CDS Position: Start and end positions are separated by ':' symbol. If there are multiple stop and start positions for a gene,
			      they should be separated by a ';'.
			      Example: 1350:2000;2500:3000;3100:3300
		Reading Frame: There should be a reading frame (0,1,2) for each start and stop position. For example, if there are four start/stops
			       then the frame should be written as 0;0;0;0
			       See Ensembl's gff/gtf description for more details about reading frames. 
	3. Output File:
		The output file will write the following tab-separated columns:
		GroupID \t Gene-ID \t Chromosome \t Synonymous-SNPs \t Non-Synonymous-SNPs

### findposition.py
Run Command: python3 findposition.py [Full Genome FASTA File]
Description: I created this as a quick assisting script. The positions of the nucleotides in the scaffolds/chromosomes MUST match the positions in the VCF file. Because of the way VCF files are made from the reference, they should match. I made this script to validate that information. Within the script, you can build a list of specfic number positions and the script will print what nucleotides are there.

### mapbtochrom.py
Run Command: python3 mapbtochrom.py [Mapping File] [Full Genome FASTA File] [Outfile Name]
Description: The 2v_countsnps.py scripts uses the chromosome/scaffold ids at the BEGINNING of the FASTA record. It is assuming there is a space (For Example: >Chrom1 loc12_infomation because that BioPython will use that as the "record.id" with the SeqIO function. For my data, the genome annotation had a different ID than the chromosome id's in the VCF files. 
So, this allowed me to created a map file for the two ID's. 

2v_countsnp.py looks for "/" in the mapfile if there are two ID's to go between. Otherwise, it ignores that and just uses the record.id as the chromosome name. 
