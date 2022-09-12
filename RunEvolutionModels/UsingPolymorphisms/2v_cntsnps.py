#Note: I am essentially starting from scratch here for this script. Although some sections
#might be copied from previous iterations.

#Note: I will need a version of this script where I need the chromosome information.
#Since validation is easier for B. anynana files. I will start with NOT needing the chromosome.

import sys
import re
import numpy as np
import Bio.Data.CodonTable
from Bio.Seq import Seq
from Bio import SeqIO
import gzip
import datetime

#############################################################
#FUNCTIONS
############################################################

#### Determine if the alternative are SNPs
#### Count SNP's from all the sample individuals
def findsnps(A,sampledata):
	altlist = []
	poslist = []
	countlist = []
	count = 1
	snpexists = False
	#check if there are multiple alts
	if re.search(",",A):
		count = 1
		list = A.split(",")
		for i in list:
			if i != "*" and i != "." and len(i) == 1:
				altlist.append(i)
				poslist.append(count)
				snpexists = True
			count += 1

	else:
		if len(A) == 1 and A != "*" and A != ".":
				altlist.append(A)
				poslist.append(1)
				snpexists = True
	for y in poslist:
		sum = 0
		for x in sampledata:
			values = x.split(':')
			GT = values[0]
			if re.search(str(y),GT):
				sum += 1
		countlist.append(sum)
	return altlist,countlist,snpexists


#################################################
#Build Alt Codon and Send to Check Function
def checkaltcodon(ref,spotincodon,altinfo):
	#flag meaning: 1 = true / synonymous, 0 = false/nonsynonymous, and 2 = None/Neither
	altcodon = ''
	syncount = 0
	noncount = 0
	flag = 2
	### Look through list of alternative alleles
	### keep track of the list index for count referencing
	i = 0
	#print("spot:"+str(spotincodon))
	#print(altinfo)
	for items in altinfo[0]:
		#### Build the alternative codon for alt allele
		if spotincodon == 1:
			altcodon = items + refcodon[1:3]
		if spotincodon == 2:
			altcodon = ref[0]+items+ref[2]
		if spotincodon == 0:
			altcodon = ref[0:2]+items
		### Return if the two codons are synonymous or non-synonymous
		flag = translatedef(ref,altcodon)
		### Sum up the counts for either option with the counts for each alternative allele
		if flag == 1:
			syncount += a[1][i]
		if flag == 0:
			noncount += a[1][i]
		#print("alt:"+altcodon)
		#print("flag (1 = syn and 0 = non): " + str(flag))
		#print("syncount: "+ str(syncount))
		#print("noncount: "+ str(noncount))
		i += 1
	return syncount,noncount


#######################################################
#Check for Synonymous vs. Non-synonymous
def translatedef(ref,alt):
	flag = 2
	#check for valid codons
	if len(ref) != 3 or len(alt) != 3:
		print("WARNING: an invalid codon was passed")
	elif re.search("[^aAgGtTcC]",ref) or re.search("[^AagGtTcC]",alt):
		print("WARNING: an invalid codon was passed")

	coding_ref = Seq(ref)
	coding_alt = Seq(alt)
	ref_protein = coding_ref.translate()
	alt_protein = coding_alt.translate()
	if ref_protein != alt_protein:
		flag = 0
	else:
		flag = 1
	return flag



#### GOAL: Read population file. Build a dictionary of SNPs
#### SNP dictionary should contain: (Global / Genome Wide) Position, (all) SNP alleles, counts for each SNP, the chromosome/scaffold location

#Open population data file
vcfname = sys.argv[1]


################################################
#Initiate Empty Variables
###############################################
header_length = 0
snp = {}
snpexists = None
all_positions = []
positioninfo = ''
#################################################
#Check if population data file is compressed
###################################################
zippedflag = False
if re.search('.gz',vcfname):
	vcf_file = gzip.open(vcfname,'rb')
	zippedflag = True
else:
	vcf_file = open(vcfname,'r+') 

######################################################
###Iterate through lines of VCF File to find SNPs. 
#########################################################

for lines in vcf_file.readlines():

#If file was compressed, decode the lines into a readable form.
	if zippedflag == True:
		lines = lines.strip().decode(encoding='UTF-8')
#Find the HEADER Line. Goal: Retrieve the number of individuals sampled in the population.
	if re.search("CHROM",lines):
		headerlist = lines.split("\t")
		header_length = len(headerlist)
#Information will only be found on lines that do not start with '##' (#CHROM is the header line)
	elif lines[1] != '#':
		fields = lines.split("\t")
		chrom = fields[0]
		pos = fields[1]
		refallele = fields[3]
		altallele = fields[4]
		sampled_individuals = fields[9:header_length]
	### GOAL 1: Find which positions are SNPS. Keep track of alternative alleles in a list.
	### GOAL 2: Count the SNPS. If there is 1 SNP, count the 1's. If there are multiple SNPs (separated by commas), count their correspoding positions starting with 1 
	### Reference alleles longer than 1 allele represent deletions NOT SNPS.
	### Alternative alleles longer than 1 represent deletetions NOT SNPS.
	### Alternative alleles with marked as "*" are the same as reference allele. NOT SNPS.
		if len(refallele) == 1:
			list_of_alt_alleles = []
			counts_of_alt_alleles = []
			tuplecall = findsnps(altallele,sampled_individuals)
			list_of_alt_alleles = tuplecall[0]
			counts_of_alt_alleles = tuplecall[1]
			snpexists = tuplecall[2]
	### Store information about SNPs in dictionary for later access
	### For some VCF file, the chromosome name matters. Others, it does not.
			if snpexists:
				snp[(chrom,pos)] = (list_of_alt_alleles, counts_of_alt_alleles)

#print(snp)
###############################################################
#BUILDING REF CODONS
###############################################################
#Concatenate a reference codon for each position of interest from the given gene info map.
#My maps made positions in the follow format:  startpos:endpos
#If there were multiple starts and ends: startpos1:endpos1;startpos2:endpos2
#We also have reading frames for each section of starts and stops. Given either as 1,2, or 0
#OR a list seperated by a semi-colon such as 1;0;0;0 that corresponds to the start and stop positions.

########
#Open tab-sepated gene info list
#The script is written assuming that the file has the following fields:
#Orthogroup\tgeneid\ttx_id\tcds_id\tchomname\tpos\tframes
#Note: For this script, only geneid, pos, and frames are used for this script. The Orthogroup field will correspond to whatever gene information you want as the output.

ingene = sys.argv[2]
genes_info_list = open(ingene,'r+')

#########
#Open genomic file for your species. 
#This particular script is assuming that all the FASTA id's were removed so that the position in the string will be the position is the VCF file.
ingenome = sys.argv[3]
genome = open(ingenome,'r')
seq = {}
for record in SeqIO.parse(genome,"fasta"):
	seq[str(record.id)] = str(record.seq)
count = 1

################################################################################
##### GOAL: Output the counts of synonymous and non-synonymous polymorphisms.
###############################################################################
#Open output file
outname = sys.argv[4]
outfile = open(outname,'w+')
#Header line:
outfile.write("#Orthogroup\tGene-ID\tChromosome\tSynonymous-SNPs\tNon-synonymous-SNPs\n")

######################################
#build codons from positions of interest
#Build alternative codons from the reference codons for SNP positions
#sum up counts of synonymous-vs-non-synonymous SNPs for each genes
######################################
for lines in genes_info_list.readlines():
	if lines[0] != "#":

	#Keep track of non-synonymous and synonymous counts for each genes
		syn_runningsum = 0
		non_runningsum = 0
	#Store the fields from gene_info_list
		lines = lines.split("\t")
		ortho = lines[0]
		geneid = lines[3]
		chrom = lines[4]
		#For one of my species there was two chromsome indicators used in different files, I put both in the map with a "/"
		#Chrom2 will only be used if that "/" is found in the map field
		#For example, in my script, chom = BANY##### and chrom2 = NW.#########.#
		if re.search("/",chrom):
			mid = re.search("/",chrom).end()
			chrom2 = chrom[0:mid-1]
			chrom = chrom[mid:]
		positioninfo = lines[5]
		frameinfo = lines[6]
	#Store the proper reading frames
		if re.search(";",frameinfo):
			frame = frameinfo.split(";")
		else:
			frame = frameinfo.strip()
		if re.search(";", positioninfo):
			spot = ''
			refcodon = ''
			altcodon = ''
                        #If above is true, there are multiple start and stops positions
			mult_list = positioninfo.split(";")

                        #Iterate through each start/stop position
			it = 0

			for pos in mult_list:
				#Get the corresponding reading frame
				f = frame[it]
				f = f.strip()
				it += 1
				startandstop = pos.split(":")
				start = startandstop[0]
				end = startandstop[1]
				#iterate through the positions of interest
				local = 0
				for y in range(int(start),int(end)+1):
					local += 1
					try:
						a = snp[(chrom,str(y))]
						#print("frame:"+str(f)
						if int(f) == 0:
							spot = local % 3
						elif int(f) == 1:
							spot = (local+1) % 3
						elif int(f) == 2:
							spot = (local+2) % 3

						#### Build the codons
						#### The numbers are shifted by one because the seq is indexed w/ 0 but the positions are indexed with 1
						if spot == 1:
							if re.search("/",lines[4]):
								refcodon = seq[chrom2][y-1:y+2]
							else:
								refcodon = seq[chrom][y-1:y+2]
							#find alt codons here
						elif spot == 2:
							if re.search("/",lines[4]):
								refcodon = seq[chrom2][y-2:y+1]
							else:
								refcodon = seq[chrom][y-2:y+1]
						elif spot == 0:
							if re.search("/",lines[4]):
								refcodon = seq[chrom2][y-3:y]
							else:
								refcodon = seq[chrom][y-3:y]
						#print(y)
						#print("ref:"+refcodon)
						count_tuple = checkaltcodon(refcodon,spot,a)
						syn_runningsum += count_tuple[0]
						non_runningsum += count_tuple[1]
					except KeyError:
						pass

		else:
                #There is only one position for this CDS
		#There should also be only one reading frame
			refcodon = ''
			altcodon = ''
			spot = ''
			startandstop = positioninfo.split(":")
			f = frame
			start = startandstop[0]
			end = startandstop[1]
			local = 0
			for y in range(int(start),int(end)+1):
				local += 1
				try:
					a = snp[(chrom,str(y))]
					if int(f) == 0:
						spot = local % 3
					elif int(f) == 1:
						spot = (local+1) % 3
					elif int(f) == 2:
						spot = (local+2) % 3
					if spot == 1:
						if re.search("/",lines[4]):
							refcodon = seq[chrom2][y-1:y+2]
						else:
							refcodon = seq[chrom][y-1:y+2]
					elif spot == 2:
						if re.search("/",lines[4]):
							refcodon = seq[chrom2][y-2:y+1]
						else:
							refcodon = seq[chrom][y-2:y+1]
					elif spot == 0:
						if re.search("/",lines[4]):
							refcodon = seq[chrom2][y-3:y]
						else:
							refcodon = seq[chrom][y-3:y]
					#print(y)
					#print("ref:"+refcodon)
					count_tuple = checkaltcodon(refcodon,spot,a)
					syn_runningsum += count_tuple[0]
					non_runningsum += count_tuple[1]
				except KeyError:
					pass

		outfile.write(ortho+"\t"+geneid+"\t"+chrom+"\t"+str(syn_runningsum)+"\t"+str(non_runningsum)+"\n")





##### CLOSE FILES
genes_info_list.close()
outfile.close()
vcf_file.close()
