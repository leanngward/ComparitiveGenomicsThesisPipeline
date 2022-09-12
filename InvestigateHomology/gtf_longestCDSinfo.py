#transhttps://www.zybuluo.com/lakesea/note/1743213


###get_the_longest_transcripts.py
#!/usr/bin/env python3
# v0.1:fixing the problem that the feature of exon/mRNA is before gene
# v0.2: Use CDS to decide which transcript is primary
# v0.3: Fix the bug that if the CDS_parent include '\r' , the lst_len will be 0
# GFF must have CDS feature
# GFF must have ID and Parent in column 9
import re
import sys
import fileinput
from collections import defaultdict

infile = sys.argv[1]
outfile = sys.argv[2]

fileinput = open(infile,'r')

gene_dict = defaultdict(list)
tx_pos_dict = defaultdict(list)
CDS_dict = defaultdict(list)
CDS_location_dict = defaultdict(list)
CDS_frame_dict = defaultdict(list)
protein_dict = defaultdict(list)
for line in fileinput.readlines():
	if line.startswith("#"):
		continue
	content = line.split("\t")
	if len(content) <= 8:
		continue
	
	if content[2] == "transcript" or content[2] == "mRNA":
		
		# use regular expression to extract the gene ID
	        # match the pattern ID=xxxxx; or ID=xxxxx
		tx_id = re.search(r'transcript_id (.*?);',content[8]).group(1)
		tx_parent =  re.search(r'gene_id (.*?)[;]',content[8]).group(1)
		tx_parent = tx_parent.strip() # remove the 'r' and '\n'

		# if the parent of transcript is not in the gene_dict, create it rather than append
		if tx_parent in gene_dict:
			gene_dict[tx_parent].append(tx_id)
		else:
			gene_dict[tx_parent] = [tx_id]
		tx_pos_dict[tx_id] = [ content[0], content[3], content[4], content[6]]
	if content[2] == 'CDS':
		width = int(content[4]) - int(content[3])
		CDS_parent = re.search(r'transcript_id (.*?)[;]',content[8]).group(1)
		CDS_parent = CDS_parent.strip() # strip the '\r' and '\n'
		CDS_dict[CDS_parent].append(width)
		startandstop = content[3]+":"+content[4]
		CDS_location_dict[CDS_parent].append(startandstop)
		CDS_frame_dict[CDS_parent].append(content[7])

		protein = re.search(r'protein_id (.*?)[;]',content[8]).group(1)
		protein = protein.strip('"')
		protein_dict[CDS_parent].append(protein)

fileoutput = open(outfile,'w')
fileoutput.write("#geneid\ttxid\tchrompos\tpos\tframes\n")

for gene, txs in gene_dict.items():
	tmp=0
	location = ''
	frameloc = ''
#	print("Gene:",gene)
#	print("IDs:",txs)

	for tx in txs:
#		print("item:",tx)
		tx_len = sum(CDS_dict[tx])
#		print("length:",tx_len)
		if tx_len == 0:
			lst_tx = tx
		elif tx_len > tmp:
			lst_tx = tx
			tmp = tx_len
	loc = CDS_location_dict[lst_tx]
	frame = CDS_frame_dict[lst_tx]
	prot = protein_dict[lst_tx]
#	print(lst_tx)
#	print(":")
#	print(prot)
	for l in loc:
		location += l
		location += ";"
	location = location[0:-1]

	for i in frame:
		frameloc += i
		frameloc += ";"
	frameloc = frameloc[0:-1]

	tx_chrom = tx_pos_dict[lst_tx][0]
	tx_start = tx_pos_dict[lst_tx][1]
	tx_end   = tx_pos_dict[lst_tx][2]
	tx_strand = tx_pos_dict[lst_tx][3]

	gene = gene.strip('"')
	tx = tx.strip('"')	
	fileoutput.write(tx+"\t"+prot[0]+"\t"+tx_chrom+"\t"+location+"\t"+frameloc+"\n")
#	print("{gene}\t{tx}\t{chrom}\t{start}\t{end}\t{strand}".format(gene=gene,tx=lst_tx,chrom=tx_chrom,start=tx_start,end=tx_end, strand=tx_strand))

fileinput.close()
fileoutput.close()
