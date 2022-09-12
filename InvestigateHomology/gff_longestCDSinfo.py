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


#NOTE: A GFF file has the following column fields: seqname \t \t source \t feature \t start \t end \t score \t strand +/- \t frame \t attributes

#the gene_dict has its key as the parent ID and its values as the as the list of transcript IDs for each mRNA
gene_dict = defaultdict(list)
#the tx_pos_dict has the ID as the key and the seqname, start, end, and strand information for each mRNA
tx_pos_dict = defaultdict(list)

#the CDS_dict has each parent ID as its keys and a list of the widths as the value
#the parent ID will match its respective mRNA
CDS_dict = defaultdict(list)

CDS_location_dict = defaultdict(list)
CDS_frame_dict = defaultdict(list)
CDS_id_dict = defaultdict(list)
#read in the file of interest
infile = sys.argv[1]
fileinput = open(infile,'r')

outfile = sys.argv[2]
fileoutput = open(outfile, 'w')

for line in fileinput.readlines():
	if line.startswith("#"):
		continue
	content = line.split("\t")
	if len(content) <= 8:
		continue
	
	if content[2] == "transcript" or content[2] == "mRNA":
		
		# use regular expression to extract the gene ID
	        # match the pattern ID=xxxxx; or ID=xxxxx

		#find the first ID instance
		tx_id = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
		#find the first parent instance
		tx_parent =  re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
		tx_parent = tx_parent.strip() # remove the 'r' and '\n'

		# if the parent of transcript is not in the gene_dict, create it rather than append
		if tx_parent in gene_dict:
			gene_dict[tx_parent].append(tx_id)
		else:
			gene_dict[tx_parent] = [tx_id]
		tx_pos_dict[tx_id] = [ content[0], content[3], content[4], content[6], content[7] ]
	if content[2] == 'CDS':
		width = int(content[4]) - int(content[3])
		CDS_parent = re.search(r'Parent=(.*?)[;\n]',content[8]).group(1)
		CDS_parent = CDS_parent.strip() # strip the '\r' and '\n'
		CDS_dict[CDS_parent].append(width)
		startandstop = content[3]+":"+content[4]
		CDS_location_dict[CDS_parent].append(startandstop)
		CDS_frame_dict[CDS_parent].append(content[7])
		cdsid = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
		CDS_id_dict[CDS_parent] = re.search(r'ID=(.*?)[;\n]',content[8]).group(1)
#loops through each parent, ID pair
#the ID's are the mRNA/transcript lines
print(len(gene_dict))
fileoutput.write("#geneid\ttx_id\tcds_id\tchromname\tpos\tframes\n")

for gene, txs in gene_dict.items():
	tmp=0
	location = ''
	frameloc = ''
	#For each transcript ID with the same parent
	for tx in txs:
		#sum the CDS width for a single transcript ID
		tx_len = sum(CDS_dict[tx])
		#see if this transcript ID is the longest so far
		if tx_len > tmp:
			lst_tx = tx
			tmp = tx_len
#	print("longest:",lst_tx)

	tx_chrom = tx_pos_dict[lst_tx][0]
#	tx_start = tx_pos_dict[lst_tx][1]
#	tx_end   = tx_pos_dict[lst_tx][2]
	startandstops = CDS_location_dict[lst_tx]
	frame = CDS_frame_dict[lst_tx]
	for i in frame:
		frameloc += i
		frameloc += ";"
	frameloc = frameloc[0:-1]

	for items in startandstops:
		location += items
		location += ';'
	location = location[0:-1]
#	print(location)
#	print("frame:\n")
#	print(frameloc)
	tx_strand = tx_pos_dict[lst_tx][3]
	tx_frame = tx_pos_dict[lst_tx][4]

	tx_cdsid = CDS_id_dict[lst_tx]
#	print(tx_cdsid)

	fileoutput.write(gene+"\t"+lst_tx+"\t"+tx_cdsid+"\t"+tx_chrom+"\t"+location+"\t"+frameloc+"\n")

#	print("{gene}\t{tx}\t{chrom}\t{start}\t{end}\t{strand}".format(gene=gene,tx=lst_tx,chrom=tx_chrom,start=tx_start,end=tx_end, strand=tx_strand))


fileinput.close()
fileoutput.close()
