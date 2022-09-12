
## Load Libraries ####

  #library(GO.db)
  #library(biomaRt)
  #library(Rgraphviz)
  #library(affy)
  #library(ALL)
  #library(data.table)
  #library(rlang)
  #library(DOSE)
  #install.packages('rlang')
  #remove.packages('rlang')
  #data(ALL)
  #affyLib <- paste(annotation(ALL), "db", sep=".")
  #library(package = affyLib, character.only = TRUE)
  #devtools::install_github("YuLab-SMU/DOSE")
  #devtools::install_github("YuLab-SMU/clusterProfiler")

  library(topGO)
  library(clusterProfiler)


## Set Working Directory ###

setwd("C:/Users/leann/OneDrive/grad research_newer/butterfly_aging/run5_asproteins_orthofinder/gomapping_singlecopy_files/Duplication Analysis")
  

## Give In File Names ####

#Annotations INFO
#The mapping of all genes are contained within the all the genomes that were run to get the gene duplication results.
#Use a _orfs_Pfam-A_GO_GOterm_mapping.tsv file from HMMER2GO
#It has the following tab-delimited columns: GENEID \t semi-colon delimited LIST OF GOIDS (GO:######)
#Structure: geneID2GO <- readMappings(file = "NAME_orfs_Pfam-A_GO_GOterm_mapping.tsv")
# This type of file is a two column table with the sequence name in the first column and the GO terms associated with that sequence. 
#geneID2GO serves as a custom annotation file for all of your genes and their possible GO terms

  fextension <- "orfs_Pfam-A_GO_GOterm_mapping.tsv"
  fname <- "run5_allspec_background"
  inFile <- paste(fname,fextension,sep="_")
  
  inFile2 <- paste(fname,"orfs_Pfam-A_GO.tsv",sep="_")

  #All the groups of interest will be visualized at once. Provide the files names here.
  #It will usually have groups. Add or take away lines to the code as needed. 
  

  group1 <- "B_N6_N3"
  group2 <-  "D_N5_N3"
  group3 <- "D_N7_N5"
  group4 <- "H_N7_N5"
  group5 <- "H_N8_N7"
  group6 <- "PA_N6_N3"
  
  
## Load in Data ###
  #Annotation Files
  
  geneID2GO <- readMappings(file = inFile)
  str(head(geneID2GO))

  #Assign a list of GENEID's with a name from the readMapping
  geneNames <- names(geneID2GO)
  geneNames
  
  #Groups of Interest Gene ID's
  group1_data <- read.delim(group1, header=FALSE)
  group2_data <- read.delim(group2, header=FALSE)
  group3_data <- read.delim(group3, header=FALSE)
  group4_data <- read.delim(group4, header=FALSE)
  group5_data <- read.delim(group5, header=FALSE)
  group6_data <- read.delim(group6, header=FALSE)

  
  #all Genes Data
  #Reload background (all genes) as a delimited file
  BG_data <- read.delim(inFile, header = FALSE)
  head(BG_data)

## PREP DATA ####
#Test which GENEIDS from your background that are in the group data (or "candidate" genes ids) using %in%
#Process the TRUE or FALSE as an integer with value 0 = FALSE and 1 = TRUE
#Set each GENEID as a factor with levels: 0 and 1 
#where 1 indicates that the GENEID is from the GROUPOFINTEREST  
  
  
  cands1 <- factor(as.integer(BG_data$V1 %in% group1_data$V1))
  cands2 <- factor(as.integer(BG_data$V1 %in% group2_data$V1))
  cands3 <- factor(as.integer(BG_data$V1 %in% group3_data$V1))
  cands4 <- factor(as.integer(BG_data$V1 %in% group4_data$V1))
  cands5 <- factor(as.integer(BG_data$V1 %in% group5_data$V1))
  cands6 <- factor(as.integer(BG_data$V1 %in% group6_data$V1))

  
  #View summary of the factors if interested in counts of 0 vs. 1
  summary(cands2)  
  
  #Name the factors (0 and 1) using the GENEID's from BACKGROUND - don't worry there's the same amount and in the same order! :)
  names(cands1) <- BG_data$V1 
  names(cands2) <- BG_data$V1
  names(cands3) <- BG_data$V1
  names(cands4) <- BG_data$V1
  names(cands5) <- BG_data$V1
  names(cands6) <- BG_data$V1

  
  #Check list by viewing as a string
  str(cands1)    

  
## USE TOPGO FEATURES ####
  #Measure enrichment of GO terms
  #gene2GO contains the GO terms for each GENEID
  #allGenes contains the file of candidate genes (OR, the genes we are interested in). They are candidates if their level indicates 1 = TRUE
  #annot indicates that it should use gene2GO object as the mapping file.

  #Group 1 - testing for ontology CC
  GOgroup_1_CC <- new("topGOdata", ontology = 'CC', allGenes = cands1,
                       annot = annFUN.gene2GO, gene2GO = geneID2GO)


  #Group 2 - testing for ontology CC
  GOgroup_2_CC <- new("topGOdata", ontology = 'CC', allGenes = cands2,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO)


  #Group 3 - testing for ontology CC
  GOgroup_3_CC <- new("topGOdata", ontology = 'CC', allGenes = cands3,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO)


  #Group 4 - testing for ontology CC, CC
  GOgroup_4_CC <- new("topGOdata", ontology = 'CC', allGenes = cands4,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)


  #Group 5 - testing for ontology CC
  GOgroup_5_CC <- new("topGOdata", ontology = 'CC', allGenes = cands5,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)  
  
  #Group 6 - testing for ontology CC
  GOgroup_6_CC <- new("topGOdata", ontology = 'CC', allGenes = cands6,
                      annot = annFUN.gene2GO, gene2GO = geneID2GO)  
  
  
## COMPUTE STASTICS FOR TOPGO OBJECT  ####
#We will need to run statical test on our TOPGO object to calculate the p-values.
#Fisher's Exact test is based on gene counts.However, since our topgo object contains gene scores representing how differentially the gene is expressed
#the Kolmogorov-Smirnov test is also an option.

#The statistical test command is runTest()
#Requires: (topgoobject, an algorithm, and test type)
#'weight01' is topgo's default algorithm

  
  f_res_1_CC = runTest(GOgroup_1_CC, algorithm='weight01', statistic='fisher')
  f_res_2_CC = runTest(GOgroup_2_CC, algorithm='weight01', statistic='fisher')
  f_res_3_CC = runTest(GOgroup_3_CC, algorithm='weight01', statistic='fisher')
  f_res_4_CC = runTest(GOgroup_4_CC, algorithm='weight01', statistic='fisher')
  f_res_5_CC = runTest(GOgroup_5_CC, algorithm='weight01', statistic='fisher')
  f_res_6_CC = runTest(GOgroup_6_CC, algorithm='weight01', statistic='fisher')
  
  
  #Define terms that are available for analysis using usedGO() so that we can access the length of our relative result list

  allGO_1_CC = usedGO(GOgroup_1_CC)
  allGO_2_CC = usedGO(GOgroup_2_CC)
  allGO_3_CC = usedGO(GOgroup_3_CC)
  allGO_4_CC = usedGO(GOgroup_4_CC)
  allGO_5_CC = usedGO(GOgroup_5_CC)
  allGO_6_CC = usedGO(GOgroup_6_CC)
  
  
## PREP DATA FOR VISUALIZATION ####
#TopGO provides the GenTable feature to interpret the results.
#This function is for analyzing the most significant GO terms and their corresponding p-values.
#This is the data we need to calculate upregulation and downregulation of our significant terms.

  ###### GENE ONTOLOGY TERMS ########
  res_1_CC = GenTable(GOgroup_1_CC, f = f_res_1_CC, orderBy='f',topNodes=length(allGO_1_CC))
  res_2_CC = GenTable(GOgroup_2_CC, f = f_res_2_CC, orderBy='f',topNodes=length(allGO_2_CC))
  res_3_CC = GenTable(GOgroup_3_CC, f = f_res_3_CC, orderBy='f',topNodes=length(allGO_3_CC))
  res_4_CC = GenTable(GOgroup_4_CC, f = f_res_4_CC, orderBy='f',topNodes=length(allGO_4_CC))
  res_5_CC = GenTable(GOgroup_5_CC, f = f_res_5_CC, orderBy='f',topNodes=length(allGO_5_CC))
  res_6_CC = GenTable(GOgroup_6_CC, f = f_res_6_CC, orderBy='f',topNodes=length(allGO_6_CC))

##### Calculate fold-change to find upregulation and downregulation ########
#Start with the significant results from res_ objects calculated above.

  sig_1_CC <- res_1_CC[res_1_CC$f < .05, ]
  sig_2_CC <- res_2_CC[res_2_CC$f < .05, ]
  sig_3_CC <- res_3_CC[res_3_CC$f < .05, ]
  sig_4_CC <- res_4_CC[res_4_CC$f < .05, ]
  sig_5_CC <- res_5_CC[res_5_CC$f < .05, ]
  sig_6_CC <- res_6_CC[res_6_CC$f < .05, ]


  
  sig_1_CC$FC <- ifelse(sig_1_CC$Expected > sig_1_CC$Significant, -1 * (sig_1_CC$Expected / sig_1_CC$Significant) ,sig_1_CC$Significant / sig_1_CC$Expected)

  sig_2_CC$FC <- ifelse(sig_2_CC$Expected > sig_2_CC$Significant, -1 * (sig_2_CC$Expected / sig_2_CC$Significant) ,sig_2_CC$Significant / sig_2_CC$Expected)     

  sig_3_CC$FC <- ifelse(sig_3_CC$Expected > sig_3_CC$Significant, -1 * (sig_3_CC$Expected / sig_3_CC$Significant) ,sig_3_CC$Significant / sig_3_CC$Expected)

  sig_4_CC$FC <- ifelse(sig_4_CC$Expected > sig_4_CC$Significant, -1 * (sig_4_CC$Expected / sig_4_CC$Significant) ,sig_4_CC$Significant / sig_4_CC$Expected)

  sig_5_CC$FC <- ifelse(sig_5_CC$Expected > sig_5_CC$Significant, -1 * (sig_5_CC$Expected / sig_5_CC$Significant) ,sig_5_CC$Significant / sig_5_CC$Expected)
  
  sig_6_CC$FC <- ifelse(sig_6_CC$Expected > sig_6_CC$Significant, -1 * (sig_6_CC$Expected / sig_6_CC$Significant) ,sig_6_CC$Significant / sig_6_CC$Expected)
  


  #Label each set's group name 
  sig_1_CC$group <- group1
  sig_2_CC$group <- group2
  sig_3_CC$group <- group3
  sig_4_CC$group <- group4
  sig_5_CC$group <- group5
  sig_6_CC$group <- group6

  #Combine groups you want to compare for the visualization
  mydf <- rbind(sig_1_CC, sig_2_CC, sig_3_CC,  sig_4_CC, sig_5_CC, sig_6_CC)   
  #Keep columns for ID, FC, and GROUP
  mydf <- mydf[ , c(1, 6, 7, 8)]
  
  
## CLUSTER PROFILER - VISUALIZE RESULTS ####
#Create a gene universe with all GOterms
  
  #term2go should be a data.frame with first column of term ID and second column of corresponding mapped gene
  terms <- read.delim(inFile2, header = FALSE)
  term2go <- data.frame(term=terms$V6, goterm=terms$V5)
  
  mydf_ordered <- mydf[order(-mydf$FC), ]
  groupNames_ordered <- mydf_ordered$group

  
  #term2names should be a data.frame with the first column of ID and the second column of corresponding term name. optional!
  filename <- "run5_allspec_background_orfs_Pfam-A_GO.tsv"
  all_group_names <- read.delim(filename,header=FALSE)
  term2names <- all_group_names[ , c(6,1)]
  
  #Transform data
  #Grab only the FC values
  geneList_clus = mydf_ordered[ , 3]
  print(geneList_clus)
  
  #Name each FC with its respective GOID
  names(geneList_clus) = as.character(mydf_ordered[,1])
  names(geneList_clus)
  
  #Sort list by decreasing FC order
  #geneList_clus = sort(geneList_clus, decreasing = TRUE)

  #Created the data the visualization will come from with the relavant modhead(el aspects as the column names.
  #Create variable "Entrez" for compareCluster model. Entrez is the FC's named with their GO terms that I am testing.
  #Include the respective "othergroup" which will function as one of the treatments in the model.
  mydf_clus = data.frame(Entrez=names(geneList_clus), FC = geneList_clus, othergroup = groupNames_ordered)
  #Define each FC as "upregulated" or "downregulated". Upregulated values are positive. Downregulated values are negative.
  #These values are defined by "group" column name and is one of the treatments for our cluster model.
  mydf_clus$group <- "upregulated"
  mydf_clus$group[mydf_clus$FC < 0] <- "downregulated"


  #Apply the compareCluster feature from clusterProfiler features
  #rm(formula_res_clus)

  formula_res_clus <- clusterProfiler::compareCluster(Entrez ~ othergroup+group, data=mydf_clus, 
                                                    fun= "enricher", TERM2GENE = term2go, TERM2NAME = NA, pAdjustMethod = "none",
                                                    pvalueCutoff = 0.05, 
                                                    universe = NA,
                                                    minGSSize = 1,
                                                    maxGSSize = 1000,
                                                    qvalueCutoff = 1)
  head(formula_res_clus)
  dotplot(formula_res_clus, x="group",showCategory = 5) + DOSE::facet_grid(~othergroup)
?dotplot
?facet_grid


  