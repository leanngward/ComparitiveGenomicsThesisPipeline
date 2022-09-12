## Load Libraries ####


## Set Working Directory ####
  setwd("C:/Users/leann/OneDrive/grad research_newer/butterfly_aging/run5_asproteins_orthofinder/branch_site_results/")

## Name input files ####
  inFile = "branchsite_output_pier_clade.tsv"
  
  BEB = "BEB_output_pier_clade.txt"
  BEBcounts = "BEB_counts_pier_clade.txt"
  
  mapfile = "unique_go_term_mapping_sc_hmel.txt"
## Load in Branchsite Output Files ####
  evodata <- read.delim(inFile, header = TRUE, sep = "\t")
  BEBdata <- read.delim(BEBcounts,header=TRUE,sep = "\t")
  head(BEBdata)  
  nrow(BEBdata)
## Find significant omega values ####
  #Calculate log-likelihood ratio
  evodata$LR <- 2 * (evodata$Alt.lnL - evodata$Null.lnL)

  #Calculate p-value
  evodata$pval <- pchisq(evodata$LR, df = 1, lower.tail = FALSE)
  #adjust p-value to control false discovery rate
  evodata$padj <- p.adjust(evodata$pval, method = 'fdr')
  head(evodata)

  #merge table with with BEB counts by Groupname 
  fulldf <- merge(evodata, BEBdata, by.x='Group', by.y='GroupName')
  head(fulldf)

  #create table with only significant results
  sigevodata <- fulldf[fulldf$padj < .05, ]
  #sigevodata <- fulldf
  nrow(sigevodata)  
  head(sigevodata)

  sum(sigevodata$Total.Sig.Sites)  

## Look at the distribution of site numbers to help choose a cut-off ####
  library(ggplot2)
  library(dplyr)
  
  # Make the histogram
  sigevodata %>%
    filter( Total.Sig.Sites< 20000 ) %>%
    ggplot( aes(x=Total.Sig.Sites)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)
  
  n <- 10
  topset <- sigevodata[sigevodata$Total.Sig.Sites > quantile(sigevodata$Total.Sig.Sites,prob=1-n/100,na.rm=TRUE),]
  testvalue <- min(topset$Total.Sig.Sites)
  print(testvalue)
  
#  filter( price<300 ) %>%
## Load TOPGO ####
  #remove.packages("topGO")
  #install.packages('rlang')
  #install.packages("devtools")
  #devtools::install_github("ycl6/topGO-feat", ref = "v2.41.0-barplot")
  library(topGO)

## Create appropriate variables and run topGO ####
  
  #Use the topGO geneSelectionFun when creating the topGO object? 
  
  #Custom annotations - for each gene identifier are listed GO terms to which this gene is specifically annotated
  
    #Create custom annotation with topGO function readMappings
    #must open the file directly. doesn't read from a table.
    geneID2GO <- readMappings(file = mapfile, sep = "\t", IDsep = ",")
    GO2geneID <- inverseList(geneID2GO)
    str(head(geneID2GO))
    length(GO2geneID)
    
    #named vector of the numerical values where the names are the gene indentifiers. makes the gene universe
    subset <- sigevodata
    
    #Here I am using the the number of selected sites
    geneuniverse <- as.numeric(subset$Total.Sig.Sites)
    names(geneuniverse) <- subset$Group
    
    #geneuniverse <- factor(as.integer(evodata$padj < .05))
    #names(geneuniverse) <- evodata$Group
    head(geneuniverse)
    
    #Write a function to choose select genes based on their scores (topGO uses adjusted p-values)
    topDiffGenes <- function(allScore) {
      return(allScore > 0 )
    }
    geneList <- topDiffGenes(geneuniverse)
    topobject_BP <- new("topGOdata", ontology = 'BP', allGenes = geneuniverse, geneSelectionFun = topDiffGenes,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    
    topobject_MF <- new("topGOdata", ontology = 'MF', allGenes = geneuniverse, geneSelectionFun = topDiffGenes,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    
    
    topobject_CC <- new("topGOdata", ontology = 'CC', allGenes = geneuniverse, geneSelectionFun = topDiffGenes,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)  

## Explore topGO results w/ statistical tests ####
    ?runTest
    resultFisher_BP <- runTest(topobject_BP, algorithm="classic",statistic="fisher") 
    
    resultFisher_MF <- runTest(topobject_MF, algorithm="classic",statistic="fisher")
    
    resultFisher_CC <- runTest(topobject_CC, algorithm="classic",statistic="fisher")
    
    
    resultKS.elim.BP <- runTest(topobject_BP, algorithm="elim", statistic="ks")
    resultKS.elim.MF <- runTest(topobject_MF, algorithm="elim", statistic="ks")
    resultKS.elim.CC <- runTest(topobject_CC, algorithm="elim", statistic="ks")
    
    resultKS.BP <- runTest(topobject_BP, algorithm = "classic", statistic = "ks")
    resultKS.MF <- runTest(topobject_MF, algorithm = "classic", statistic = "ks")
    resultKS.CC <- runTest(topobject_CC, algorithm = "classic", statistic = "ks")
    
    allRes_BP <- GenTable(topobject_BP, classicFisher = resultFisher_BP, elimKS = resultKS.elim.BP, classicKS = resultKS.BP,
                            orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)
    
    
    allRes_MF <- GenTable(topobject_MF, classicFisher = resultFisher_MF, elimKS = resultKS.elim.MF,classicKS = resultKS.MF,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)

    allRes_CC <- GenTable(topobject_CC, classicFisher = resultFisher_CC, elimKS = resultKS.elim.CC,classicKS = resultKS.CC,
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 100)    
## Explore GO Network ####
    showSigOfNodes(topobject_BP,score(resultKS.elim.BP), firstSigNodes = 10, useInfo="all")
    showSigOfNodes(topobject_MF,score(resultKS.elim.MF), firstSigNodes = 10, useInfo="all")
    showSigOfNodes(topobject_CC,score(resultKS.elim.CC), firstSigNodes = 10, useInfo="all")
## Explore Results with Bar Plot
    
    
    # KS Elim Method
    enrichment_barplot(topobject_BP, resultKS.elim.BP, showTerms = 20, numChar = 40,
                       orderBy = "Scores", y = "Count", xlab = NULL, ylab = NULL, title = NULL)
    
    enrichment_barplot(topobject_MF, resultKS.elim.MF, showTerms = 20, numChar = 40,
                       orderBy = "Scores", y = "Count", xlab = NULL, ylab = NULL, title = NULL)
    enrichment_barplot(topobject_CC, resultKS.elim.CC, showTerms = 20, numChar = 40,
                       orderBy = "Scores", y = "Count", xlab = NULL, ylab = NULL, title = NULL)
    
    # Fisher Method
    #enrichment_barplot(topobject_BP, resultFisher_BP, showTerms = 10, numChar = 40,
                       orderBy = "Scores", y = "Count", xlab = NULL, ylab = NULL, title = NULL)
    
    #enrichment_barplot(topobject_MF, resultFisher_MF, showTerms = 10, numChar = 40,
                       orderBy = "Scores", y = "Count", xlab = NULL, ylab = NULL, title = NULL)
    #enrichment_barplot(topobject_CC, resultFisher_CC, showTerms = 10, numChar = 40,
                       orderBy = "Scores", y = "Count", xlab = NULL, ylab = NULL, title = NULL)
    