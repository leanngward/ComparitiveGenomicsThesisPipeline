
## Set Up Environment and Variables ####
  #Set working directory
  setwd("C:/Users/leann/OneDrive/grad research_newer/butterfly_aging/run5_asproteins_orthofinder/")

#Set evodata file name
  infile <- "longevitygenesonly_results/branch_output_heliconius_clade.tsv"
  
  #Set go file name
  goterms <- "longevitygenesonly_results/run5_sc_output_orfs_Pfam-A_GO.tsv"
  
  #Set mapping file name

  #The mapfile is the go term mapping with ONLY heliconius (you can choose any species) and the the speciesname_ removed
  #Edit in excel first
  mapfile <- "longevitygenesonly_results/unique_go_term_mapping_sc_hmel.txt"
  list <- "longevitygenesonly_results/genage_go_terms.txt"
  
  ## Find the subset of orthogroups in the longevity list - IF relevant####
  
  listdf <- read.delim(list,header=TRUE,sep="\t")
  longevitylist <- listdf$V5
  head(longevitylist)
  
  mapdf <- read.delim(mapfile,header=FALSE,sep="\t")
  
  inlist <- rep("NA",1)
  for (row in 1:nrow(mapdf)){
    go <- unlist(strsplit(mapdf[row,2],","))
    for (i in go) {
      if (i %in% longevitylist) {
        inlist <- append(inlist,mapdf[row,1])
      }
    }
    
  }
  
  head(inlist)
  inlist <- na.omit(inlist)
  inlist <- unique(inlist)
  print(inlist)
  
## Load in the Evolution Model Output ####
  evodata <- read.delim(infile, header = TRUE, sep = "\t")
  nrow(evodata)
  head(evodata)
## grab only the longevity subset, if applicable ####

  evodata <- evodata[evodata$Group %in% inlist, ]
  
## Find significant omega values ####
  #Calculate log-likelihood ratio
  evodata$LR <- 2 * (evodata$Alt.lnL - evodata$Null.lnL)

  #Calculate p-value
  evodata$pval <- pchisq(evodata$LR, df = 1, lower.tail = FALSE)
  #adjust p-value to control false discovery rate
  evodata$padj <- p.adjust(evodata$pval, method = 'fdr')
  
  #create table with only significant results
  sigevodata <- evodata[evodata$padj < .05, ]
  nrow(sigevodata)  
  
  head(sigevodata)
## Look at omega distribution for alternative model ####
  # Libraries
  library(ggplot2)
  library(dplyr)
  
   
  # Make the histogram
  sigevodata %>%
    #filter( price<300 ) %>%
    ggplot( aes(x=Alt.Omega.1)) +
    geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8)
  
  head(sigevodata)
  n <- 90
  topset <- sigevodata[sigevodata$Alt.Omega.1 > quantile(sigevodata$Alt.Omega.1,prob=1-n/100,na.rm=TRUE),]
  testvalue <- min(topset$Alt.Omega.1)
  print(testvalue)
  #find a percentage of the top and bottom values
  
## Load TOPGO ####
  remove.packages("topGO")
  install.packages('rlang')
  install.packages("devtools")
  devtools::install_github("ycl6/topGO-feat", ref = "v2.41.0-barplot")
  library(topGO)
  
## Create appropriate variables and run topGO ####
  
  #Use the topGO geneSelectionFun when creating the topGO object? 
  
  #Custom annotations - for each gene identifier are listed GO terms to which this gene is specifically annotated
  
    #Create custom annotation with topGO function readMappings
    #must open the file directly. doesn't read from a table.
    geneID2GO <- readMappings(file = mapfile, sep = "\t", IDsep = ",")
    GO2geneID <- inverseList(geneID2GO)
    str(head(geneID2GO))
    length(geneID2GO)
    #named vector of the numerical values where the names are the gene indentifiers. makes the gene universe
    
    subset <- sigevodata
    #Here I am using the foreground omega values
    geneuniverse <- as.numeric(subset$Alt.Omega.1)
    names(geneuniverse) <- subset$Group
    #geneuniverse <- factor(as.integer(evodata$padj < .05))
    #names(geneuniverse) <- evodata$Group
    head(geneuniverse)
    length(geneuniverse)

    #Write a function to choose select genes based on their scores (topGO uses adjusted p-values)
    topDiffGenes <- function(allScore) {
      return(allScore < testvalue )
    }
    geneList <- topDiffGenes(geneuniverse)
    length(geneList)
    topobject_BP <- new("topGOdata", ontology = 'BP', allGenes = geneuniverse, geneSelectionFun = topDiffGenes,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    
    topobject_MF <- new("topGOdata", ontology = 'MF', allGenes = geneuniverse, geneSelectionFun = topDiffGenes,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    
    
    topobject_CC <- new("topGOdata", ontology = 'CC', allGenes = geneuniverse, geneSelectionFun = topDiffGenes,
                     annot = annFUN.gene2GO, gene2GO = geneID2GO)
    
## Explore topGO results w/ statistical tests ####
    resultFisher_BP <- runTest(topobject_BP, algorithm="classic", statistic="fisher")
    #resultFisher   
    
    resultFisher_MF <- runTest(topobject_MF, algorithm="classic", statistic="fisher")
    resultFisher_CC <- runTest(topobject_CC, algorith="classic", statistic="fisher")
    
    #resultKS <- runTest(topobject, algorithm="classic", statistic="ks")
    resultKS.elim.BP <- runTest(topobject_BP, algorithm="elim", statistic="ks")
    resultKS.elim.MF <- runTest(topobject_MF, algorithm="elim", statistic="ks")
    resultKS.elim.CC <- runTest(topobject_CC, algorithm="elim", statistic="ks")
    #resultKS    

    
    allRes_BP <- GenTable(topobject_BP, classicFisher = resultFisher_BP, elimKS = resultKS.elim.BP, 
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 50)  
    
    allRes_MF <- GenTable(topobject_MF, classicFisher = resultFisher_MF, elimKS = resultKS.elim.MF, 
                          orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 50)
    
    allRes_CC <- GenTable(topobject_CC, classicFisher = resultFisher_CC, elimKS = resultKS.elim.CC, 
                       orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 50)    
    #view(allRes)

    
    ## Explore GO network
    #BiocManager::install("Rgraphviz")
#    library(Rgraphviz)
    showSigOfNodes(topobject_BP,score(resultKS.elim.BP), firstSigNodes = 5, useInfo = 'all')
    showSigOfNodes(topobject_MF,score(resultKS.elim.MF), firstSigNodes = 5, useInfo = 'all')
    showSigOfNodes(topobject_CC,score(resultKS.elim.CC), firstSigNodes = 5, useInfo = 'all')
    
## Explore results with a Bar Plot ####


   enrichment_barplot(topobject_BP, resultKS.elim.BP, showTerms = 20, numChar = 40, 
                       orderBy = "Scores", y = "Count", 
                       xlab = NULL, ylab = NULL, title = NULL)
   enrichment_barplot(topobject_MF, resultKS.elim.MF, showTerms = 20 , numChar = 40, 
                      orderBy = "Scores", y = "Count", 
                       xlab = NULL, ylab = NULL, title = NULL)
  
   enrichment_barplot(topobject_CC, resultKS.elim.CC, showTerms = 20, numChar = 40, 
                       orderBy = "Scores", y = "Count", 
                       xlab = NULL, ylab = NULL, title = NULL)
    