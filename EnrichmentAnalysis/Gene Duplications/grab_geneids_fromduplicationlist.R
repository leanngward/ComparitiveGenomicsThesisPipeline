
## Load Libraries ####
  library(tidyverse)

## Set Working Directory ####
  setwd("C:/Users/leann/OneDrive/grad research_newer/butterfly_aging/run5_asproteins_orthofinder/gomapping_singlecopy_files")

## Set File Names ####
  
  inFile = "Duplications.tsv"
  
  
## Load in duplication data from Orthofinder ####
  
  alldups <- read.delim(inFile, header = TRUE)
  head(alldups)
  nrow(alldups)
  
  
## Set the nodes of interest ####
  node1 = "N6"
  node2 = "N3"
  
  
  
## Set species of interest ####
  
  species = "PA"
  
## Make Outfile Name ####
  
  outFile = paste(species,node1,node2,sep="_")
  outFile  
## Get the gene ID's for interest ####
  
  #Lets get a group of duplications for our species of interest. 
  #Also, we're only going to get keep nodes with at least 50% support
  
  
  #Grab only this group
    dups <- alldups[alldups$"Species.Tree.Node" == node1 | alldups$"Species.Tree.Node" == node2, ]
  #Keep groups with at least 50% support
    dups <- dups[dups$"Support" >= 0.5, ]
    dups[, 2]
    nrow(dups)
  #Get a list of Gene.1 and Gene.2 columns
    genes1 <- unlist(strsplit(as.character(dups[1,6]),","))
    genes2 <- unlist(strsplit(as.character(dups[1,7]),","))
  #Get the first term in Gene.1 and Gene.2
    genes1 <- trimws(genes1, which = c("both"))
    genes2 <- trimws(genes2, which = c("both"))
  #initiate list of gene ids for this group
    genes <- genes1
    genes <- append(genes, genes2)

  #This section gets the rest of Gene.1 and Gene.2 for this group
  #Loop through every row in within the group
    for (row in 2:nrow(dups)) {
      #splits genes 1 column by the comma
      genes1 <- unlist(strsplit(as.character(dups[row,6]),","))
      genes1 <- trimws(genes1, which = c("both"))
      #add genes1 to the group list
      genes <- append(genes, genes1)
      #splits genes 2 by the comma
      genes2 <- unlist(strsplit(as.character(dups[row,7]),","))
      genes2 <- trimws(genes2, which = c("both"))
      #add genes 2 to list of H genes
      genes <- append(genes, genes2)
    }
    
    head(genes)
    length(genes)
    
## Remove non-species if interest, if necessary ####
    
    nonspecies = "paegeria"
    #For this group, I need to exclude the diulia. I will keep only items in the list without diulia
    #initate empty list
    newlist <- rep(NA, 10)
    for (item in genes) {
      if (str_detect(item,nonspecies) == FALSE) {
        newlist <- append(newlist, item)
      }
    }
    newlist
    final_list <- na.omit(newlist)
    final_list
    length(final_list)
    
    genes <- final_list
    
    
## Option: Only include one type of species ####
    nonspecies = "paegeria"
     newlist <- rep(NA, 10)
      for (item in genes) {
        if (str_detect(item, nonspecies) == TRUE) {
        newlist <- append(newlist, item)
        }
      }
      length(newlist)
      final_list <- na.omit(newlist)
      final_list
      length(final_list)
      genes <- final_list
  
## Write outfile ####
  write.table(genes, outFile, row.names=FALSE, quote = FALSE)
