
################################################################################################
## Set Up Environment and Variables ####
  #Set working directory
    setwd("C:/Users/leann/OneDrive/grad research_newer/butterfly_aging/run7/longevitygenes")
  #Set up strings for output tables
    species1 = "Heliconius"
    category1 = "Clade"
    species2 = "Pieridae"
    category2 = "Clade"

# OMEGA ANALYSIS
  #Set evolution data (from branch model) file name

  #CTL-F REPLACE EITHER 'pier' or 'pier' with the other

    heli_infile <- "hclade_edited_branch_output.tsv"
    pier_infile <- "edited_pclade_branch_output.tsv"

#Provide lists of longevity genes mapped to the groups of interest (OrthoGroups)

    group1list <- "genelists/longevity_and_metabolism_hits.txt"
    #group1list <- "genelists/control_runhits.tsv"

#########################################################################################################################    
## Open Files and Find the subset of orthogroups in the longevity list - IF relevant ####

    #Open longevity genes list
    #NOTE: The species and orthogroups need to be split by the "_" ... I did this in Excel

      group1 <- read.delim(file=group1list,header=FALSE,sep="\t")

    #OMEGA ANALYSIS
    #branch model
      evodata <- read.delim(heli_infile,header=TRUE,sep="\t")
      heli_evodata <- read.delim(heli_infile, header=TRUE,sep="\t")
      pier_evodata <- read.delim(pier_infile, header=TRUE,sep="\t")

      head(heli_evodata)
      head(pier_evodata)

###########################################################################################################################    

#OMEGA ANALYSIS
  ## Find significant omega values ####
  # Calculate log-likelihood ratio

    heli_evodata$LR <- -2 * (heli_evodata$Null.lnL - heli_evodata$Alt.lnL)
    pier_evodata$LR <- -2 * (pier_evodata$Null.lnL - pier_evodata$Alt.lnL)

  #Calculate p-value
    heli_evodata$pval <- pchisq(heli_evodata$LR, df = 1, lower.tail = FALSE)
    pier_evodata$pval <- pchisq(pier_evodata$LR, df = 1, lower.tail = FALSE)

  #adjust p-value to control false discovery rate
    heli_evodata$padj <- p.adjust(heli_evodata$pval, method = 'fdr')
    pier_evodata$padj <- p.adjust(pier_evodata$pval, method = 'fdr')

  #create table with only significant results
    heli_sigevodata <- heli_evodata[heli_evodata$padj < .05, ]
    pier_sigevodata <- pier_evodata[pier_evodata$padj < .05, ]
    nrow(heli_sigevodata)
    nrow(pier_sigevodata)



## Find when the Foreground omega > Background Omega and vice versa

  #OMEGA ANALYSIS
      heli_sigevodata$testcol <- heli_sigevodata$Alt.Omega.1 > heli_sigevodata$Alt.Omega.0.Background
      pier_sigevodata$testcol <- pier_sigevodata$Alt.Omega.1 > pier_sigevodata$Alt.Omega.0.Background


## GET THE NUMBER OF LONGEVITY GENES in ALL OF THE GENE SET ####

    unique.group1 <- unique(group1[ , c(1,3)])

    group1_allgenes <- merge(x=unique.group1,y=evodata, by.x=c("V3"), by.y=c("Group"))
    print(group1_allgenes)

    cnt_group1_allgenes <- unique(group1_allgenes[c(1,2)])



  #list of all longevity genes detected for each genes for all genes
    df_allgenes <- cnt_group1_allgenes
    cnt_allgenes <- nrow(df_allgenes)

  #vector of TRUE/FALSE if genes is a longevity genes
    allgenes_longevity_vector <- evodata$Group %in% df_allgenes$V3
    length(allgenes_longevity_vector)

## Explore longevity genes for OMEGA ANALYSIS .... ####

  #Using the top and bottom 10% of values, cross reference which genes are in these lists

  top10_heli_omega <- subset(heli_sigevodata, heli_sigevodata$testcol == TRUE)
  top10_pier_omega <- subset(pier_sigevodata, pier_sigevodata$testcol == TRUE)

  bottom10_heli_omega <- subset(heli_sigevodata,heli_sigevodata$testcol == FALSE)
  bottom10_pier_omega <- subset(pier_sigevodata,pier_sigevodata$testcol == FALSE)

  ## group1 ... ####

    ## unique hits of genes/orthogroup

    unique.group1 <- unique(group1[ , c(1,3)])

  ## mapped hits
    group1_hits10_heli <- merge(x=unique.group1,y=top10_heli_omega,by.x=c("V3"), by.y=c("Group"))
    group1_hits90_heli <- merge(x=unique.group1,y=bottom10_heli_omega,by.x=c("V3"), by.y=c("Group"))

    group1_hits10_pier <- merge(x=unique.group1,y=top10_pier_omega,by.x=c("V3"), by.y=c("Group"))
    group1_hits90_pier <- merge(x=unique.group1,y=bottom10_pier_omega,by.x=c("V3"), by.y=c("Group"))

    print(group1_hits10_heli)
    print(group1_hits90_heli)
    print(group1_hits10_pier)
    print(group1_hits90_pier)

    cnt_group1_top10_heli <- unique(group1_hits10_heli[c(2)])
    cnt_group1_bot10_heli <- unique(group1_hits90_heli[c(2)])
    cnt_group1_top10_pier <- unique(group1_hits10_pier[c(2)])
    cnt_group1_bot10_pier <- unique(group1_hits90_pier[c(2)])


## Count how many longevity genes were found in each category for each groups: ####

  # FG > BG HELICONIUS

    df_top10heli <- cnt_group1_top10_heli
    cnt_top10heli <- nrow(df_top10heli)
    print(cnt_top10heli)
    
    outname_top10heli = paste("fghigher_longevitygenes_found",species1,category1,sep="_")
    write.table(df_top10heli,file=outname_top10heli,sep="\t",quote = FALSE, row.names = FALSE)
    

  # BG > FG HELICONIUS

    df_bot10heli <- cnt_group1_bot10_heli
    cnt_bot10heli <- nrow(df_bot10heli)
    print(cnt_bot10heli)
    
    outname_bot10heli = paste("bggreater_longevitygenes_found",species1,category1,sep="_")
    write.table(df_bot10heli,file=outname_bot10heli,sep="\t",quote = FALSE, row.names = FALSE)
    

  # FG > BG pierIDAE

    df_top10pier <- cnt_group1_top10_pier
    cnt_top10pier <- nrow(df_top10pier)
    print(cnt_top10pier)
    
    
    outname_top10pier = paste("fggreater_longevitygenes_found",species2,category2,sep="_")
    write.table(df_top10pier,file=outname_top10pier,sep="\t",quote = FALSE, row.names = FALSE)
    

  # BG > FG of pierIDAE
    df_bot10pier <- cnt_group1_bot10_pier
    cnt_bot10pier <- nrow(df_bot10pier)
    print(cnt_bot10pier)
    
    outname_bot10pier = paste("bggreater_longevitygenes_found",species2,category2,sep="_")
    write.table(df_bot10pier,file=outname_bot10pier,sep="\t",quote = FALSE, row.names = FALSE)
    



## What's the probability of sampling a certain number of longevity genes in our groups of interest versus sampling them from the total number of genes?


### FOR ALL TESTS
#Vector to sample from
    vec <- allgenes_longevity_vector

#Number of simulations
    NB_simulations <- 3000000


### TEST OF FG > BG OMEGA IN HELICONIUS

  #Get sample size - Actual number of genes in the top10% of omega rates - this will be the size of our random samples
    nbhelitop10 = nrow(top10_heli_omega)

    print("Actual number of longevity genes sampled for Heliconius ")
    nbtop_heli = cnt_top10heli
    print(nbtop_heli)

#vector to store number of randomly sampled genes
    long_found = rep(NA,NB_simulations)

    for (i in (1:NB_simulations)) {
      vecsample = sample(vec,nbhelitop10,replace = F) #sample longevity genes from all genes (TRUE = longevity, FALSE = NOT LONGEVITY)
      long_found[i] = length(which(vecsample == TRUE))
    }

    hist(long_found, nclass = 50, col = "blueviolet", main = "Distribution of Longevity Genes Sampled",xlab = "Number of longevity genes Found")

    nbtop_heli = cnt_top10heli
    print("Actual number of longevity genes sampled for Heliconius ")
    print(nbtop_heli)

# p-value == the number of observations above/below the test statistic / the number of all observations

  #right-tail
    right_nbSuccess = length(which(long_found>=nbtop_heli))
    left_nbSuccess = length(which(long_found<=nbtop_heli))

    pvalue <- right_nbSuccess / NB_simulations
    print(pvalue)
    
    l_pvalue <- left_nbSuccess / NB_simulations
    print(1-pvalue)

#### TEST OF BG > FG OMEGA IN HELICONIUS

  #Get sample size - the number of genes with significantly low omegas in Heliconius
    nbhelibot10 = nrow(bottom10_heli_omega)

  #vector to store number of randomly sampled genes
    long_found = rep(NA,NB_simulations)

  #Simulate random samples 
    for (i in (1:NB_simulations)) {
      vecsample = sample(vec,nbhelibot10,replace = F) #sample longevity genes from all genes (TRUE = longevity, False = NOT LONGEVITY)
      long_found[i] = length(which(vecsample == TRUE))
    }

    hist(long_found, nclass = 50, col = "blueviolet", main = "Distribution of Longevity Genes Sampled",xlab = "Number of longevity genes Found")


    bot10_long_heli = cnt_bot10heli
    print("Actual number of longevity genes sampled for Heliconius")
    print(bot10_long_heli)

    # p-value == the number of observations above/below the test statistic / the number of all observations

    #left-tail
      left_nbSuccess = length(which(long_found<bot10_long_heli))
    #right-tail
      right_nbSuccess = length(which(long_found>=bot10_long_heli))

      pvalue <- right_nbSuccess / NB_simulations
      print(pvalue)
      l_pvalue = left_nbSuccess / NB_simulations

### TEST OF fg > bg OMEGA IN pierIDAE

    #Get sample size - Actual number of genes in the top10% of omega rates - this will be the size of our random samples
      nbpiertop10 = nrow(top10_pier_omega)

      print("Actual number of longevity genes sampled for pieridae")
      nbtop_pier = cnt_top10pier
      print(nbtop_pier)

    #vector to store number of randomly sampled genes
      long_found = rep(NA,NB_simulations)

      for (i in (1:NB_simulations)) {
        vecsample = sample(vec,nbpiertop10,replace = F) #sample longevity genes from all genes (TRUE = longevity, False = NOT LONGEVITY)
        long_found[i] = length(which(vecsample == TRUE))
      }

      hist(long_found, nclass = 50, col = "blueviolet", main = "Distribution of Longevity Genes Sampled",xlab = "Number of longevity genes Found")

      nbtop_pier = cnt_top10pier
      print("Actual number of longevity genes sampled for pieridae ")
      print(nbtop_pier)

  #calculate p-value

  #right-tail
    right_nbSuccess = length(which(long_found>=nbtop_pier))
    left_nbSuccess = length(which(long_found<=nbtop_pier))
    pvalue <- right_nbSuccess / NB_simulations
    print(pvalue)

  #if you need the left-tail p-value
    l_pvalue <- left_nbSuccess / NB_simulations
    print(l_pvalue)

### TEST OF BG > fg IN pierIDAE

  #Get sample size - Actual number of genes in the top10% of omega rates - this will be the size of our random samples
    nbpierbot10 = nrow(bottom10_pier_omega)

    print("Actual number of longevity genes sampled for pieridae")
    nbbot_pier = cnt_bot10pier
    print(nbbot_pier)

    #vector to store number of randomly sampled genes
    long_found = rep(NA,NB_simulations)

    for (i in (1:NB_simulations)) {
      vecsample = sample(vec,nbpierbot10,replace = F) #sample longevity genes from all genes (TRUE = longevity, False = NOT LONGEVITY)
      long_found[i] = length(which(vecsample == TRUE))
    }

      hist(long_found, nclass = 50, col = "blueviolet", main = "Distribution of Longevity Genes Sampled",xlab = "Number of longevity genes Found")

    nbbot_pier = cnt_bot10pier
    print("Actual number of longevity genes sampled for pieridae ")
    print(nbbot_pier)

  #calculate p-value

  #right-tail
    right_nbSuccess = length(which(long_found>=nbbot_pier))
    left_nbSuccess = length(which(long_found<=nbbot_pier))
    pvalue <- right_nbSuccess / NB_simulations
    l_pvalue <- left_nbSuccess / NB_simulations
    print(pvalue)

  #if you need the left-tail p-value
    print(l_pvalue)

