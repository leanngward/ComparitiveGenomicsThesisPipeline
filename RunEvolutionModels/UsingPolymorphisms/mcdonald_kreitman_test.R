
## set working directory ####
  setwd("C:/Users/leann/OneDrive/grad research_newer/butterfly_aging/run5_asproteins_orthofinder/mcdonald_kreitman_data/run2/")

## Set input file names ####

  kaks_name <- "NandSoutput.tsv"
  h1 <- "h_erato_1to5.out"
  h2 <- "h_erato_6to10.out"
  h3 <- "h_erato_11to15.out"
  h4 <- "h_erato_16to19.out"
  h5 <- "h_erato_20to22.out"
  bany_name <- "bany.out"
  
  
## Parse and Clean data
  
  # Herato data needs to be. Only some orthogroups were able to found in some of the data.
  herato1 <- read.delim(h1,sep="\t",header=TRUE)
  herato2 <- read.delim(h2,sep="\t",header=TRUE)
  herato3 <- read.delim(h3,sep="\t",header=TRUE)
  herato4 <- read.delim(h4,sep="\t",header=TRUE)
  herato5 <- read.delim(h5,sep="\t",header=TRUE)
  
  herato1[herato1 == 0] <- NA
  herato2[herato2 == 0] <- NA 
  herato3[herato3 == 0] <- NA 
  herato4[herato4 == 0] <- NA 
  herato5[herato5 == 0] <- NA 

  herato1 <- na.omit(herato1)
  herato2 <- na.omit(herato2)
  herato3 <- na.omit(herato3)
  herato4 <- na.omit(herato4)
  herato5 <- na.omit(herato5)    
  
  #Checked if orthogroup values were unique - PASSED
  #Looks like there's three group orthogroups missing
  hdf <- rbind(herato1,herato2,herato3,herato4,herato5)

  #Remove zeros from banynana data
  ba_df <- read.delim(bany_name, header = TRUE, sep = "\t")
  ba_df[ba_df == 0] <- NA
  ba_df <- na.omit(ba_df)
  nrow(ba_df)
  
  #Read in ka/ks data
  kaks <- read.delim(kaks_name, header = TRUE, sep = "\t")

  
  #Explore which groups are missing between data files
  keep <-  kaks$Group %in% hdf$X.Orthogroup 
  kaks_H <- kaks[which(keep == TRUE), ]
  
  keep2 <- kaks$Group %in% ba_df$X.Orthogroup
  kaks_B<- kaks[which(keep2 == TRUE), ]
  
  #fix for non-numeric error, not sure why it occurred
  kaks_B$B.vs..D.ks <- as.numeric(kaks_B$B.vs..D.ks)

  
## Calculations for each groups ####

  # data frame for Herato vs. Dmelanogaster test
  htest <- data.frame(group = kaks_H$Group, ka = kaks_H$H.vs..D.ka, ks = kaks_H$H.vs..D.kS, piN = hdf$Non.synonymous.SNPs, piS = hdf$Synonymous.SNPs)
  
  # data frame for Banynana vs. Dmelanogaster test
  btest <- data.frame(group = kaks_B$Group, ka = kaks_B$B.vs..D.ka, ks = kaks_B$B.vs..D.ks, piN = ba_df$Non.synonymous.SNPs, piS = ba_df$Synonymous.SNPs)
  
  ## HEY CHECK PN VS PS AND GET DIFFERENT N VS S values
  # Calculate the neutrality index ( [pN/pS] / [ka/ks])
  
  htest$NI <- (htest$piN / htest$piS) / (htest$ka / htest$ks)
  btest$NI <- (btest$piN / btest$piS) / (btest$ka / btest$ks)
  

  
  #?fisher.test
  
  #data <- c(htest[1,2], htest[1,4], htest[1,3], htest[1,5])
  #print(data)
  #mat <- matrix(data,nrow=2,ncol=2,byrow=TRUE)
  #print(mat)
  #ft <- fisher.test(mat)
  #?fisher.test
  #ft$p.value

 for (row in 1:nrow(htest)) {
    data <- c(htest[row,2], htest[row,4], htest[row,3], htest[row,5])
    mat <- matrix(data,nrow=2,ncol=2,byrow=TRUE)
    
    ft <- fisher.test(mat)
    htest[row,7] <- as.numeric(ft$p.value)
 }


  for (row in 1:nrow(btest)) {
    #create contigency matrix
    data <- c(btest[row,2], btest[row,4], btest[row,3], btest[row,5])
    mat <- matrix(data,nrow=2,ncol=2,byrow=TRUE)
    
    #test if the ratios are significantly different from each other
    ft <- fisher.test(mat)
    btest[row,7] <- as.numeric(ft$p.value)
  }  
  
  colnames(htest)[7] <- "p.value"
  colnames(btest)[7] <- "p.value"

  #Calculate Direction of Selection (an alternative to the Neutrality Index)
  # DOS > 1 is positive selection, = 0 is neutral selection, < 1 purifying selection
  htest$dos <- htest$ka / (htest$ka+htest$ks) - htest$piN / (htest$piN + htest$piS) 
  btest$dos <- btest$ka / (btest$ka + btest$ks) - btest$piN / (btest$piN + btest$piS)
  
  htest$alpha <- 1 - ((htest$piN/htest$piS)/(htest$ka/htest$ks))
  btest$alpha <- 1 - ((btest$piN/btest$piS)/(btest$ka/btest$ks))
  
  sig_htest <- htest[htest$p.value < .05, ]  
  sig_btest <- btest[btest$p.value < .05, ]
  
  nrow(sig_btest)  
  nrow(sig_htest)  
  print(sig_htest[sig_htest$dos > 1, ])
  print(sig_btest[sig_btest$dos > 1, ])
  