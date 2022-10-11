
#install.packages("ggpubr")
library(ggpubr)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
hrbrthemes::import_roboto_condensed()

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
  
  #how many total genes?
  nrow(hdf)
  
  #Remove zeros from banynana data
  ba_df <- read.delim(bany_name, header = TRUE, sep = "\t")
  nrow(ba_df)
  ba_df[ba_df == 0] <- NA
  ba_df <- na.omit(ba_df)
  
  #how many total genes?
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
  
  #htest$NI <- (htest$piN / htest$piS) / (htest$ka / htest$ks)
  #btest$NI <- (btest$piN / btest$piS) / (btest$ka / btest$ks)
  htest$NI <- (htest$ks*htest$piN) / (htest$ka*htest$piS)
  btest$NI <- (btest$ks*btest$piN) / (btest$ka*btest$piS)
  
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
  
  
  
## Plot the DOS value for both H and B

  
  htest$Species <- "H.erato"
  btest$Species <- "B.anynana"
  
  ##Plot all genes
  nrow(htest)
  nrow(btest)
  
  b_and_h <- rbind(htest,btest)

  
  # Histogram with all genes and their selective pressures
  p <- b_and_h %>%
    ggplot( aes(x=dos,y=stat(density*width), fill= Species, color = Species)) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',binwidth=.02) +
    scale_fill_manual(values=c("darkblue","cyan")) +
    labs(fill="") +
    geom_vline(xintercept=0, linetype="dotted", color="black")
  
  p
  
  #What percentage of genes are significant for both species?
  
  sig_htest <- htest[htest$p.value < .05, ]  
  sig_btest <- btest[btest$p.value < .05, ]
  
  nrow(sig_htest)
  nrow(sig_htest) / nrow(htest)
  
  nrow(sig_btest)
  nrow(sig_btest) / nrow(btest)
  
  #How many gene are NOT significant?
  nrow(htest[htest$p.value > .05, ])
  nrow(btest[btest$p.value > .05, ])
  nrow(htest[htest$p.value > .05, ]) / nrow(htest)
  nrow(btest[btest$p.value > .05, ]) / nrow(btest)
  
  #How many genes are under positive selection?
  
  nrow(sig_htest[sig_htest$dos > 0, ]) #/ nrow(htest)
  nrow(sig_btest[sig_btest$dos > 0, ]) #/ nrow(btest)
  
  #How many genes are under purifying selection?
  
  nrow(sig_htest[sig_htest$dos < 0, ]) #/ nrow(htest)
  nrow(sig_btest[sig_btest$dos < 0, ]) #/ nrow(btest)
  
  #Relaxed selection?
  #nrow(sig_htest[sig_htest$dos < 0.02 & sig_htest$dos > -0.01, ]) / nrow(htest)
  #nrow(sig_btest[sig_btest$dos < 0.02 & sig_btest$dos > -0.01, ]) / nrow(btest)
  
  
  ##Plot significant genes
  
  #Using y=..density.. scales the histograms so the area under each is 1, or sum(binwidth*y)=1.
  #As a result, you would use y = binwidth*..density..
  #to have y represent the fraction of the total in each bin
  #. In your case, binwidth=0.5.
  sig_b_and_h <- rbind(sig_htest,sig_btest)
  sig_p <- sig_b_and_h %>%
    ggplot( aes(x=dos,y=stat(density*width), fill= Species, color = Species)) +
    geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity',binwidth=.03) +
    scale_fill_manual(values=c("darkblue","cyan")) +
    #theme_ipsum() +
    labs(fill="") +
    geom_vline(xintercept=0, linetype="dotted", color="darkslategrey")
  
  sig_p

  #KS test for distributions of DOS between all genes of both groups
  b = subset(b_and_h, Species == "H.erato")$dos
  h = subset(b_and_h, Species == "B.anynana")$dos
  ks.test(b, h)  

  #KIS test for distribution of DOS between sig genes of both groups
  sigb = subset(sig_b_and_h,Species == "H.erato")$dos
  sigh = subset(sig_b_and_h,Species == "B.anynana")$dos
  ks.test(sigb,sigh)
  
  #FISHER's exact test for distrubution of the number of non-significant (relaxed) vs
  #significant genes experience positive selection
  
  h_relaxed = nrow(htest[htest$p.value > .05, ])
  b_relaxed = nrow(btest[btest$p.value > .05, ])
  
  h_pos = nrow(sig_htest[sig_htest$dos > 0, ])
  b_pos = nrow(sig_btest[sig_btest$dos > 0, ])
  
  #for h. erato, is there a significant difference in genes in relaxed vs. positive selection?
  data_h_relpos <- data.frame(
    "Relaxed Selection" = c(h_relaxed,nrow(htest)),
    "Positive Selection" = c(h_pos,nrow(htest)),
    row.names = c("# of Genes in Category","# of Total Genes")
  )
  data_h_relpos
  
    #chi-squared test
    chisq.test(data_h_relpos)$expected
    chisq.test(data_h_relpos)
    
    
    #fisher's exact test
    ft = fisher.test(data_h_relpos)
    ft
  
  #for b.anynana, is there a significant difference in genes in relaxed vs. positive selection?
  data_b_relpos <- data.frame(
    "Relaxed Selection" = c(b_relaxed,nrow(btest)),
    "Positive Selection" = c(b_pos,nrow(btest)),
    row.names = c("# of Genes in Category","# of Total Genes")
  )
  data_b_relpos
  
    #chi-squared test
    chisq.test(data_b_relpos)$expected
    chisq.test(data_b_relpos)
    #fisher's exact test
    ft = fisher.test(data_b_relpos)
    ft
  
  
  
  
  
  
  #How do these proportions compare between H. erato and B. anynana?
    data <- data.frame(
      "relaxed" = c(h_relaxed,b_relaxed),
      "positive.sel" = c(h_pos,b_pos),
      row.names = c("Herato","Banynana"),
      stringsAsFactors = FALSE
      
    )
    data

  mosaicplot(data,main = "Proportion of ",color = TRUE)
  ft <- fisher.test(data)
  ft
  
  
  #How does the proportion of genes under positive selection compare between H. erato and B. anynana?
  
    data <- data.frame(
      "H.erato" = c(h_pos,nrow(htest)),
      "B.anynana" = c(b_pos,nrow(btest)),
      row.names = c("Positive Selection","All Genes"),
      stringsAsFactors = False
      
    )
    data
    
    ft <- fisher.test(data)
  ft    
  
