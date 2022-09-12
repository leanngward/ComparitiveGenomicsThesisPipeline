setwd("F:/run2_codemlmodels")
library("readxl")
library(ggplot2)
library(reshape2)
library(plyr)
library(qvalue)
#library(svglite)
library(tidyr)


ld = read.delim("branch_crown_output.txt", header = TRUE, sep = "\t")
head(ld)
#Calculate log ratio 
ld$LR <- ld$Alt.lnL - ld$Null.lnL
head(ld)
#Change column names
colnames(ld) <- c("Group", "nulllnl", "nullomega", "altlnl", "bg", "fg", "LR")
#Calculate significance of the log ratio
ld$pval <- pchisq(ld$LR, df = 1, lower.tail = FALSE)
ld$padj <- p.adjust(ld$pval, method = 'fdr')

print("# of genes with omega higher in Heliconius")
head(ld)

#number of omegas where the background is higher than the foreground
nrow(subset(ld,bg  > fg))


#Sort Results by LRT
sorted_ld = ld[order(-ld$LR),]
head(sorted_ld)

#Convert dataframe into long format
wide_ld = subset(sorted_ld, select = c(Group,bg,fg,LR,pval,padj))
head(wide_ld)
long_ld = gather(wide_ld, level, omega, bg, fg)
head(long_ld)

#remove "omega" form bg and fg
long_ld$level <- replace(long_ld$level, long_ld$level=='bg', 'Background')
long_ld$level <- replace(long_ld$level, long_ld$level=='fg', 'Heliconius')


#bin the omegas
c = 3
long_ld$omega_bin <- long_ld$omega
long_ld$omega_bin[long_ld$omega_bin > c] <- c
head(long_ld)

#Graph the results.
mu <- ddply(subset(long_ld), "level", summarise, grp.mean=mean(omega_bin))
mu

avg_bg = mu[1,2]
favg_bg = format(round(avg_bg,2), nsmall=2)

avg_fg = mu[2,2]
favg_fg = format(round(avg_fg,2), nsmall=2)

myplot = ggplot(subset(long_ld), aes(x=omega_bin, fill = level, color=level)) +
  ggtitle("Density Distribution of Omega Estimates for Single-Copy Genes in Heliconius and Background Species") +
  geom_density(alpha = 0.6) +
  xlab('dN/dS') + 
  ylab('Density') +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=level),
             linetype="dashed") +
  #geom_text(aes(x = 1, label = paste0("omega = ",favg_bg), y=10), colour = "red") +
  #geom_text(aes(x = .7, label = paste0("omega = ",favg_fg), y=10), colour = "blue") +
  scale_color_manual(values=c(Background="red", Heliconius="darkblue")) +
  scale_fill_manual(values=c(Background="red", Heliconius="darkblue")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line =   element_line(colour = "black"))+
  xlim(0,.75)

myplot 

print(favg_bg)
#Perform stats - Kolgomorov-Smirnoff test (D=0.10349, P = 4.515e-06).
bg_omega = subset(long_ld, level == "Background")$omega_bin
fg_omega = subset(long_ld, level == "Heliconius")$omega_bin
ks.test(fg_omega, bg_omega)
