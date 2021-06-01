library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots)  
library(tidyr)
library(stringi)


# RNAseq data from https://drive.google.com/drive/folders/186W8gwrOXspJpYNuS5iWvEkdHv8-mPKm?usp=sharing

human <- read.table("~/EPFL/Geneomics and bioinformatics/data/human/Human.CPM.txt",row.names = 1)

colnames_human_origin <- paste("human", colnames(human), sep = "_")
human_samples <- data.frame(strsplit(colnames(human),".",fixed = TRUE))
human_samples2 <- data.frame(t(human_samples))
rownames(human_samples2) <- colnames_human_origin
colnames(human_samples2) <- c("Organ","DevStage","Repetition")
human_samples2$Specie <- "human"

install.packages("Hmisc")

library(Hmisc)

res2<-rcorr(as.matrix(human),type="spearman")  # this takes some time


# from http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software 
# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

fcm <- flattenCorrMatrix(res2$r, res2$P)

samples <- data.frame(transpose(strsplit(fcm$row,".",fixed = TRUE)))
colnames(samples) <- c("organ","age", "sample")

samples2 <- data.frame(transpose(strsplit(fcm$column,".",fixed = TRUE)))
colnames(samples2) <- c("organ2","age2", "sample2")

fulldf <- cbind(fcm,samples,samples2)

# to keep only the same dev stage comparisons:
fulldf <- fulldf[fulldf$age==fulldf$age2,]

fulldf <- fulldf[fulldf$organ=="Brain",]
fulldf <- fulldf[fulldf$organ2 != "Brain",] # to remove brain-brain corr

unique(fulldf$age)


fulldf$ageNum <- NA

fulldf[fulldf$age=="4wpc",]$ageNum <- 1
fulldf[fulldf$age=="5wpc",]$ageNum <- 2
fulldf[fulldf$age=="7wpc",]$ageNum <- 3
fulldf[fulldf$age=="8wpc",]$ageNum <- 4
fulldf[fulldf$age=="9wpc",]$ageNum <- 5
fulldf[fulldf$age=="10wpc",]$ageNum <- 6
fulldf[fulldf$age=="11wpc",]$ageNum <- 7
fulldf[fulldf$age=="12wpc",]$ageNum <- 8
fulldf[fulldf$age=="13wpc",]$ageNum <- 9
fulldf[fulldf$age=="16wpc",]$ageNum <- 10

fulldf[fulldf$age=="18wpc",]$ageNum <- 11
fulldf[fulldf$age=="19wpc",]$ageNum <- 12
fulldf[fulldf$age=="20wpc",]$ageNum <- 13

fulldf[fulldf$age=="newborn",]$ageNum <- 14
fulldf[fulldf$age=="infant",]$ageNum <- 15
fulldf[fulldf$age=="toddler",]$ageNum <- 16
fulldf[fulldf$age=="school",]$ageNum <- 17
fulldf[fulldf$age=="teenager",]$ageNum <- 18
fulldf[fulldf$age=="youngAdult",]$ageNum <- 19
fulldf[fulldf$age=="youngMidAge",]$ageNum <- 20
fulldf[fulldf$age=="olderMidAge",]$ageNum <- 21
fulldf[fulldf$age=="senior",]$ageNum <- 22



library(ggplot2)

ggplot(fulldf, aes(ageNum,cor, col=organ2)) + ylim(0.7,1) + 
  geom_point() + 
  labs(x ="Developmental stages", y="Transcriptome similarity (Spearman's correlation)") +
  coord_fixed(ratio=80) +
  scale_x_continuous(breaks=c(1,5,9,14,19,22), labels=c("4w","9w","13w","new","ya","sen")) +
  labs(color="Organ")


