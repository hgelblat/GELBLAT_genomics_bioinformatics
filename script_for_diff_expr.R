#https://f1000research.com/articles/5-1408/v3

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")

# load pkgs
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(dplyr)
library(gplots)  
library(tidyr)
library(stringi)

library(limma) 
library(edgeR) 
human <- read.table("~/EPFL/Geneomics and bioinformatics/data/human/Human.CPM.txt", quote="\"", comment.char="")
rhesus <- read.table("~/EPFL/Geneomics and bioinformatics/data/macaque/Rhesus.CPM.txt", quote="\"", comment.char="")


# Gene mapping from supp xlsx file
#install.packages("readxl") #run it only once. It installs pkg to read xlsx file (or you can copy content to txt and use read.table of fread):
library(readxl)

supp_info_14 <- read_excel("41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 14",skip=2)

supp_info_15 <- read_excel("41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 15",skip=2)

supp_info_16 <- read_excel("41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 16",skip=2)

supp_info_17 <- read_excel("41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 17",skip=2)

supp_info_18 <- read_excel("41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 18",skip=2)


supp_map <- rbind(supp_info_14[supp_info_14$classification=='Conserved',],
                  supp_info_15[supp_info_15$classification=='Conserved',],
                  supp_info_16[supp_info_16$classification=='Conserved',],
                  supp_info_17[supp_info_17$classification=='Conserved',],
                  supp_info_18[supp_info_18$classification=='Conserved',])

supp_map <- supp_map[,c("Human_ID","Mouse_ID","Rat_ID","Rabbit_ID","Opossum_ID")]
supp_map <- unique(supp_map)

human_mapping_ids <- read.table("data/ids_mapping/human_mapping_ids.txt",sep='\t',header = 1)
chicken_ids <- read.table("data/ids_mapping/chicken_ids_gn.txt",sep='\t',header = 1)
rhesus_ids <- read.table("data/ids_mapping/rhesus_ids_gn.txt",sep='\t',header = 1)
chicken_rhesus <- merge(chicken_ids,rhesus_ids, by="Gene.name")
chicken_rhesus_human <- merge(chicken_rhesus,human_mapping_ids, by="Gene.name")

#here is our complete mapping of genes
all_genes <- merge(supp_map,chicken_rhesus_human, by.x="Human_ID",by.y="Gene.stable.ID")

hum_vs_rhesus <- all_genes[,c("Human_ID","Gene.stable.ID.y")]
colnames(hum_vs_rhesus) <- c("Gene stable ID", "Rhesus_ID")

rhesus <- merge(rhesus, hum_vs_rhesus, by.x = "row.names",by.y = "Rhesus_ID",sort = TRUE)
human <- merge(human, hum_vs_rhesus, by.x = "row.names",by.y = "Gene stable ID",sort = TRUE)

colnames(human) <- paste("human", colnames(human), sep = "_")
colnames(rhesus) <- paste("rhesus", colnames(rhesus), sep = "_")

all_data <- merge(human,rhesus, by.x = "human_Row.names",by.y ='rhesus_Gene stable ID',sort = TRUE)

human_samples <- colnames(human)[2:298] # if you choose other species, change the number in [2:296]
rhesus_samples <- colnames(rhesus)[2:166] # if you choose other species, change the number in [2:296] 

brain_data <- all_data[,(grepl( "Brain" , names( all_data ) ))|(grepl( "human_Row.names" , names( all_data ) ))]

rnaseq_mat <- as.matrix(brain_data[,grepl( "Brain" , names( brain_data ) )]) 

row.names(rnaseq_mat) <- brain_data$human_V1

colnames(rnaseq_mat) 

sample_types <- cbind(data.table(colnames(rnaseq_mat)), data.table(colnames(rnaseq_mat)) %>%
  separate(V1, c("species", "metadata"), "_") %>%
  separate(metadata, c("organ", "age", "sample")))

x <- DGEList(counts=rnaseq_mat,samples=sample_types, group=sample_types$species )
class(x)

dim(x)

samplenames <- colnames(x) #substring(colnames(x), 1, nchar(colnames(x))-3)  # if you want shorter names
#colnames(x) <- samplenames


group <- x$samples$group 

##########################
#Data pre-processing
##########################

#Popular transformations include counts per million (CPM), log2-counts per million (log-CPM), 
#reads per kilobase of transcript per million (RPKM), and fragments per kilobase of transcript per million (FPKM).

#In our analyses, CPM and log-CPM transformations are used regularly although 
#they do not account for gene length differences as RPKM and FPKM values do. 

#gene lengths remain constant for comparisons of interest and any observed differences 
#are a result of changes in condition rather than changes in gene length.

#cpm function in edgeR
cpm <- x # we have CPM in original data
lcpm <- cpm(x, log=TRUE)

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M) 

summary(lcpm)


################
#Removing genes that are lowly expressed
################

keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

#for Fig1 from Tutorial
lcpm.cutoff <- log2(10/M + 2/L)
library(RColorBrewer)
nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.5), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
#legend("topright", samplenames, text.col=col, bty="n")

###############
#Normalising gene expression distributions
###############

#calcNormFactors function in edgeR


x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

#########

##############
#Unsupervised clustering of samples
##############

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,1))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
plotMDS(lcpm, labels=group, col=col.group)
title(main="Sample groups")


#####################
#Differential expression analysis
#Creating a design matrix and contrasts
#####################

design <- model.matrix(~0+group)
colnames(design) <- gsub("group", "", colnames(design))
design
typeof(human)

contr.matrix <- makeContrasts(
  HvsR = human-rhesus,
  levels = colnames(design))
contr.matrix

##################
#Removing heteroscedascity from count data
##################
v <- voom(x, design, plot=TRUE)  # generates fig 4a voom mean-variance trend
#v

##################
#Fitting linear models for comparisons of interest
##################

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
dte <- decideTests(efit)
plotSA(efit) # generates fig 4b final  mean-variance trend

################
#Examining the number of DE genes
################

summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
#summary(dt)

de.common <- which(dt[,1]!=0)
length(de.common)


###############
#Examining individual DE genes from top to bottom
###############

h.vs.r <- topTreat(tfit, coef=1, n=Inf)
#head(h.vs.r)

###############
#Useful graphical representations of differential expression results
###############
plotMD(efit, column=1, status=dte[,1], main=colnames(efit)[1], xlim=c(0,13)) 
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(0,13)) 

#dt[,1]

df <- data.table(h.vs.r)
df$genenames <- rownames(h.vs.r)

head(df[with(df, order(adj.P.Val,logFC)), ])

df.down <- df[with(df, order(adj.P.Val,logFC)),]

downgenes <- df.down[df.down$adj.P.Val<0.001 & df.down$logFC < -1 ]$genenames

df.up <- df[with(df, order(adj.P.Val,-logFC)),]
head(df.up)

upgenes <- df.up[df.up$adj.P.Val<0.001 & df.up$logFC > 1 ]$genenames
length(upgenes)
print(upgenes[0:100])

length(downgenes)
print(
  downgenes
)

write.table(downgenes,"data/ids_mapping/downgenes.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)



write.table(upgenes,"data/ids_mapping/downgenes.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)

