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


# RNAseq data from https://drive.google.com/drive/folders/186W8gwrOXspJpYNuS5iWvEkdHv8-mPKm?usp=sharing

chicken <- read.table("~/EPFL/Geneomics and bioinformatics/data/chicken/Chicken.CPM.txt",row.names = 1)
human <- read.table("~/EPFL/Geneomics and bioinformatics/data/human/Human.CPM.txt",row.names = 1)
mouse <- read.table("~/EPFL/Geneomics and bioinformatics/data/mouse/Mouse.CPM.txt",row.names = 1)
opossum <- read.table("~/EPFL/Geneomics and bioinformatics/data/opossum/Opossum.CPM.txt",row.names = 1)
rabbit <- read.table("~/EPFL/Geneomics and bioinformatics/data/rabbit/Rabbit.CPM.txt",row.names = 1)
rat <- read.table("~/EPFL/Geneomics and bioinformatics/data/rat/Rat.CPM.txt",row.names = 1)
rhesus <- read.table("~/EPFL/Geneomics and bioinformatics/data/macaque/Rhesus.CPM.txt",row.names = 1)

# Create the table with samples metadata for figures
colnames_chicken_origin <- paste("chicken", colnames(chicken), sep = "_")
chicken_samples <- data.frame(strsplit(colnames(chicken),".",fixed = TRUE))
chicken_samples2 <- data.frame(t(chicken_samples))
rownames(chicken_samples2) <- colnames_chicken_origin
colnames(chicken_samples2) <- c("Organ","DevStage","Repetition")
chicken_samples2$Specie <- "chicken"

colnames_human_origin <- paste("human", colnames(human), sep = "_")
human_samples <- data.frame(strsplit(colnames(human),".",fixed = TRUE))
human_samples2 <- data.frame(t(human_samples))
rownames(human_samples2) <- colnames_human_origin
colnames(human_samples2) <- c("Organ","DevStage","Repetition")
human_samples2$Specie <- "human"

#mouse info need some cleaning because the separator . is the same as decimal symbol for its developmental stages
colnames(mouse) <- gsub('(.*)(\\.e)(.*)(\\.)(.\\..)', '\\1\\2\\3\\,\\5', colnames(mouse))
colnames_mouse_origin <- paste("mouse", colnames(mouse), sep = "_")
mouse_samples <- data.frame(strsplit(colnames(mouse),".",fixed = TRUE))
mouse_samples2 <- data.frame(t(mouse_samples))
rownames(mouse_samples2) <- colnames_mouse_origin
colnames(mouse_samples2) <- c("Organ","DevStage","Repetition")
mouse_samples2$Specie <- "mouse"

colnames(opossum) <- gsub('(.*)(\\.)(.*)(\\.)(.\\..)', '\\1\\2\\3\\,\\5', colnames(opossum))
colnames_opossum_origin <- paste("opossum", colnames(opossum), sep = "_")
opossum_samples <- data.frame(strsplit(colnames(opossum),".",fixed = TRUE))
opossum_samples2 <- data.frame(t(opossum_samples))
rownames(opossum_samples2) <- colnames_opossum_origin
colnames(opossum_samples2) <- c("Organ","DevStage","Repetition")
opossum_samples2$Specie <- "opossum"

colnames(rabbit) <- gsub('(.*)(\\.e)(.*)(\\.)(.\\..)', '\\1\\2\\3\\,\\5', colnames(rabbit))
colnames_rabbit_origin <- paste("rabbit", colnames(rabbit), sep = "_")
rabbit_samples <- data.frame(strsplit(colnames(rabbit),".",fixed = TRUE))
rabbit_samples2 <- data.frame(t(rabbit_samples))
rownames(rabbit_samples2) <- colnames_rabbit_origin
colnames(rabbit_samples2) <- c("Organ","DevStage","Repetition")
rabbit_samples2$Specie <- "rabbit"

rat_samples <- data.frame(strsplit(colnames(rat),".",fixed = TRUE))
colnames_rat_origin <- paste("rat", colnames(rat), sep = "_")
rat_samples2 <- data.frame(t(rat_samples))
rownames(rat_samples2) <- colnames_rat_origin
colnames(rat_samples2) <- c("Organ","DevStage","Repetition")
rat_samples2$Specie <- "rat"

rhesus_samples <- data.frame(strsplit(colnames(rhesus),".",fixed = TRUE))
colnames_rhesus_origin <- paste("rhesus", colnames(rhesus), sep = "_")
rhesus_samples2 <- data.frame(t(rhesus_samples))
rownames(rhesus_samples2) <- colnames_rhesus_origin
colnames(rhesus_samples2) <- c("Organ","DevStage","Repetition")
rhesus_samples2$Specie <- "rhesus"

samples <- rbind(chicken_samples2,human_samples2,mouse_samples2,opossum_samples2,rabbit_samples2,rat_samples2,rhesus_samples2)


# Gene mapping from supp xlsx file
#install.packages("readxl")
library(readxl)

supp_info_14 <- read_excel("~/EPFL/Geneomics and bioinformatics/41586_2019_1338_MOESM2_ESM.xlsx",
                         sheet = "SI Table 14",skip=2)

supp_info_15 <- read_excel("~/EPFL/Geneomics and bioinformatics/41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 15",skip=2)

supp_info_16 <- read_excel("~/EPFL/Geneomics and bioinformatics/41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 16",skip=2)

supp_info_17 <- read_excel("~/EPFL/Geneomics and bioinformatics/41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 17",skip=2)

supp_info_18 <- read_excel("~/EPFL/Geneomics and bioinformatics/41586_2019_1338_MOESM2_ESM.xlsx",
                           sheet = "SI Table 18",skip=2)


supp_map <- rbind(supp_info_14[supp_info_14$classification=='Conserved',],
                  supp_info_15[supp_info_15$classification=='Conserved',],
                  supp_info_16[supp_info_16$classification=='Conserved',],
                  supp_info_17[supp_info_17$classification=='Conserved',],
                  supp_info_18[supp_info_18$classification=='Conserved',])

supp_map <- supp_map[,c("Human_ID","Mouse_ID","Rat_ID","Rabbit_ID","Opossum_ID")]
supp_map <- unique(supp_map)

# Get corresponding gene ids for chicken and rhesus manually on http://www.ensembl.org/biomart/martview uploading the file "data/ids_mapping/human_genes_conserved_from_paper.txt"
# get gene names id first
human_mapping_ids <- read.table("~/EPFL/Geneomics and bioinformatics/data/ids_mapping/human_mapping_ids.txt",sep='\t',header = 1)

chicken_ids <- read.table("~/EPFL/Geneomics and bioinformatics/data/ids_mapping/chicken_ids_gn.txt",sep='\t',header = 1)
rhesus_ids <- read.table("~/EPFL/Geneomics and bioinformatics/data/ids_mapping/rhesus_ids_gn.txt",sep='\t',header = 1)

chicken_rhesus <- merge(chicken_ids,rhesus_ids, by="Gene.name")

chicken_rhesus_human <- merge(chicken_rhesus,human_mapping_ids, by="Gene.name")

#here is our complete mapping of genes
all_genes <- merge(supp_map,chicken_rhesus_human, by.x="Human_ID",by.y="Gene.stable.ID")

colnames(all_genes)

chicken_with_common_ids <- merge(all_genes[,c("Human_ID","Gene.stable.ID.x")],
                                 chicken,by.x = "Gene.stable.ID.x",by.y = "row.names")
chicken_with_common_ids$Gene.stable.ID.x <- NULL
chicken_with_common_ids <- unique(chicken_with_common_ids)
chicken_with_common_ids <- chicken_with_common_ids[!duplicated(chicken_with_common_ids$Human_ID),]
rownames(chicken_with_common_ids) <- chicken_with_common_ids$Human_ID
chicken_with_common_ids$Human_ID <- NULL

mouse_with_common_ids <- merge(all_genes[,c("Human_ID","Mouse_ID")],
                               mouse,by.x = "Mouse_ID",by.y = "row.names")
mouse_with_common_ids$Mouse_ID <- NULL
mouse_with_common_ids <- unique(mouse_with_common_ids)
mouse_with_common_ids <- mouse_with_common_ids[!duplicated(mouse_with_common_ids$Human_ID),]
rownames(mouse_with_common_ids) <- mouse_with_common_ids$Human_ID
mouse_with_common_ids$Human_ID <- NULL

opossum_with_common_ids <- merge(all_genes[,c("Human_ID","Opossum_ID")],
                                 opossum,by.x = "Opossum_ID",by.y = "row.names")
opossum_with_common_ids$Opossum_ID <- NULL
opossum_with_common_ids <- unique(opossum_with_common_ids)
opossum_with_common_ids <- opossum_with_common_ids[!duplicated(opossum_with_common_ids$Human_ID),]
rownames(opossum_with_common_ids) <- opossum_with_common_ids$Human_ID
opossum_with_common_ids$Human_ID <- NULL

rabbit_with_common_ids <- merge(all_genes[,c("Human_ID","Rabbit_ID")],
                                rabbit,by.x = "Rabbit_ID",by.y = "row.names")
rabbit_with_common_ids$Rabbit_ID <- NULL
rabbit_with_common_ids <- unique(rabbit_with_common_ids)
rabbit_with_common_ids <- rabbit_with_common_ids[!duplicated(rabbit_with_common_ids$Human_ID),]
rownames(rabbit_with_common_ids) <- rabbit_with_common_ids$Human_ID
rabbit_with_common_ids$Human_ID <- NULL


rat_with_common_ids <- merge(all_genes[,c("Human_ID","Rat_ID")],
                             rat,by.x = "Rat_ID",by.y = "row.names")
rat_with_common_ids$Rat_ID <- NULL
rat_with_common_ids <- unique(rat_with_common_ids)
rat_with_common_ids <- rat_with_common_ids[!duplicated(rat_with_common_ids$Human_ID),]
rownames(rat_with_common_ids) <- rat_with_common_ids$Human_ID
rat_with_common_ids$Human_ID <- NULL

rhesus_with_common_ids <- merge(all_genes[,c("Human_ID","Gene.stable.ID.y")],
                                rhesus,by.x = "Gene.stable.ID.y",by.y = "row.names")
rhesus_with_common_ids$Gene.stable.ID.y <- NULL
rhesus_with_common_ids <- unique(rhesus_with_common_ids)
rhesus_with_common_ids <- rhesus_with_common_ids[!duplicated(rhesus_with_common_ids$Human_ID),]
rownames(rhesus_with_common_ids) <- rhesus_with_common_ids$Human_ID
rhesus_with_common_ids$Human_ID <- NULL

# samples in order: chicken_samples2,human_samples2,mouse_samples2,opossum_samples2,rabbit_samples2,rat_samples2,rhesus_samples2
colnames(chicken_with_common_ids) <- paste("chicken", colnames(chicken_with_common_ids), sep = "_")
colnames(human) <- paste("human", colnames(human), sep = "_")
all_data_1 <- merge(chicken_with_common_ids,human,by="row.names")
rownames(all_data_1) <- all_data_1$Row.names
all_data_1$Row.names <- NULL

colnames(mouse_with_common_ids) <- paste("mouse", colnames(mouse_with_common_ids), sep = "_")
colnames(opossum_with_common_ids) <- paste("opossum", colnames(opossum_with_common_ids), sep = "_")

all_data_2 <- merge(mouse_with_common_ids,opossum_with_common_ids,by="row.names")
rownames(all_data_2) <- all_data_2$Row.names
all_data_2$Row.names <- NULL

colnames(rabbit_with_common_ids) <- paste("rabbit", colnames(rabbit_with_common_ids), sep = "_")
colnames(rat_with_common_ids) <- paste("rat", colnames(rat_with_common_ids), sep = "_")

all_data_3 <- merge(rabbit_with_common_ids,rat_with_common_ids,by="row.names")
rownames(all_data_3) <- all_data_3$Row.names
all_data_3$Row.names <- NULL

colnames(rhesus_with_common_ids) <- paste("rhesus", colnames(rhesus_with_common_ids), sep = "_")

all_data_4 <- merge(all_data_3,rhesus_with_common_ids,by="row.names")
rownames(all_data_4) <- all_data_4$Row.names
all_data_4$Row.names <- NULL

all_data_12 <- merge(all_data_1,all_data_2,by="row.names")
rownames(all_data_12) <- all_data_12$Row.names
all_data_12$Row.names <- NULL

all_data <- merge(all_data_12,all_data_4,by="row.names")
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

all(colnames(all_data) == rownames(samples) )




library(FactoMineR)
library("factoextra")


dt <- data.matrix(all_data)
pca_res <- PCA(t(dt), graph = FALSE)


fviz_pca_ind(pca_res, geom="point", alpha=I(0)) +
  geom_point(aes(shape=samples$Specie, colour=samples$Organ)) + 
  scale_shape_manual(values = c("human"=15, "mouse"=16,"opossum"=17, "rabbit"=18, "chicken"=8, "rhesus"=22, "rat"=21)) +
  scale_colour_manual(breaks = c("Brain","Cerebellum", "Heart", "Kidney", "Liver", "Ovary", "Testis"),
                      values=c("#0494c6","#00cfff", "#e10000", "#db9600", "#009c00", "#e6009d", "#ff5900")) +
  scale_y_continuous(breaks=c(-50,0, 50), limits=c(-55, 55)) + scale_x_continuous(breaks=c(-40,0, 40)) +
  labs(color="Organ", shape="Species",size="Developmental stage", x="PC1", y="PC2",
       title="PCA")

pdf("a4r_output_pca.pdf", paper="a4r")
fviz_pca_ind(pca_res, geom="point", alpha=I(0)) +
  geom_point(aes(shape=samples$Specie, colour=samples$Organ,size=samples$DevStage)) + 
  scale_size_manual(values = c("e10"=0.5,"e12"=0.5,"e14"=1,"e17"=1,"P0"=1,"P7"=1.5,"P35"=1.5,"P70"=2,"P155"=2,
                               "4wpc"=0.25,"5wpc"=0.5,"6wpc"=0.5,"7wpc"=0.5,"8wpc"=1,"9wpc"=1,"10wpc"=1,"11wpc"=1.75,"12wpc"=2,"13wpc"=2,"16wpc"=2,"18wpc"=2.5,"19wpc"=2.5,"20wpc"=3,
                               "newborn"=3,"infant"=3,"toddler"=3,"school"=3,"teenager"=3.5,"youngAdult"=3.5,"youngMidAge"=4,"olderMidAge"=4,"senior"=4,
                               "youngTeenager"=3,"oldTeenager"=3.5,"Senior"=3.5,
                               "e10,5"=1,"e11,5"=1,"e12,5"=1,"e13,5"=1,"e14,5"=1,"e15,5"=2,"e16,5"=2,"e17,5"=3,
                               "e18,5"=3,"P3"=1,"P14"=1,"P28"=1,"P63"=1,
                               "13,5"=1,"14"=2,"16"=3,"18"=2,"20"=3,"24"=3,"28"=3,"35"=4,"42"=4,"56"=4,"74"=5,"104"=5,
                               "134"=4,"164"=4,"194"=4,"e13"=1,"e18"=2,"e19,5"=2,"e21"=2,"e23"=3,"e24"=3,"e27"=3,"P84"=3,"P186"=4,
                               "e11"=0.5,"e15"=0.5,"e16"=1,"e19"=2,
                               "e20"=1,"P42"=3,"P112"=4,"e93"=4,"e108"=4,"e112"=4,"e123"=4,"e130"=5,"P23"=2,"P152"=3,"P365"=3,"P1095"=3,"P3285"=4,"P5475"=4,"P8030"=5,"P183"=3)) +
  scale_shape_manual(values = c("human"=15, "mouse"=16,"opossum"=17, "rabbit"=18, "chicken"=8, "rhesus"=22, "rat"=21)) +
  scale_colour_manual(breaks = c("Brain","Cerebellum", "Heart", "Kidney", "Liver", "Ovary", "Testis"),
                      values=c("#0494c6","#00cfff", "#e10000", "#db9600", "#009c00", "#e6009d", "#ff5900")) +
  scale_y_continuous(breaks=c(-50,0, 50), limits=c(-55, 55)) + scale_x_continuous(breaks=c(-40,0, 40)) +
  labs(color="Organ", shape="Species",size="Developmental stage", x="PC1", y="PC2",
       title="PCA")

dev.off()


