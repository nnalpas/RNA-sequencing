#######################################################
# RNA-seq data analysis of sense counts (paired data) #
#######################################################

# Analysis nearly identical to previous one, with the exception that for article purpose all CN0H samples have been removed as well as all TB samples

#############################
# List of required packages #
#############################

library(edgeR)
library(biomaRt)
library(MASS)

######################################################
# Use featureCounts output files as input files in R #
######################################################

###############
# Preparation #
###############

# Move to the appropriate folder
setwd("C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Sense_gene")
getwd()
workDir <- getwd()
workDir

################################################
# Read in and concatenate input files within R #
################################################

# Create vector of all files name
fileDir <- "C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Sense_gene/Counts"
files <- list.files(path=fileDir, pattern="*(CN|MB)_(2|6|24|48)H$", all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE)
files

# Reads and merges a set of files containing counts
Count <- readDGE(files=files, path=fileDir, columns=c(1,3))
names(Count)
head(Count$samples)
head(Count$counts)

# Ouptut samples data
write.table(x=Count$samples, file="Alv_samples.txt",  sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Count$counts, file="Alv_rawcounts.txt",  sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

######################################################################
# Gene annotation using custom sense annotation info and biomaRt #
######################################################################

# Read in custom annotation information for all sense genes
sense_info <- read.table(file="gene_info_251113.txt", sep="\t", header=TRUE, fill=TRUE)
sense_info <- sense_info[grep(pattern="^ENSBTAG", x=sense_info[,1], perl=TRUE),]
colnames(sense_info)[1] <- "ensembl_gene_id"
colnames(sense_info)[7] <- "ensembl_transcript_id"
dim(sense_info)
head(sense_info)

# Select biomaRt database and dataset (use biomaRt from Ensembl 71: Apr 2013)
listMarts(host='apr2013.archive.ensembl.org')
mart <- useMart(host='apr2013.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL")
#listMarts()
#mart <- useMart(biomart="ensembl")
listDatasets(mart)
mart <- useMart(host='apr2013.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL", dataset="btaurus_gene_ensembl")
#mart <- useMart(biomart="ensembl", dataset="btaurus_gene_ensembl")

# Check the arguments required for biomaRt query
listFilters(mart)
listAttributes(mart)
head(listAttributes(mart), n=500)
attributePages(mart)

# Build a biomarRt query to obtain ID, name and description for each gene
bioMart_query <- getBM(attributes=c("ensembl_gene_id", "external_gene_id", "description"), filters="ensembl_gene_id", values=rownames(Count$counts), mart=mart, uniqueRows=TRUE, bmHeader=TRUE)
bioMart_query_entrezgene <- getBM(attributes=c("ensembl_gene_id", "entrezgene"), filters="ensembl_gene_id", values=rownames(Count$counts), mart=mart, uniqueRows=TRUE, bmHeader=TRUE)

# Format the biomart_query so that there is only one row per ensembl gene ID
bioMart_entrezgene <- aggregate(entrezgene ~ ensembl_gene_id, FUN = "as.vector", data=bioMart_query_entrezgene, na.action="as.vector")
head(bioMart_entrezgene, n=100)

# Merge the different bioMart queries and custom annotation
bioMart_annotation <- merge(x=bioMart_query, y=sense_info, by="ensembl_gene_id", all=TRUE)
bioMart_annotation <- merge(x=bioMart_annotation, y=bioMart_entrezgene, by="ensembl_gene_id", all=TRUE)
head(bioMart_annotation)
dim(bioMart_annotation)

# Check that the biomaRt output contains all requested genes
length(rownames(Count$counts))
dim(bioMart_annotation)
head(bioMart_annotation)
table(duplicated(bioMart_annotation[,1]))
length(unique(bioMart_annotation[,1], incomparables=FALSE))

# Merge the biomart_query with the count table
Annotated_count <- merge(x=bioMart_annotation, y=Count$counts, by.x="ensembl_gene_id", by.y=0, all=TRUE)
head(Annotated_count)

# Check the merged matrix in terms of size and content
dim(Annotated_count)
table(duplicated(Annotated_count[,1]))
length(unique(Annotated_count[,1], incomparables=FALSE))

#############################
# Filter out all rRNA genes #
#############################

# Read in list of all rRNA genes
list_rRNA <- scan(file="rRNA_gene.txt", what="character", sep="\n")

# Create filter pattern by matching gene ID in the query matrix with the vector list of gene ID which passed low expression filtering
rRNA_to_keep <- Annotated_count$ensembl_gene_id %in% list_rRNA
rRNA_filter <- !(Annotated_count$ensembl_gene_id %in% list_rRNA)
summary(rRNA_to_keep)
summary(rRNA_filter)

# Filter out rRNA genes using pattern
rRNA <- DGEList(counts=Annotated_count[rRNA_to_keep,(ncol(bioMart_annotation)+1):ncol(Annotated_count)], lib.size=NULL, norm.factors=NULL, genes=Annotated_count[rRNA_to_keep,1:ncol(bioMart_annotation)], remove.zeros=FALSE)
rownames(rRNA$counts) <- rownames(rRNA$genes) <- rRNA$genes$ensembl_gene_id
rRNA$genes$ensembl_gene_id <- NULL
Alv_no_rRNA <- Annotated_count[rRNA_filter,]

# Output txt file of rRNA genes raw count
write.table(x=rRNA$counts, file="Alv_rRNA_rawcount.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

##########################################
# Create a DGElist of all samples counts #
##########################################

# Create a target matrix and a experimental group vector and animal block
target <- read.delim(file="target.txt", header=TRUE,sep="\t")
group <- factor(paste(target$Treatment, target$Time_point, sep="."))
animal <- factor(target$Animal)

# Create a DGElist containing the group information
Alv_dgelist <- DGEList(counts=Alv_no_rRNA[,(ncol(bioMart_annotation)+1):ncol(Annotated_count)], lib.size=NULL, norm.factors=NULL, group=group, genes=Alv_no_rRNA[,1:ncol(bioMart_annotation)], remove.zeros=FALSE)
rownames(Alv_dgelist$counts) <- rownames(Alv_dgelist$genes) <- Alv_dgelist$genes$ensembl_gene_id
Alv_dgelist$genes$ensembl_gene_id <- NULL
names(Alv_dgelist)
head(Alv_dgelist$samples)
head(Alv_dgelist$counts)
head(Alv_dgelist$genes)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Log10 transform the count data for better visualization
count_log10 <- log10(x=(Alv_dgelist$counts[,1:ncol(Alv_dgelist$counts)]+1))

# Plot density of count for all libraries
png(filename="Density_Alv_mac.png", width=1366, height=768, units="px")
plot(x=density(count_log10[,1]), main="Density plot of count per gene", lty=1, xlab="Log10 of count per gene", ylab="Density", col="black", ylim=c(0.0,1.25))
for (i in 2:ncol(count_log10)) {
  lines(density(count_log10[,i]), lty=1, col="black")
}
dev.off()

#####################################
# Filtering of lowly expressed tags #
#####################################

# Filter non expressed tags (all samples have zero counts)
Alv_no_zeros <- Alv_dgelist[rowSums(Alv_dgelist$counts) > 0,]
dim(Alv_no_zeros$counts)
head(Alv_no_zeros$samples)

# Filter lowly expressed tags, retaining only tags with at least 1 count per million in 10 or more libraries (10 libraries correspond to one group of treatment)
Alv_filt <- Alv_no_zeros[rowSums(cpm(Alv_no_zeros)>1) >=10,]
dim(Alv_filt$counts)
head(Alv_filt$samples)

# Output txt file of raw count
write.table(x=Alv_filt$counts, file="Alv_filt_rawcount.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Recompute the library size
Alv_filt$samples$lib.size <- colSums(Alv_filt$counts)
head(Alv_filt$samples)
head(Alv_dgelist$samples)

###################################################################################################################
# Normalization of data using trimmed mean of M-values (normalization based on RNA composition between libraries) #
###################################################################################################################

# Calculate normalisation factor for our DGElist (with edgeR, counts are not transformed in any way after normalization, instead normalization will modify library size)
Alv_norm <- calcNormFactors(Alv_filt)
head(Alv_norm$samples)

####################################################################
# Quality check of filtered libraries by plotting density of count #
####################################################################

# Log10 transform the filtered count data for better visualization
count_filt_log10 <- log10(Alv_norm$counts[,1:ncol(Alv_norm$counts)]+1)

# Plot density of count for all libraries
png(filename="Density_Alv_mac_post_filter.png", width=1366, height=768, units="px")
plot(density(count_filt_log10[,1]), main="Density plot of count per gene", lty=1, xlab="Log10 of count per gene", ylab="Density", col="black", ylim=c(0.0,0.6))
for (i in 2:ncol(count_filt_log10)) {
  lines(density(count_filt_log10[,i]), lty=1, col="black")
}
dev.off()

################################################
# Multidimensional scaling plot on all samples #
################################################

# Code below will output the MDS plot in directory as a png file
png(filename="MDS_Alv_mac.png", width=1366, height=768, units="px")
MDS <- plotMDS(x=Alv_norm, top=1000000, gene.selection="pairwise", xlab="Dimension 1", ylab="Dimension 2", col=, cex=)
#points(x=MDS$cmdscale.out[,1], y=MDS$cmdscale.out[,2], type = "p", pch=c(15,17), col=rep(c(4, 2), times=6))
dev.off()

# More visually friendly version of the MDS plot
symbol <- vector()
color <- vector()
for (i in 1:ncol(Alv_norm$counts)) {
  if (length(grep(pattern="CN", x=colnames(Alv_norm$counts)[i]))==1) {
    symbol <- c(symbol, "o")
  }
  if (length(grep(pattern="MB", x=colnames(Alv_norm$counts)[i]))==1) {
    symbol <- c(symbol, "+")
  }
  if (length(grep(pattern="0H", x=colnames(Alv_norm$counts)[i]))==1) {
    color <- c(color, "black")
  }
  if (length(grep(pattern="2H", x=colnames(Alv_norm$counts)[i]))==1) {
    color <- c(color, "gold")
  }
  if (length(grep(pattern="6H", x=colnames(Alv_norm$counts)[i]))==1) {
    color <- c(color, "green")
  }
  if (length(grep(pattern="24H", x=colnames(Alv_norm$counts)[i]))==1) {
    color <- c(color, "blue")
  }
  if (length(grep(pattern="48H", x=colnames(Alv_norm$counts)[i]))==1) {
    color <- c(color, "red")
  }
}
png(filename="MDS_Alv_mac(D1_vs_D3).png", width=1600, height=1000, units="px")
MDS <- plotMDS(x=Alv_norm, top=1000000, gene.selection="pairwise", xlab="Dimension 1", ylab="Dimension 3", col=color, labels=symbol, dim.plot=c(1,3), cex=2.0)
legend(x="bottomright", legend=c("CN", "MB", "2H", "6H", "24H", "48H"), col=c("black", "black", "gold", "green", "blue", "red"), pch=c(1,3,15,15,15,15), cex=2.0)
dev.off()

# Write into a table the coordinates of each library for the MDS plot
write.table(x=MDS$cmdscale.out, file="MDS_xy_all_samples.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

##############################################
# Create a design matrix for paired analysis #
##############################################

# Create a design matrix
design <- model.matrix(~animal+group)
rownames(design) <- rownames(Alv_norm$samples)
colnames(design) <- gsub(pattern="(animal)|(group)", replacement="", x=colnames(design), perl=TRUE)
design

# Write into a table the design matrix
write.table(x=design, file="Alv_design.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

################################################################################################
# Estimate the dispersion parameter for each tag using Cox-Reid method (for multi-factor data) #
################################################################################################

Alv_disp <- estimateGLMCommonDisp(y=Alv_norm, design=design, verbose=TRUE)
Alv_disp <- estimateGLMTrendedDisp(y=Alv_disp, design=design)
Alv_disp <- estimateGLMTagwiseDisp(y=Alv_disp, design=design)
names(Alv_disp)

# Plot the dispersion
png(filename="BCV_Alv_mac.png", width=1366, height=768, units="px")
plotBCV(Alv_disp)
dev.off()

# Show the calculated dispersion
Alv_disp$common.dispersion

# And its square root, the coefficient of biological variation
sqrt(Alv_disp$common.dispersion)

# Create a matrix of the tagwise dispersion associated to ensembl gene
Tagwisedisp <- cbind(Alv_disp$genes, Alv_disp$tagwise.dispersion)
head(Tagwisedisp)
dim(Tagwisedisp)

# Write into a table the calculated tagwise dispersion
write.matrix(x=Tagwisedisp, file="Tagwise_dispersion.txt", sep = "\t")

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using the design matrix and calculated dispersion
Alv_fit <- glmFit(y=Alv_disp, design=design)
names(Alv_fit)

###########################################
# Differential expression call MB2H-CN2H #
###########################################

# Carry out the likeli-hood ratio test
Alv_lrt_MB.2H <- glmLRT(glmfit=Alv_fit, contrast=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0, 0))
names(Alv_lrt_MB.2H)
Alv_lrt_MB.2H$comparison
head(Alv_lrt_MB.2H$table)

# Summarise the number of up- and down-regulated genes
summary(decideTestsDGE(Alv_lrt_MB.2H, p.value = 0.05))

# Ajust for multiple testing
DE_MB.2H <- topTags(Alv_lrt_MB.2H, n="inf", adjust.method="BH")
names(DE_MB.2H)
head(DE_MB.2H$table)

# Output data
write.table(x=DE_MB.2H$table[,c(1,2,10:14)], file="DE_MB_sense.2H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Plot the FC versus CPM per tagwise
png(filename="Smear_FC_CPM_MB.2H.png", width=1366, height=768, units="px")
plotSmear(object=Alv_lrt_MB.2H, de.tags=(rownames(Alv_lrt_MB.2H$table)[as.logical(decideTestsDGE(Alv_lrt_MB.2H, p.value = 0.05))]), ylim=c(-1.9, 2.1))
abline(h=c(-1, 1), col="blue")
dev.off()

###########################################
# Differential expression call MB6H-CN6H #
###########################################

# Carry out the likeli-hood ratio test
Alv_lrt_MB.6H <- glmLRT(glmfit=Alv_fit, contrast=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1))
names(Alv_lrt_MB.6H)
Alv_lrt_MB.6H$comparison
head(Alv_lrt_MB.6H$table)

# Summarise the number of up- and down-regulated genes
summary(decideTestsDGE(Alv_lrt_MB.6H, p.value = 0.05))

# Ajust for multiple testing
DE_MB.6H <- topTags(Alv_lrt_MB.6H, n="inf", adjust.method="BH")
names(DE_MB.6H)
head(DE_MB.6H$table)

# Output data
write.table(x=DE_MB.6H$table[,c(1,2,10:14)], file="DE_MB_sense.6H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Plot the FC versus CPM per tagwise
png(filename="Smear_FC_CPM_MB.6H.png", width=1366, height=768, units="px")
plotSmear(object=Alv_lrt_MB.6H, de.tags=(rownames(Alv_lrt_MB.6H$table)[as.logical(decideTestsDGE(Alv_lrt_MB.6H, p.value = 0.05))]), ylim=c(-3.0, 4.0))
abline(h=c(-1, 1), col="blue")
dev.off()

###########################################
# Differential expression call MB24H-CN24H #
###########################################

# Carry out the likeli-hood ratio test
Alv_lrt_MB.24H <- glmLRT(glmfit=Alv_fit, coef=14)
names(Alv_lrt_MB.24H)
Alv_lrt_MB.24H$comparison
head(Alv_lrt_MB.24H$table)

# Summarise the number of up- and down-regulated genes
summary(decideTestsDGE(Alv_lrt_MB.24H, p.value = 0.05))

# Ajust for multiple testing
DE_MB.24H <- topTags(Alv_lrt_MB.24H, n="inf", adjust.method="BH")
names(DE_MB.24H)
head(DE_MB.24H$table)

# Output data
write.table(x=DE_MB.24H$table[,c(1,2,10:14)], file="DE_MB_sense.24H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Plot the FC versus CPM per tagwise
png(filename="Smear_FC_CPM_MB.24H.png", width=1366, height=768, units="px")
plotSmear(object=Alv_lrt_MB.24H, de.tags=(rownames(Alv_lrt_MB.24H$table)[as.logical(decideTestsDGE(Alv_lrt_MB.24H, p.value = 0.05))]), ylim=c(-4.3, 7.1))
abline(h=c(-1, 1), col="blue")
dev.off()

###########################################
# Differential expression call MB48H-CN48H #
###########################################

# Carry out the likeli-hood ratio test
Alv_lrt_MB.48H <- glmLRT(glmfit=Alv_fit, contrast=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 1, 0))
names(Alv_lrt_MB.48H)
Alv_lrt_MB.48H$comparison
head(Alv_lrt_MB.48H$table)

# Summarise the number of up- and down-regulated genes
summary(decideTestsDGE(Alv_lrt_MB.48H, p.value = 0.05))

# Ajust for multiple testing
DE_MB.48H <- topTags(Alv_lrt_MB.48H, n="inf", adjust.method="BH")
names(DE_MB.48H)
head(DE_MB.48H$table)

# Output data
write.table(x=DE_MB.48H$table[,c(1,2,10:14)], file="DE_MB_sense.48H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Plot the FC versus CPM per tagwise
png(filename="Smear_FC_CPM_MB.48H.png", width=1366, height=768, units="px")
plotSmear(object=Alv_lrt_MB.48H, de.tags=(rownames(Alv_lrt_MB.48H$table)[as.logical(decideTestsDGE(Alv_lrt_MB.48H, p.value = 0.05))]), ylim=c(-6.0, 8.0))
abline(h=c(-1, 1), col="blue")
dev.off()

#########################################################
# Merge all DE call data from the different time points #
#########################################################

# Merge all DE table for the different time points into a single dataframe
Full_DE <- merge(x=DE_MB.2H$table[,1:ncol(DE_MB.2H$table)], y=DE_MB.6H$table[,(ncol(DE_MB.6H$table)-4):ncol(DE_MB.6H$table)], by="row.names")
rownames(Full_DE) <- Full_DE[,1]
Full_DE <- Full_DE[,-1]
colnames(Full_DE) <- gsub(pattern=".x$", replacement="_2H", x=colnames(Full_DE), perl=TRUE)
colnames(Full_DE) <- gsub(pattern=".y$", replacement="_6H", x=colnames(Full_DE), perl=TRUE)
Full_DE <- merge(x=Full_DE, y=DE_MB.24H$table[,(ncol(DE_MB.24H$table)-4):ncol(DE_MB.24H$table)], by="row.names")
rownames(Full_DE) <- Full_DE[,1]
Full_DE <- Full_DE[,-1]
Full_DE <- merge(x=Full_DE, y=DE_MB.48H$table[,(ncol(DE_MB.48H$table)-4):ncol(DE_MB.48H$table)], by="row.names")
rownames(Full_DE) <- Full_DE[,1]
colnames(Full_DE) <- gsub(pattern=".x$", replacement="_24H", x=colnames(Full_DE), perl=TRUE)
colnames(Full_DE) <- gsub(pattern=".y$", replacement="_48H", x=colnames(Full_DE), perl=TRUE)
colnames(Full_DE)[1] <- "ensembl_gene_id"
head(Full_DE)

# Write into a table the full DE call data
write.matrix(x=Full_DE, file="Full_DE_sense.txt", sep = "\t")

########################
# Time series analysis #
########################

# Select biomaRt database and dataset (use biomaRt from Ensembl 71: Apr 2013)
listMarts(host='Apr2013.archive.ensembl.org')
GO_mart <- useMart(host='Apr2013.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL")
listDatasets(GO_mart)
GO_mart <- useMart(host='Apr2013.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL", dataset = "btaurus_gene_ensembl")

# Build a biomarRt query to obtain ID, name and description for each GO categories associated with genes
GO_bioMart_query <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006", "definition_1006"), filters="ensembl_gene_id", values=rownames(Count$counts), mart=GO_mart, uniqueRows=TRUE, bmHeader=TRUE)
head(GO_bioMart_query)
GO_description <- GO_bioMart_query[,-1]
dim(GO_description)
GO_description <- unique(GO_description)
dim(GO_description)
head(GO_description)

# Format the GO annotation for compatibility with STEM
GO_annotation <- GO_bioMart_query[,c(1,2)]
GO_annotation <- aggregate(go_id ~ ensembl_gene_id, FUN = "as.vector", data=GO_annotation, na.action="as.vector")
head(GO_annotation)

# Ouptut GO annotation and description data
write.matrix(x=GO_annotation, file="GO_annotation.txt", sep = "\t")
write.table(x=GO_description, file="GO_description.txt", sep = "\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

# Merge all DE table for the different time points into a single dataframe
time_serie <- merge(x=DE_MB.2H$table[,c(3,7)], y=DE_MB.6H$table[,c(3,7)], by="row.names")
rownames(time_serie) <- time_serie[,1]
time_serie <- time_serie[,-1]
colnames(time_serie) <- c("logFC_2H", "FDR_2H", "logFC_6H", "FDR_6H")
time_serie <- merge(x=time_serie, y=DE_MB.24H$table[,c(3,7)], by="row.names")
rownames(time_serie) <- time_serie[,1]
time_serie <- time_serie[,-1]
colnames(time_serie) <- c("logFC_2H", "FDR_2H", "logFC_6H", "FDR_6H", "logFC_24H", "FDR_24H")
time_serie <- merge(x=time_serie, y=DE_MB.48H$table[,c(3,7)], by="row.names")
rownames(time_serie) <- time_serie[,1]
time_serie <- time_serie[,-1]
colnames(time_serie) <- c("logFC_2H", "FDR_2H", "logFC_6H", "FDR_6H", "logFC_24H", "FDR_24H", "logFC_48H", "FDR_48H")
head(time_serie)

# Transform dataframe into matrix and add one new column
time_serie <- data.matrix(time_serie)
class(time_serie)

# Filter out genes which are not differentially expressed at any of the time point
list_DEG <- scan(file="list_DEG.txt", what="character", sep="\n")
DE_time_serie <- rownames(time_serie) %in% list_DEG
time_serie <- time_serie[DE_time_serie,]
head(time_serie)

# Remove the FDR value from matrix in order to be compatible with STEM software
STEM_time_serie <- time_serie[,c(1,3,5,7)]
head(STEM_time_serie)

# Ouptut DEG data over time
write.table(x=STEM_time_serie, file="STEM_time_serie.txt",  sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=time_serie, file="DEG_time_serie.txt",  sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#####################################################
# Bovine gene annotation by orthology to H. sapiens #
#####################################################

# Select the full list of bovine genes ID
Bovine_to_human <- Annotated_count[,1:2]
head(Bovine_to_human)

# Select biomaRt database and dataset (use biomaRt from Ensembl 71: Apr 2013)
mart <- useMart(host='apr2013.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL", dataset="btaurus_gene_ensembl")

# Build a biomaRt query to obtain ortholog H. sapiens Ensembl ID for each bovine gene
bioMart_query_ortholog <- getBM(attributes=c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene"), filters="ensembl_gene_id", values=Bovine_to_human, mart=mart, uniqueRows=TRUE, bmHeader=TRUE)
head(bioMart_query_ortholog, n=10)
dim(bioMart_query_ortholog)

# Build a biomaRt query to obtain H. sapiens gene ID and gene symbol
mart_hsapiens <- useMart(host='apr2013.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
bioMart_query_hsapiens <- getBM(attributes=c("ensembl_gene_id", "external_gene_id"), mart=mart_hsapiens, uniqueRows=TRUE, bmHeader=TRUE)
colnames(bioMart_query_hsapiens) <- c("hsapiens_homolog_ensembl_gene", "hsapiens_homolog_external_gene_id")
head(bioMart_query_hsapiens, n=10)
dim(bioMart_query_hsapiens)

# Merge the different biomaRt queries into a single matrix
ortholog <- merge(x=bioMart_query_ortholog, y=bioMart_query_hsapiens, by="hsapiens_homolog_ensembl_gene", all=TRUE)
ortholog <- merge(x=Annotated_count[,1:2], y=ortholog, by="ensembl_gene_id", all.x=TRUE)
head(ortholog, n=500)
tail(ortholog)

# Use a loop and hashes to remove all duplicate bovine gene ID and so obtain a unique human ortholog per bovine gene
ortholog_value <- new.env()
ortholog_id <- new.env()
for (i in 1:nrow(ortholog)) {
  if ((length(ortholog[i,3]) == 0) || is.na(ortholog[i,3])) {
    score <- 0
  }
  else if ((length(ortholog[i,2]) == 0) || (is.na(ortholog[i,2])) || (length(ortholog[i,4]) == 0) || (is.na(ortholog[i,4])) || (ortholog[i,2] != ortholog[i,4])) {
    score <- 1
  }
  else if (ortholog[i,2] == ortholog[i,4]) {
    score <- 2
  }
  else {
    stop("Error!")
  }
  if (exists(x=ortholog[i,1], envir=ortholog_value, inherits = FALSE)) {
    if (score > ortholog_value[[ ortholog[i,1] ]]) {
      ortholog_value[[ ortholog[i,1] ]] <- score
      ortholog_id[[ ortholog[i,1] ]] <- ortholog[i,3]
    }
    else if (score == ortholog_value[[ ortholog[i,1] ]]) {
      ortholog_id[[ ortholog[i,1] ]] <- ""
    }
  }
  else {
    if (score == 0) {
      ortholog_value[[ ortholog[i,1] ]] <- score
    }
    else {
      ortholog_value[[ ortholog[i,1] ]] <- score
      ortholog_id[[ ortholog[i,1] ]] <- ortholog[i,3]
    }
  }
}

# Create the matrix containing the unique human ortholog per bovine gene
final_ortholog <- matrix(nrow=nrow(Annotated_count), ncol=2)
for (i in 1:nrow(Annotated_count)) {
  final_ortholog[i,1] <- Annotated_count[i,1]
  final_ortholog[i,2] <- ortholog_id[[ Annotated_count[i,1] ]]
}
colnames(final_ortholog) <- c("ensembl_gene_id", "hsapiens_homolog_ensembl_gene")
head(final_ortholog)
summary(duplicated(final_ortholog[,1]))
summary(duplicated(final_ortholog[,2], incomparables=c("", "NA")))

# Add the human and bovine gene symbol to the matrix
final_ortholog <- merge(x=final_ortholog, y=Annotated_count[,1:2], by="ensembl_gene_id")
final_ortholog <- merge(x=final_ortholog, y=bioMart_query_hsapiens, by="hsapiens_homolog_ensembl_gene")
final_ortholog <- final_ortholog[,c(1,4,2,3)]
final_ortholog <- as.matrix(final_ortholog)
final_ortholog[sort.list(final_ortholog[,1]), ]
head(final_ortholog, n=500)
tail(final_ortholog)

# Use a loop and hashes to remove all duplicate human gene ID and so obtain a unique bovine gene per human ortholog
human_value <- new.env()
human_id <- new.env()
for (i in 1:nrow(final_ortholog)) {
  if ((length(final_ortholog[i,3]) == 0) || is.na(final_ortholog[i,3])) {
    score <- 0
  }
  else if ((length(final_ortholog[i,2]) == 0) || (is.na(final_ortholog[i,2])) || (length(final_ortholog[i,4]) == 0) || (is.na(final_ortholog[i,4])) || (final_ortholog[i,2] != final_ortholog[i,4])) {
    score <- 1
  }
  else if (final_ortholog[i,2] == final_ortholog[i,4]) {
    score <- 2
  }
  else {
    stop("Error!")
  }
  if (exists(x=final_ortholog[i,1], envir=human_value, inherits = FALSE)) {
    if (score > human_value[[ final_ortholog[i,1] ]]) {
      human_value[[ final_ortholog[i,1] ]] <- score
      human_id[[ final_ortholog[i,1] ]] <- final_ortholog[i,3]
    }
    else if (score == human_value[[ final_ortholog[i,1] ]]) {
      human_id[[ final_ortholog[i,1] ]] <- ""
    }
  }
  else {
    if (score == 0) {
      human_value[[ final_ortholog[i,1] ]] <- score
    }
    else {
      human_value[[ final_ortholog[i,1] ]] <- score
      human_id[[ final_ortholog[i,1] ]] <- final_ortholog[i,3]
    }
  }
}

# Create the matrix containing the bovine-human gene orthology (with no duplication)
human_ortholog <- matrix(nrow=length(human_id), ncol=2)
human_ortholog[,1] <- ls(human_id)
for (i in 1:length(human_id)) {
  human_ortholog[i,2] <- human_id[[ human_ortholog[i,1] ]]
}
colnames(human_ortholog) <- c("hsapiens_homolog_ensembl_gene", "ensembl_gene_id")
head(human_ortholog, n=10)
summary(duplicated(human_ortholog[,1]))
summary(duplicated(human_ortholog[,2], incomparables=c("", "NA")))

#############################################################
# Annotate all DE genes by orthology and add all novel gene #
#############################################################

# Merge the H. sapiens ortholog gene ID annotation with DE results
ortholog_full_DE <- merge(x=human_ortholog, y=Full_DE[,c(1,11:30)], by="ensembl_gene_id", all.y=TRUE)
colnames(ortholog_full_DE)[2] <- "Hsapiens_ensembl_gene_id"
ortholog_full_DE <- ortholog_full_DE[(!is.na(x=ortholog_full_DE$Hsapiens_ensembl_gene_id)),-1]
head(ortholog_full_DE)
dim(ortholog_full_DE)

# Read in the novel DE genes
novel_genes <- read.table(file="Full_DE_novel.txt", sep="\t", header=TRUE, fill=TRUE)
novel_genes <- novel_genes[,c(1,2,11:30)]
novel_genes <- novel_genes[(!is.na(x=novel_genes$Hsapiens_ensembl_gene_id)),-1]
head(novel_genes)
dim(novel_genes)

# Merge the H. sapiens ortholog with the novel gene
ortholog_novel_DE <- rbind(a=ortholog_full_DE, b=novel_genes)
head(ortholog_novel_DE)
dim(ortholog_novel_DE)
summary(duplicated(ortholog_novel_DE[,1]))

# Split data per time points
ortholog_novel_DE_MB_2H <- ortholog_novel_DE[,1:6]
colnames(ortholog_novel_DE_MB_2H) <- c("Hsapiens_ensembl_gene_id", "logFC", "logCPM", "LR", "PValue", "FDR")
head(ortholog_novel_DE_MB_2H)
ortholog_novel_DE_MB_6H <- ortholog_novel_DE[,c(1,7:11)]
colnames(ortholog_novel_DE_MB_6H) <- c("Hsapiens_ensembl_gene_id", "logFC", "logCPM", "LR", "PValue", "FDR")
head(ortholog_novel_DE_MB_6H)
ortholog_novel_DE_MB_24H <- ortholog_novel_DE[,c(1,12:16)]
colnames(ortholog_novel_DE_MB_24H) <- c("Hsapiens_ensembl_gene_id", "logFC", "logCPM", "LR", "PValue", "FDR")
head(ortholog_novel_DE_MB_24H)
ortholog_novel_DE_MB_48H <- ortholog_novel_DE[,c(1,17:21)]
colnames(ortholog_novel_DE_MB_48H) <- c("Hsapiens_ensembl_gene_id", "logFC", "logCPM", "LR", "PValue", "FDR")
head(ortholog_novel_DE_MB_48H)

#############################################################################
# Use of Sigora for over-representation of gene-pair signatures in pathways #
#############################################################################

# Create a function to perform Sigora analysis based on user-defined parameters
sigora_analysis <- function(input, fdr, logfc, direction, output) {
  # Load the library Sigora
  library(sigora)
  # Output parameters of the analysis
  write.table(x=paste("Sigora analysis on input: ", deparse(substitute(input)), " with gene below FDR: ", fdr, " and direction of expression: ", direction, " of ", logfc, "log fold-change."), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Create the list of target gene according to user specified parameters
  if (direction == "up") {
    target <- input$table[(input$table$FDR < fdr),]
    target <- rownames(target[(target$logFC > logfc),])
  }
  else if (direction == "down") {
    target <- input$table[(input$table$FDR < fdr),]
    target <- rownames(target[(target$logFC < logfc),])
  }
  else if (direction == "all") {
    target <- input$table[(input$table$FDR < fdr),]
    target <- rownames(target[((target$logFC > logfc) || (target$logFC < logfc)),])
  }
  target <- as.matrix(target)
  colnames(target) <- "ensembl_gene_id"
  write.table(x=paste("Number of bovine target genes: ", length(target)), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Obtain the human gene ID for the target bovine gene
  homo_id <- merge(x=target, y=human_ortholog, by="ensembl_gene_id")
  homo_id <- as.vector(homo_id[,-1])
  write.table(x=paste("Number of human ortholog target genes: ", length(homo_id)), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Convert human gene ID to Sigora ID for analysis purpose
  sig_id <<- ens_converter(homo_id)
  write.table(x=paste("Number of Sigora corresponding ID: ", length(sig_id)), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x="------------------------------------------------------------------------------------------------------------------------------------------------------", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Perform Sigora pathway analysis using Reactome database and outpout results
  sink(file=output, append=TRUE, split=TRUE)
  sigs(samplename=sig_id, archive="R", markers=1, level=4)
  sink()
  write.table(x="Summary results (based on Reactome):", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x=summary_results, file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x="------------------------------------------------------------------------------------------------------------------------------------------------------", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Perform Sigora pathway analysis using Kegg database and outpout results
  sink(file=output, append=TRUE, split=TRUE)
  sigs(samplename=sig_id, archive="k", markers=1, level=2)
  sink()
  write.table(x="Summary results (based on Kegg):", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x=summary_results, file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
}

# Perform Sigora analysis using different parameters
sigora <- sigora_analysis(input=DE_MB.2H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB2H.txt")
sigora <- sigora_analysis(input=DE_MB.2H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB2H.txt")
sigora <- sigora_analysis(input=DE_MB.2H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB2H.txt")
sigora <- sigora_analysis(input=DE_MB.2H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB2H.txt")
sigora <- sigora_analysis(input=DE_MB.6H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB6H.txt")
sigora <- sigora_analysis(input=DE_MB.6H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB6H.txt")
sigora <- sigora_analysis(input=DE_MB.6H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB6H.txt")
sigora <- sigora_analysis(input=DE_MB.6H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB6H.txt")
sigora <- sigora_analysis(input=DE_MB.24H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB24H.txt")
sigora <- sigora_analysis(input=DE_MB.24H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB24H.txt")
sigora <- sigora_analysis(input=DE_MB.24H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB24H.txt")
sigora <- sigora_analysis(input=DE_MB.24H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB24H.txt")
sigora <- sigora_analysis(input=DE_MB.48H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB48H.txt")
sigora <- sigora_analysis(input=DE_MB.48H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB48H.txt")
sigora <- sigora_analysis(input=DE_MB.48H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB48H.txt")
sigora <- sigora_analysis(input=DE_MB.48H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB48H.txt")

# Create a function to perform Sigora analysis based on user-defined parameters
sigora_analysis_human <- function(input, fdr, logfc, direction, output) {
  # Load the library Sigora
  library(sigora)
  # Output parameters of the analysis
  write.table(x=paste("Sigora analysis on input: ", deparse(substitute(input)), " with gene below FDR: ", fdr, " and direction of expression: ", direction, " of ", logfc, "log fold-change."), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Create the list of target gene according to user specified parameters
  if (direction == "up") {
    target <- input[(input$FDR < fdr),]
    target <- target[(target$logFC > logfc),1]
  }
  else if (direction == "down") {
    target <- input[(input$FDR < fdr),]
    target <- target[(target$logFC < logfc),1]
  }
  else if (direction == "all") {
    target <- input[(input$FDR < fdr),]
    target <- target[((target$logFC > logfc) || (target$logFC < logfc)),1]
  }
  target <- as.vector(target)
  write.table(x=paste("Number of human target genes: ", length(target)), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Convert human gene ID to Sigora ID for analysis purpose
  sig_id <<- ens_converter(target)
  write.table(x=paste("Number of Sigora corresponding ID: ", length(sig_id)), file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x="------------------------------------------------------------------------------------------------------------------------------------------------------", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Perform Sigora pathway analysis using Reactome database and outpout results
  sink(file=output, append=TRUE, split=TRUE)
  sigs(samplename=sig_id, archive="R", markers=1, level=4)
  sink()
  write.table(x="Summary results (based on Reactome):", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x=summary_results, file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x="------------------------------------------------------------------------------------------------------------------------------------------------------", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  # Perform Sigora pathway analysis using Kegg database and outpout results
  sink(file=output, append=TRUE, split=TRUE)
  sigs(samplename=sig_id, archive="k", markers=1, level=2)
  sink()
  write.table(x="Summary results (based on Kegg):", file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
  write.table(x=summary_results, file=output, sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE, append=TRUE)
}

# Perform Sigora analysis using different parameters for the human ortholog and the novel genes
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_2H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB2H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_2H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB2H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_2H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB2H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_2H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB2H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_6H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB6H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_6H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB6H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_6H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB6H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_6H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB6H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_24H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB24H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_24H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB24H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_24H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB24H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_24H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB24H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_48H, fdr=0.05, logfc=0, direction="up", output="Sigora_Up_MB48H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_48H, fdr=0.05, logfc=0, direction="down", output="Sigora_Down_MB48H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_48H, fdr=0.05, logfc=0, direction="all", output="Sigora_MB48H(novel).txt")
sigora <- sigora_analysis_human(input=ortholog_novel_DE_MB_48H, fdr=0.001, logfc=0, direction="all", output="Sigora_fdr0.001_MB48H(novel).txt")

# Clean up the unrequired variables
remove(sig_id, sigora, summary_results, detailed_results)

########################################################################
# Use of Pathway-guide (SPIA) for identification of disturbed pathways #
########################################################################

# Install required libraries (and make sure to install them in library folder where Pathway-guide is looking)
source("http://bioconductor.org/biocLite.R")
biocLite("ROntoTools")
source("http://bioconductor.org/biocLite.R")
biocLite("KEGGREST")
source("http://bioconductor.org/biocLite.R")
biocLite("KEGGgraph")
source("http://bioconductor.org/biocLite.R")
biocLite("Rgraphviz")
install.packages("httr")
install.packages("stringr")
install.packages("digest")
install.packages("png")
install.packages("")
install.packages("")
install.packages("")
install.packages("")
install.packages("")
install.packages("")
install.packages("")
install.packages("")

# Output the different data to use in Pathway-guide
write.table(x=DE_MB.2H$table[(DE_MB.2H$table$FDR < 0.05),c(10,14)], file="Sig_MB_sense.2H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=DE_MB.6H$table[(DE_MB.6H$table$FDR < 0.05),c(10,14)], file="Sig_MB_sense.6H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=DE_MB.24H$table[(DE_MB.24H$table$FDR < 0.05),c(10,14)], file="Sig_MB_sense.24H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=DE_MB.48H$table[(DE_MB.48H$table$FDR < 0.05),c(10,14)], file="Sig_MB_sense.48H.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=ortholog_novel_DE_MB_2H[(ortholog_novel_DE_MB_2H$FDR < 0.05),c(1,2,5)], file="Sig_MB_sense.2H(novel).txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(x=ortholog_novel_DE_MB_6H[(ortholog_novel_DE_MB_6H$FDR < 0.05),c(1,2,5)], file="Sig_MB_sense.6H(novel).txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(x=ortholog_novel_DE_MB_24H[(ortholog_novel_DE_MB_24H$FDR < 0.05),c(1,2,5)], file="Sig_MB_sense.24H(novel).txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(x=ortholog_novel_DE_MB_48H[(ortholog_novel_DE_MB_48H$FDR < 0.05),c(1,2,5)], file="Sig_MB_sense.48H(novel).txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
write.table(x=ortholog_novel_DE, file="DE_MB_sense_novel(Hsapiens).txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

#########################################################################
# Compare the different pathways obtained with Sigora and Pathway-guide #
#########################################################################

# Load required library
library(gdata)

# Read in the lists of pathways obtained with Pathway-guide from kegg
Path_MB2H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_2H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2, )
colnames(Path_MB2H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB2H_novel), perl=TRUE)
head(Path_MB2H_novel)
dim(Path_MB2H_novel)
Path_MB6H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_6H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB6H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB6H_novel), perl=TRUE)
head(Path_MB6H_novel)
dim(Path_MB6H_novel)
Path_MB24H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_24H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB24H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB24H_novel), perl=TRUE)
head(Path_MB24H_novel)
dim(Path_MB24H_novel)
Path_MB48H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_48H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB48H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB48H_novel), perl=TRUE)
head(Path_MB48H_novel)
dim(Path_MB48H_novel)

# Read in the lists of pathways obtained with Pathway-guide from Reactome
Path_MB2H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_2H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2, )
colnames(Path_MB2H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB2H_novel), perl=TRUE)
head(Path_MB2H_novel)
dim(Path_MB2H_novel)
Path_MB6H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_6H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB6H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB6H_novel), perl=TRUE)
head(Path_MB6H_novel)
dim(Path_MB6H_novel)
Path_MB24H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_24H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB24H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB24H_novel), perl=TRUE)
head(Path_MB24H_novel)
dim(Path_MB24H_novel)
Path_MB48H_novel <- read.xls(xls="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_48H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB48H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB48H_novel), perl=TRUE)
head(Path_MB48H_novel)
dim(Path_MB48H_novel)

# Create a function to read in the lists of Kegg pathways obtained with Sigora
sigora_input <- function(file) {
  sigora_path <- read.table(file=file, header=FALSE, sep="\t", quote="", dec=".", na.strings="NA", col.names=c("row", "count",  "weightsum",  "pwyid",  "pwyname",  "p_value",  "adj.p_value",	"round_sum",	"m",	"n",	"k"), fill=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE)
  sigora_path <- sigora_path[(grep(pattern="TRUE", x=(sigora_path[,2]=="Summary results (based on Kegg):")))+2:nrow(sigora_path),-1]
  sigora_path <- sigora_path[complete.cases(sigora_path[,1]),]
  row.names(sigora_path) <- sigora_path$pwyname
  sigora_path <- sigora_path[,-4]
  sigora_path$p_value <- as.numeric(as.character(sigora_path$p_value))
  sigora_path$adj.p_value <- as.numeric(as.character(sigora_path$adj.p_value))
  sigora_path <- sigora_path[order(sigora_path$p_value, na.last=TRUE),]
  sigora_path <- cbind(x=sigora_path, y=rank(sigora_path$adj.p_value, na.last="keep", ties.method="first"))
  colnames(sigora_path)[10] <- "rank"
  colnames(sigora_path) <- gsub(pattern="$", replacement="_sigora", x=colnames(sigora_path), perl=TRUE)
  colnames(sigora_path) <- gsub(pattern="^x.", replacement="", x=colnames(sigora_path), perl=TRUE)
  return(sigora_path)
}

# Create a function to read in the lists of Reactome pathways obtained with Sigora
sigora_input <- function(file) {
  sigora_path <- read.table(file=file, header=FALSE, sep="\t", quote="", dec=".", na.strings="NA", col.names=c("row", "count",  "weightsum",  "pwyid",  "pwyname",  "p_value",  "adj.p_value",  "round_sum",	"m",	"n",	"k"), fill=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE)
  sigora_path <- sigora_path[((grep(pattern="TRUE", x=(sigora_path[,2]=="Summary results (based on Reactome):")))+2):((grep(pattern="TRUE", x=(sigora_path[,1]=='[1] "KEGG"')))-4),-1]
  sigora_path <- sigora_path[complete.cases(sigora_path[,1]),]
  row.names(sigora_path) <- sigora_path$pwyname
  sigora_path <- sigora_path[,-4]
  sigora_path$p_value <- as.numeric(as.character(sigora_path$p_value))
  sigora_path$adj.p_value <- as.numeric(as.character(sigora_path$adj.p_value))
  sigora_path <- sigora_path[order(sigora_path$p_value, na.last=TRUE),]
  sigora_path <- cbind(x=sigora_path, y=rank(sigora_path$adj.p_value, na.last="keep", ties.method="first"))
  colnames(sigora_path)[10] <- "rank"
  colnames(sigora_path) <- gsub(pattern="$", replacement="_sigora", x=colnames(sigora_path), perl=TRUE)
  colnames(sigora_path) <- gsub(pattern="^x.", replacement="", x=colnames(sigora_path), perl=TRUE)
  return(sigora_path)
}

# Read in the lists of pathways obtained with Sigora
Sig_MB2H_novel <- sigora_input(file="./Sigora/Sigora_MB2H(novel).txt")
head(Sig_MB2H_novel)
dim(Sig_MB2H_novel)
Sig_MB6H_novel <- sigora_input(file="./Sigora/Sigora_MB6H(novel).txt")
head(Sig_MB6H_novel)
dim(Sig_MB6H_novel)
Sig_MB24H_novel <- sigora_input(file="./Sigora/Sigora_MB24H(novel).txt")
head(Sig_MB24H_novel)
dim(Sig_MB24H_novel)
Sig_MB48H_novel <- sigora_input(file="./Sigora/Sigora_MB48H(novel).txt")
head(Sig_MB48H_novel)
dim(Sig_MB48H_novel)

# Create a function to compare the pathways identified with Pathway-guide and Sigora
common_path <- function(file1, file2) {
  common_path <- as.data.frame(merge(x=file1, y=file2, by="row.names"))
  rownames(common_path) <- common_path$Row.names
  colnames(common_path)[1] <- "significant"
  common_path <- common_path[,c(2:ncol(common_path),1)]
  for (i in 1:nrow(common_path)) {
    if (is.na(common_path$pG.FDR_Pathguide[i]) || is.na(common_path$adj.p_value_sigora)) {
      common_path$significant[i] <- "FALSE"
    }
    else if (common_path$pG.FDR_Pathguide[i]<0.05 && common_path$adj.p_value_sigora[i]<0.05) {
      common_path$significant[i] <- "TRUE"
    }
    else {
      common_path$significant[i] <- "FALSE"
    }
  }
  common_path <- common_path[order(common_path$significant, decreasing=TRUE, na.last=TRUE),]
  return(common_path)
}

# Identify the common significant pathways between Pathway-guide and Sigora
Pathway_MB2H <- common_path(file1=Path_MB2H_novel, file2=Sig_MB2H_novel)
head(Pathway_MB2H)
summary(Pathway_MB2H$significant=="TRUE")
Pathway_MB6H <- common_path(file1=Path_MB6H_novel, file2=Sig_MB6H_novel)
head(Pathway_MB6H)
summary(Pathway_MB6H$significant=="TRUE")
Pathway_MB24H <- common_path(file1=Path_MB24H_novel, file2=Sig_MB24H_novel)
head(Pathway_MB24H)
summary(Pathway_MB24H$significant=="TRUE")
Pathway_MB48H <- common_path(file1=Path_MB48H_novel, file2=Sig_MB48H_novel)
head(Pathway_MB48H)
summary(Pathway_MB48H$significant=="TRUE")

# Output the kegg common pathways
write.table(x=Pathway_MB2H, file="Pathway_MB2H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB6H, file="Pathway_MB6H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB24H, file="Pathway_MB24H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB48H, file="Pathway_MB48H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

# Output the reactome common pathways
write.table(x=Pathway_MB2H, file="Pathway_MB2H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB6H, file="Pathway_MB6H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB24H, file="Pathway_MB24H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB48H, file="Pathway_MB48H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#######
# END #
#######