############################################################
# RNA-seq data analysis of novel gene counts (paired data) #
############################################################

# Analysis nearly identical to previous one, with the exception that for article purpose all CN0H samples have been removed as well as all TB samples

#############################
# List of required packages #
#############################

# Create a function to load or install (then load) the required packages
loadpackage <- function(package) {
  if (require(package=deparse(substitute(package)), character.only=TRUE, quietly=TRUE)) {
    print(paste(deparse(substitute(package)), " is loaded correctly!", sep=""))
  }
  else {
    print(paste("Trying to install ", deparse(substitute(package)), sep=""))
    install.packages(pkgs=deparse(substitute(package)), quiet=TRUE)
    if(require(package=deparse(substitute(package)), character.only=TRUE, quietly=TRUE)) {
      print(paste(deparse(substitute(package)), " is correctly installed and loaded from CRAN!", sep=""))
    }
    else {
      source(file="http://bioconductor.org/biocLite.R", verbose=FALSE)
      biocLite(pkgs=deparse(substitute(package)), suppressUpdates=TRUE)
      if(require(package=deparse(substitute(package)), character.only=TRUE, quietly=TRUE)) {
        print(paste(deparse(substitute(package)), " is correctly installed and loaded from Bioconductor!", sep=""))
      }
      else {
        stop(paste('"', "Could not install ", deparse(substitute(package)), '"', sep=""))
      }
    }
  }
  print(paste(deparse(substitute(package)), " version: ", packageVersion(pkg=deparse(substitute(package))), sep=""))
}

# Load the required packages
loadpackage(package=edgeR)
loadpackage(package=biomaRt)
loadpackage(package=MASS)

######################################################
# Use featureCounts output files as input files in R #
######################################################

###############
# Preparation #
###############

# Move to the appropriate folder
setwd("F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Novel_gene")
getwd()
workDir <- getwd()
workDir

################################################
# Read in and concatenate input files within R #
################################################

# Create vector of all files name
fileDir <- "F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Novel_gene/Counts"
files <- list.files(path=fileDir, pattern="*_(2|6|24|48)H$", all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE)
files

# Reads and merges a set of files containing counts
Count <- readDGE(files=files, path=fileDir, columns=c(1,3))
names(Count)
head(Count$samples)
head(Count$counts)

# Ouptut samples data
write.table(x=Count$samples, file="Alv_samples.txt",  sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Count$counts, file="Alv_rawcounts.txt",  sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

########################################################################################
# BiomaRt annotation of novel gene using human ortholog and using Cufflinks annotation #
########################################################################################

# Read in custom annotation information for all novel genes
novel_info <- read.table(file="gene_info_251113.txt", sep="\t", header=TRUE, fill=TRUE)RBH_annotation
novel_info <- novel_info[grep(pattern="^NOVBTAG", x=novel_info[,1], perl=TRUE),]
colnames(novel_info)[1] <- "novel_gene_id"
colnames(novel_info)[7] <- "novel_transcript_id"
dim(novel_info)
head(novel_info)

# Read the input RBH file
RBH <- read.table(file="RBH_181113_noduplicate.txt", header=TRUE, colClasses=c(rep(x="character", 4), "numeric", "numeric", "numeric"), sep="\t")
table(duplicated(RBH[,1]))
table(duplicated(RBH[,3]))
head(RBH)

# Select biomaRt database and dataset (use biomaRt from current release Ensembl 73)
listMarts()
mart <- useMart(biomart="ensembl")
listDatasets(mart)
mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# Check the arguments required for biomaRt query
listFilters(mart)
head(listAttributes(mart), n=500)
attributePages(mart)

# Build a biomarRt query to obtain ID, name and description for each gene
bioMart_query <- getBM(attributes=c("ensembl_gene_id", "external_gene_id", "description"), filters="ensembl_gene_id", values=RBH$Hsapiens_ensembl_gene_id, mart=mart, uniqueRows=TRUE, bmHeader=TRUE)

# Merge the biomart_query with the RBH table
RBH_annotation <- merge(x=RBH, y=bioMart_query, by.x="Hsapiens_ensembl_gene_id", by.y="ensembl_gene_id", all=TRUE)
head(RBH_annotation)

# Check the merged matrix in terms of size and content
dim(RBH_annotation)
table(duplicated(RBH_annotation[,1]))
table(duplicated(RBH_annotation[,3]))

# Read in custom annotation information for all reference genes in order to check that RBH gene name are not already in reference annotation
sense_info <- read.table(file="gene_info_251113.txt", sep="\t", header=TRUE, fill=TRUE)
sense_info <- sense_info[grep(pattern="^ENSBTAG", x=sense_info[,1], perl=TRUE),1]
length(sense_info)
head(sense_info)

# Select biomaRt database and dataset (use biomaRt from Ensembl 71: Apr 2013)
mart_bos_ref <- useMart(host='apr2013.archive.ensembl.org', biomart="ENSEMBL_MART_ENSEMBL", dataset="btaurus_gene_ensembl")

# Build a biomarRt query to obtain ID, name and description for each gene
bioMart_query_bos_ref <- getBM(attributes=c("ensembl_gene_id", "external_gene_id"), filters="ensembl_gene_id", values=sense_info, mart=mart_bos_ref, uniqueRows=TRUE, bmHeader=TRUE)
head(bioMart_query_bos_ref)

# Check that RBH gene name are not already present in the Bos taurus reference gene name
pseudogene_filt <- RBH_annotation$external_gene_id %in% bioMart_query_bos_ref$external_gene_id
summary(pseudogene_filt)

# All novel genes were kept wether they are putative novel pseudogene or novel gene

# Merge the two different annotation
Annotation <- merge(x=RBH_annotation[,c(1,3,8,9)], y=novel_info, by="novel_gene_id", all=TRUE)
head(Annotation)
dim(Annotation)
table(duplicated(Annotation[,1]))

# Merge the biomart_query with the count table
Annotated_count <- merge(x=Annotation, y=Count$counts, by.x="novel_gene_id", by.y=0, all=TRUE)
head(Annotated_count)

# Check the merged matrix in terms of size and content
dim(Annotated_count)
table(duplicated(Annotated_count[,1]))
length(unique(Annotated_count[,1], incomparables=FALSE))

##########################################
# Create a DGElist of all samples counts #
##########################################

# Create a target matrix and a experimental group vector and animal block
target <- read.delim(file="target.txt", header=TRUE,sep="\t")
group <- factor(paste(target$Treatment, target$Time_point, sep="."))
animal <- factor(target$Animal)

# Create a DGElist containing the group information
Alv_dgelist <- DGEList(counts=Annotated_count[,(ncol(Annotation)+1):ncol(Annotated_count)], lib.size=NULL, norm.factors=NULL, group=group, genes=Annotated_count[,1:ncol(Annotation)], remove.zeros=FALSE)
rownames(Alv_dgelist$counts) <- rownames(Alv_dgelist$genes) <- Alv_dgelist$genes$novel_gene_id
Alv_dgelist$genes$novel_gene_id <- NULL
names(Alv_dgelist)
head(Alv_dgelist$samples)
head(Alv_dgelist$counts)
head(Alv_dgelist$genes)

###########################################################
# Quality check of libraries by plotting density of count #
###########################################################

# Log2 transform the count data for better visualization
count_log2 <- log2(x=(Alv_dgelist$counts[,1:ncol(Alv_dgelist$counts)]+1))

# Plot density of count for all libraries
png(filename="Density_Alv_mac.png", width=1366, height=768, units="px")
plot(x=density(count_log2[,1]), main="Density plot of count per gene", lty=1, xlab="Log2 of count per gene", ylab="Density", col="black", ylim=c(0.0, 0.6))
for (i in 2:ncol(count_log2)) {
  lines(density(count_log2[,i]), lty=1, col="black")
}
dev.off()

#####################################
# Filtering of lowly expressed tags #
#####################################

# Filter non expressed tags (all samples have zero counts)
Alv_no_zeros <- Alv_dgelist[rowSums(Alv_dgelist$counts) > 0,]
dim(Alv_no_zeros$counts)
head(Alv_no_zeros$samples)

# Filter lowly expressed tags, retaining only tags with at least 30 counts (30 counts because it would correspond roughly to the sense DEG filtering) in 10 or more libraries (10 libraries correspond to one group of treatment)
Alv_filt <- Alv_no_zeros[rowSums(cpm(Alv_no_zeros)>30) >=10,]
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

# Log2 transform the filtered count data for better visualization
count_filt_log2 <- log2(Alv_norm$counts[,1:ncol(Alv_norm)]+1)

# Plot density of count for all libraries
png(filename="Density_Alv_mac_post_filter.png", width=1366, height=768, units="px")
plot(density(count_filt_log2[,1]), main="Density plot of count per gene", lty=1, xlab="Log2 of count per gene", ylab="Density", col="black", ylim=c(0.0,0.26))
for (i in 2:ncol(count_filt_log2)) {
  lines(density(count_filt_log2[,i]), lty=1, col="black")
}
dev.off()

################################################
# Multidimensional scaling plot on all samples #
################################################

# Code below will output the MDS plot in directory as a png file
png(filename="MDS_Alv_mac.png", width=1366, height=768, units="px")
MDS <- plotMDS(x=Alv_norm, top=1000000, gene.selection="pairwise", xlab="Dimension 1", ylab="Dimension 2")
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
legend(x="topleft", legend=c("CN", "MB", "2H", "6H", "24H", "48H"), col=c("black", "black", "gold", "green", "blue", "red"), pch=c(1,3,15,15,15,15), cex=2.0)
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

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using the design matrix and calculated dispersion
Alv_fit <- glmFit(y=Alv_disp, design=design)
names(Alv_fit)

################################
# Differential expression call #
################################

# Create a function to perform the differential expression within edgeR according to provided parameters
diff_expr_edgeR <- function(treat1, treat2, data, design, group, adjpvalue, method, LRTdata, DEdata, DEfile, Smearfile) {
  if(length(grep(pattern=treat1, x=colnames(design)))==1 && length(grep(pattern=treat2, x=colnames(design)))==1) {
    contr <- rep(x=0, times=length(colnames(design)))
    contr[c(grep(pattern=treat1, x=colnames(design)), grep(pattern=treat2, x=colnames(design)))] <- c(1, -1)
    lrt <- glmLRT(glmfit=data, contrast=contr)
  }
  else if(length(grep(pattern=treat1, x=colnames(design)))==0 && length(grep(pattern=treat1, x=levels(group)))==1) {
    contr <- grep(pattern=treat2, x=colnames(design))
    lrt <- glmLRT(glmfit=data, coef=contr)
  }
  else if(length(grep(pattern=treat2, x=colnames(design)))==0 && length(grep(pattern=treat2, x=levels(group)))==1) {
    contr <- grep(pattern=treat1, x=colnames(design))
    lrt <- glmLRT(glmfit=data, coef=contr)
  }
  else {
    stop("Error: Check that the treatments provided are in group table!")
  }
  de <- topTags(object=lrt, n="inf", adjust.method=method)
  print("Names of the edgeR likeli-hood ratio test dataframe:")
  print(names(lrt))
  print("Comparison perfomed in the edgeR likeli-hood ratio test:")
  print(lrt$comparison)
  print("Heading of the edgeR likeli-hood ratio test dataframe:")
  print(head(lrt$table))
  print("Summary of the number od edgeR DEG:")
  print(summary(decideTestsDGE(lrt, p.value=adjpvalue)))
  print("Names of the edgeR multiple correction test dataframe:")
  print(names(de))
  print("Heading of the edgeR multiple correction test dataframe:")
  print(head(de$table))
  write.table(x=de$table[,c("Hsapiens_ensembl_gene_id", "external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  png(filename=paste(Smearfile, "png", sep="."), width=1366, height=768, units="px")
  plotSmear(object=lrt, de.tags=(rownames(lrt$table)[as.logical(decideTestsDGE(lrt, p.value = 0.05))]))
  abline(h=c(-1, 1), col="blue")
  dev.off()
  assign(x=LRTdata, value=lrt, envir=.GlobalEnv)
  assign(x=DEdata, value=de, envir=.GlobalEnv)
}

# Test for differential expression between the different time points/treatments
diff_expr_edgeR(treat1="MB.2H", treat2="CN.2H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.2H", DEdata="DE_MB.2H", DEfile="DE_MB_novel.2H", Smearfile="Smear_FC_CPM_MB.2H")
diff_expr_edgeR(treat1="MB.6H", treat2="CN.6H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.6H", DEdata="DE_MB.6H", DEfile="DE_MB_novel.6H", Smearfile="Smear_FC_CPM_MB.6H")
diff_expr_edgeR(treat1="MB.24H", treat2="CN.24H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.24H", DEdata="DE_MB.24H", DEfile="DE_MB_novel.24H", Smearfile="Smear_FC_CPM_MB.24H")
diff_expr_edgeR(treat1="MB.48H", treat2="CN.48H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.48H", DEdata="DE_MB.48H", DEfile="DE_MB_novel.48H", Smearfile="Smear_FC_CPM_MB.48H")

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
colnames(Full_DE)[1] <- "novel_gene_id"
head(Full_DE)

# Write into a table the full DE call data
write.matrix(x=Full_DE, file="Full_DE_novel.txt", sep = "\t")

#######
# END #
#######