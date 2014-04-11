###########################################################################################################################
# RNA-seq data analysis of antisense gene counts (paired data) using annotation Bos_taurus.UMD3.1.71_antisense_251113.gtf #
###########################################################################################################################

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
setwd("F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Antisense_gene")
getwd()
workDir <- getwd()
workDir

################################################
# Read in and concatenate input files within R #
################################################

# Create vector of all files name
fileDir <- "F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Antisense_gene/Counts"
files <- list.files(path=fileDir, pattern="*(2|6|24|48)H$", all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE)
files

# Reads and merges a set of files containing counts
Count <- readDGE(files=files, path=fileDir, columns=c(1,3))
names(Count)
head(Count$samples)
head(Count$counts)

# Ouptut samples data
write.table(x=Count$samples, file="Alv_samples.txt",  sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

######################################################################
# Gene annotation using custom antisense annotation info and biomaRt #
######################################################################

# Read in custom annotation information for all antisense genes
antisense_info <- read.table(file="gene_info_251113.txt", sep="\t", header=TRUE, fill=TRUE)
antisense_info <- antisense_info[grep(pattern="NATENSBTAG", x=antisense_info[,1]),]
antisense_info[,ncol(antisense_info)] <- gsub(pattern="NATENSBTAG", replacement="ENSBTAG", x=antisense_info[,1])
colnames(antisense_info)[ncol(antisense_info)] <- "sense_ensembl_gene_id"
dim(antisense_info)
head(antisense_info)

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
bioMart_query <- getBM(attributes=c("ensembl_gene_id", "external_gene_id", "description"), filters="ensembl_gene_id", values=gsub(pattern="NATENSBTAG", replacement="ENSBTAG", x=rownames(Count$counts), perl=TRUE), mart=mart, uniqueRows=TRUE, bmHeader=TRUE)
bioMart_query_entrezgene <- getBM(attributes=c("ensembl_gene_id", "entrezgene"), filters="ensembl_gene_id", values=gsub(pattern="NATENSBTAG", replacement="ENSBTAG", x=rownames(Count$counts), perl=TRUE), mart=mart, uniqueRows=TRUE, bmHeader=TRUE)

# Format the biomart_query so that there is only one row per ensembl gene ID
bioMart_entrezgene <- aggregate(entrezgene ~ ensembl_gene_id, FUN = "as.vector", data=bioMart_query_entrezgene, na.action="as.vector")
head(bioMart_entrezgene, n=100)

# Merge the different bioMart queries and custom annotation
bioMart_annotation <- merge(x=bioMart_query, y=antisense_info, by.x="ensembl_gene_id", by.y="sense_ensembl_gene_id", all=TRUE)
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
Annotated_count <- merge(x=bioMart_annotation, y=Count$counts, by.x="bostaurus_gene_id", by.y=0, all=TRUE)
colnames(Annotated_count)[2] <- "sense_ensembl_gene_id"
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
Alv_dgelist <- DGEList(counts=Annotated_count[,(ncol(bioMart_annotation)+1):ncol(Annotated_count)], lib.size=NULL, norm.factors=NULL, group=group, genes=Annotated_count[,1:ncol(bioMart_annotation)], remove.zeros=FALSE)
rownames(Alv_dgelist$counts) <- rownames(Alv_dgelist$genes) <- Alv_dgelist$genes$bostaurus_gene_id
Alv_dgelist$genes$bostaurus_gene_id <- NULL
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
plot(x=density(count_log2[,1]), main="Density plot of count per gene", lty=1, xlab="Log2 of count per gene", ylab="Density", col="black", ylim=c(0.0, 1.05))
for (i in 2:ncol(count_log2)) {
  lines(density(count_log2[,i]), lty=1, col="black")
}
dev.off()

##########################################################################
# Read in and concatenate input files within R to perform ratio analyses #
##########################################################################

# Create vector of all files name
fileDir <- "F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Antisense_gene/ratio_sense"
files <- list.files(path=fileDir, pattern="*(2|6|24|48)H$", all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE)
files

# Reads and merges a set of files containing counts for sense data over full gene length (as per antisense annotation)
sense_gene <- readDGE(files=files, path=fileDir, columns=c(1,3))
sense_gene <- as.data.frame(x=sense_gene$counts)
head(sense_gene)
dim(sense_gene)

# Create vector of all files name
fileDir <- "F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Antisense_gene/ratio_antisense"
files <- list.files(path=fileDir, pattern="*(2|6|24|48)H$", all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE)
files

# Reads and merges a set of files containing counts for antisense data over exon (as per sense reference annotation)
antisense_exon <- readDGE(files=files, path=fileDir, columns=c(1,3))
antisense_exon <- as.data.frame(x=antisense_exon$counts)
head(antisense_exon)
dim(antisense_exon)

# Read in the sense raw counts data file
sense_exon <- read.table(file="Alv_rawcounts_sense_for_ratio.txt", header=TRUE, colClasses=c("character", rep(x="numeric", 78)), sep="\t")
head(sense_exon)
dim(sense_exon)

############################################
# Sense-antisense raw counts data merging  #
############################################

# Merge the sense raw counts with the antisense raw counts calculated based on full gene length
AS_S_gene <- merge(x=Count$counts, y=sense_gene, by="row.names")
colnames(AS_S_gene) <- gsub(pattern=".x$", replacement="_AS", x=colnames(AS_S_gene), perl=TRUE)
colnames(AS_S_gene) <- gsub(pattern=".y$", replacement="_S", x=colnames(AS_S_gene), perl=TRUE)
rownames(AS_S_gene) <- AS_S_gene[,1]
AS_S_gene <- AS_S_gene[,-1]
rownames(AS_S_gene) <- gsub(pattern="NAT", replacement="", x=rownames(AS_S_gene))
head(AS_S_gene)
dim(AS_S_gene)

# Merge the sense raw counts with the antisense raw counts calculated based on exon only
AS_S_exon <- merge(x=antisense_exon, y=sense_exon, by.x="row.names", by.y="ensembl_gene_id")
colnames(AS_S_exon) <- gsub(pattern=".x$", replacement="_AS", x=colnames(AS_S_exon), perl=TRUE)
colnames(AS_S_exon) <- gsub(pattern=".y$", replacement="_S", x=colnames(AS_S_exon), perl=TRUE)
rownames(AS_S_exon) <- AS_S_exon[,1]
AS_S_exon <- AS_S_exon[,-1]
gene_filter <- rownames(AS_S_exon) %in% rownames(AS_S_gene)
AS_S_exon <- AS_S_exon[gene_filter,]
head(AS_S_exon)
dim(AS_S_exon)

# Perform calculation of counts for sense and antisense based on intron only
#AS_S_intron <- merge(x=AS_S_gene, y=AS_S_exon, by="row.names")
#colnames(AS_S_intron) <- gsub(pattern=".x$", replacement="_gene", x=colnames(AS_S_intron), perl=TRUE)
#colnames(AS_S_intron) <- gsub(pattern=".y$", replacement="_exon", x=colnames(AS_S_intron), perl=TRUE)
#rownames(AS_S_intron) <- AS_S_intron[,1]
#AS_S_intron <- AS_S_intron[,-1]
#head(AS_S_intron)
#dim(AS_S_intron)
#AS_S_intron <- as.data.frame(AS_S_intron)
#for (i in 1:(ncol(AS_S_gene))) {
#  if (gsub(pattern="_gene$", replacement="", x=colnames(AS_S_intron)[i], perl=TRUE) == gsub(pattern="_exon$", replacement="", x=colnames(AS_S_intron)[i+(ncol(AS_S_gene))], perl=TRUE)) {
#    AS_S_intron <- cbind(AS_S_intron, (AS_S_intron[,i]-AS_S_intron[,(i+(ncol(AS_S_gene)))]))
#    colnames(AS_S_intron)[i+(ncol(AS_S_gene)*2)] <- gsub(pattern="_gene$", replacement="", x=colnames(AS_S_intron)[i], perl=TRUE)
#  }
#  else {
#    stop("Column names do not match! Abort!")
#  }
#}
#AS_S_intron <- AS_S_intron[,c(-(1:(ncol(AS_S_gene)*2)))]
#head(AS_S_intron)
#dim(AS_S_intron)

##################################################
# Filter out all lowly expressed antisense genes #
##################################################

# Filter lowly expressed tags, retaining only tags with at least 20 count per million in 10 or more libraries (10 libraries correspond to one group of treatment)
Ratio_filt <- Alv_dgelist[rowSums(cpm(Alv_dgelist)>20) >=10,]
Ratio_filt <- Ratio_filt$counts[,c(-(1:78))]
rownames(Ratio_filt) <- gsub(pattern="NAT", replacement="", x=rownames(Ratio_filt))
head(Ratio_filt)
dim(Ratio_filt)

# Keep only the genes which passed the low expression filtering from the different matrices
AS_S_gene_filt <- merge(x=Ratio_filt, y=AS_S_gene, by="row.names")
rownames(AS_S_gene_filt) <- AS_S_gene_filt[,1]
AS_S_gene_filt <- AS_S_gene_filt[,-1]
head(AS_S_gene_filt)
dim(AS_S_gene_filt)
AS_S_exon_filt <- merge(x=Ratio_filt, y=AS_S_exon, by="row.names")
rownames(AS_S_exon_filt) <- AS_S_exon_filt[,1]
AS_S_exon_filt <- AS_S_exon_filt[,-1]
head(AS_S_exon_filt)
dim(AS_S_exon_filt)
#AS_S_intron_filt <- merge(x=Ratio_filt, y=AS_S_intron, by="row.names")
#rownames(AS_S_intron_filt) <- AS_S_intron_filt[,1]
#AS_S_intron_filt <- AS_S_intron_filt[,-1]
#head(AS_S_intron_filt)
#dim(AS_S_intron_filt)

###############################################################
# Antisense-sense ratio calculation over the full gene length #
###############################################################

# Define the matrix which will contain the ratio values
ratio_gene <- matrix(nrow=nrow(AS_S_gene_filt), ncol=((ncol(AS_S_gene_filt)/2)*3))
rownames(ratio_gene) <- rownames(AS_S_gene_filt)
ratio_gene <- as.data.frame(ratio_gene)

# Create a loop to calculate ratios for sense and antisense
for (i in 1:(ncol(AS_S_gene_filt)/2)) {
  if (gsub(pattern="_AS$", replacement="", x=colnames(AS_S_gene_filt)[i], perl=TRUE) == gsub(pattern="_S$", replacement="", x=colnames(AS_S_gene_filt)[i+(ncol(AS_S_gene_filt)/2)], perl=TRUE)) {
    ratio_gene[,(i*3-2)] <- (AS_S_gene_filt[,i]*100/rowSums(AS_S_gene_filt[,c(i,(i+(ncol(AS_S_gene_filt)/2)))]))
    ratio_gene[,(i*3-1)] <- (AS_S_gene_filt[,(i+(ncol(AS_S_gene_filt)/2))]*100/rowSums(AS_S_gene_filt[,c(i,(i+(ncol(AS_S_gene_filt)/2)))]))
    ratio_gene[,(i*3)] <- rowSums(AS_S_gene_filt[,c(i,(i+(ncol(AS_S_gene_filt)/2)))])
    colnames(ratio_gene)[i*3-2] <- gsub(pattern="_AS$", replacement="_AS_ratio", x=colnames(AS_S_gene_filt)[i], perl=TRUE)
    colnames(ratio_gene)[i*3-1] <- gsub(pattern="_AS$", replacement="_S_ratio", x=colnames(AS_S_gene_filt)[i], perl=TRUE)
    colnames(ratio_gene)[i*3] <- gsub(pattern="_AS$", replacement="_AS_S_total", x=colnames(AS_S_gene_filt)[i], perl=TRUE) 
  }
  else {
    stop("Column names do not match! Abort!")
  }
}
head(ratio_gene)
dim(ratio_gene)

# Calculate the average and median of antisense ratio
ratio_gene_average <- vector()
for (i in 1:ncol(ratio_gene)) {
  if (length(grep(pattern="_AS_ratio$", x=colnames(ratio_gene)[i], perl=TRUE))==1) {
    ratio_gene_average <- c(ratio_gene_average, ratio_gene[,i])
  }
}
length(ratio_gene_average)
mean(x=ratio_gene_average,  trim=0, na.rm=TRUE)
median(x=ratio_gene_average, na.rm=TRUE)

# Output data
write.table(x=ratio_gene, file="Antisense_sense_ratio(over_gene).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

########################################################
# Antisense-sense ratio calculation over the exon only #
########################################################

# Define the matrix which will contain the ratio values
ratio_exon <- matrix(nrow=nrow(AS_S_exon_filt), ncol=((ncol(AS_S_exon_filt)/2)*3))
rownames(ratio_exon) <- rownames(AS_S_exon_filt)
ratio_exon <- as.data.frame(ratio_exon)

# Create a loop to calculate ratios for sense and antisense
for (i in 1:(ncol(AS_S_exon_filt)/2)) {
  if (gsub(pattern="_AS$", replacement="", x=colnames(AS_S_exon_filt)[i], perl=TRUE) == gsub(pattern="_S$", replacement="", x=colnames(AS_S_exon_filt)[i+(ncol(AS_S_exon_filt)/2)], perl=TRUE)) {
    ratio_exon[,(i*3-2)] <- (AS_S_exon_filt[,i]*100/rowSums(AS_S_exon_filt[,c(i,(i+(ncol(AS_S_exon_filt)/2)))]))
    ratio_exon[,(i*3-1)] <- (AS_S_exon_filt[,(i+(ncol(AS_S_exon_filt)/2))]*100/rowSums(AS_S_exon_filt[,c(i,(i+(ncol(AS_S_exon_filt)/2)))]))
    ratio_exon[,(i*3)] <- rowSums(AS_S_exon_filt[,c(i,(i+(ncol(AS_S_exon_filt)/2)))])
    colnames(ratio_exon)[i*3-2] <- gsub(pattern="_AS$", replacement="_AS_ratio", x=colnames(AS_S_exon_filt)[i], perl=TRUE)
    colnames(ratio_exon)[i*3-1] <- gsub(pattern="_AS$", replacement="_S_ratio", x=colnames(AS_S_exon_filt)[i], perl=TRUE)
    colnames(ratio_exon)[i*3] <- gsub(pattern="_AS$", replacement="_AS_S_total", x=colnames(AS_S_exon_filt)[i], perl=TRUE) 
  }
  else {
    stop("Column names do not match! Abort!")
  }
}
head(ratio_exon)
dim(ratio_exon)

# Calculate the average and median of antisense ratio
ratio_exon_average <- vector()
for (i in 1:ncol(ratio_exon)) {
  if (length(grep(pattern="_AS_ratio$", x=colnames(ratio_exon)[i], perl=TRUE))==1) {
    ratio_exon_average <- c(ratio_exon_average, ratio_exon[,i])
  }
}
length(ratio_exon_average)
mean(x=ratio_exon_average,  trim=0, na.rm=TRUE)
median(x=ratio_exon_average, na.rm=TRUE)

# Output data
write.table(x=ratio_exon, file="Antisense_sense_ratio(over_exon).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

###############################################################################
# Ratio summarisation between Control and M. bovis samples within time points #
###############################################################################

# Create a function to summarise the ratio across treatments and time points
ratio_summary <- function(ratio, arg) {
  feature <- gsub(pattern="ratio_", replacement="", x=deparse(substitute(ratio)))
  DE_ratio <- matrix(nrow=nrow(ratio), ncol=4)
  rownames(DE_ratio) <- rownames(ratio)
  DE_ratio <- as.data.frame(DE_ratio)
  MB <- paste("MB_", arg, "_AS_ratio", sep="")
  MB_tot <- paste("MB_", arg, "_AS_S_total", sep="")
  CN <- paste("CN_", arg, "_AS_ratio", sep="")
  CN_tot <- paste("CN_", arg, "_AS_S_total", sep="")
  DE_ratio[,1] <- rowMeans(x=ratio[,grep(pattern=MB, x=colnames(ratio), perl=TRUE)], na.rm=TRUE)
  DE_ratio[,2] <- rowMeans(x=ratio[,grep(pattern=MB_tot, x=colnames(ratio), perl=TRUE)], na.rm=TRUE)
  DE_ratio[,3] <- rowMeans(x=ratio[,grep(pattern=CN, x=colnames(ratio), perl=TRUE)], na.rm=TRUE)
  DE_ratio[,4] <- rowMeans(x=ratio[,grep(pattern=CN_tot, x=colnames(ratio), perl=TRUE)], na.rm=TRUE)
  MBmean <- paste("mean_", MB, sep="")
  CNmean <- paste("mean_", CN, sep="")
  MB_tot_mean <- paste("mean_", MB_tot, sep="")
  CN_tot_mean <- paste("mean_", CN_tot, sep="")
  colnames(DE_ratio)[1:4] <- c(MBmean, MB_tot_mean, CNmean, CN_tot_mean)
  png(filename=paste("Ratio_vs_total_", feature, "_", arg, ".png", sep=""), width=500, height=500, units="px")
  plot(x=DE_ratio[,2], y=DE_ratio[,1], xlab="Total count per gene", ylab="% of antisense per gene", log="x", col=rgb(210,75,60,100,maxColorValue=255))
  points(x=DE_ratio[,4], y=DE_ratio[,3], col=rgb(70,180,90,100,maxColorValue=255))
  abline(h=c(20, 80), col="blue")
  legend(x="topright", legend=c("M. bovis-infected", "Control-uninfected"), fill=c(rgb(210,75,60,100,maxColorValue=255), rgb(70,180,90,100,maxColorValue=255)))
  dev.off()
  return(DE_ratio)
}

# Calculate the mean ratio for each conditions and time points
DE_ratio_gene <- ratio_summary(ratio=ratio_gene, arg="2H")
DE_ratio_gene <- cbind(DE_ratio_gene, ratio_summary(ratio=ratio_gene, arg="6H"))
DE_ratio_gene <- cbind(DE_ratio_gene, ratio_summary(ratio=ratio_gene, arg="24H"))
DE_ratio_gene <- cbind(DE_ratio_gene, ratio_summary(ratio=ratio_gene, arg="48H"))
head(DE_ratio_gene)
dim(DE_ratio_gene)
DE_ratio_exon <- ratio_summary(ratio=ratio_exon, arg="2H")
DE_ratio_exon <- cbind(DE_ratio_exon, ratio_summary(ratio=ratio_exon, arg="6H"))
DE_ratio_exon <- cbind(DE_ratio_exon, ratio_summary(ratio=ratio_exon, arg="24H"))
DE_ratio_exon <- cbind(DE_ratio_exon, ratio_summary(ratio=ratio_exon, arg="48H"))
head(DE_ratio_exon)
dim(DE_ratio_exon)

# Perform wilcoxon test to identify genes with significantly different ratio between exon and gene
DE_ratio_gene[,"p_value"] <- as.vector(x=rep(x=1, 4830), mode="numeric")
for (i in 1:nrow(DE_ratio_gene)) {
  if ((rownames(ratio_gene)[i]==rownames(ratio_exon)[i])&&(rownames(DE_ratio_gene)[i]==rownames(ratio_gene)[i])) {
    if (is.na(mean(x=as.vector(x=ratio_gene[i,grep(pattern="_AS_ratio", x=colnames(ratio_gene))], mode="numeric"), na.rm=TRUE)) | is.na(mean(x=as.vector(x=ratio_exon[i,grep(pattern="_AS_ratio", x=colnames(ratio_exon))], mode="numeric"), na.rm=TRUE))) {
      DE_ratio_gene[i,"p_value"] <- "NaN"
    }
    else {
      wilcox <- wilcox.test(x=(as.vector(x=ratio_gene[i,grep(pattern="_AS_ratio", x=colnames(ratio_gene))], mode="numeric")), y=(as.vector(x=ratio_exon[i,grep(pattern="_AS_ratio", x=colnames(ratio_exon))], mode="numeric")), alternative="greater", paired=TRUE, exact=FALSE)
      DE_ratio_gene[i,"p_value"] <- wilcox$p.value
    }
  }
  else {
    stop("Non matching rownames!")
  }
}
DE_ratio_gene$p_value <- as.numeric(DE_ratio_gene$p_value)
DE_ratio_gene[,"adj_p_value"] <- as.vector(x=p.adjust(p=as.numeric(DE_ratio_gene$p_value), method="BH"), mode="numeric")
head(DE_ratio_gene)
tail(DE_ratio_gene)
dim(DE_ratio_gene)
summary(DE_ratio_gene$adj_p_value < 0.05)

# Output data
write.table(x=DE_ratio_gene, file="DE_ratio_(over_gene).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=DE_ratio_exon, file="DE_ratio_(over_exon).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#############################################
# Plotting of antisense versus sense counts #
#############################################

# Create the color filter
color <- vector()
for (i in 1:nrow(AS_S_gene_filt)) {
  if (DE_ratio_gene[row.names(AS_S_gene_filt[i,]),"adj_p_value"] == "NaN") {
    color <- c(color, "grey")
  }
  else if (DE_ratio_gene[row.names(AS_S_gene_filt[i,]),"adj_p_value"] >= 0.05) {
    color <- c(color, "black")
  }
  else if (DE_ratio_gene[row.names(AS_S_gene_filt[i,]),"adj_p_value"] < 0.05) {
    color <- c(color, "red")
  }
  else {
    stop("Unauthorised value!")
  }
}
length(color)

# Plot the sense cpm versus the antisense cpm
png(filename="CPM_sense_antisense.1.png", width=1000, height=1000, units="px")
par(mfrow=c(4,4))
for (i in 1:16) {
  col_antisense <- gsub(pattern="_AS", replacement="", x=colnames(AS_S_gene_filt)[i], perl=TRUE)
  col_sense <- gsub(pattern="$", replacement="_S", x=col_antisense, perl=TRUE)
  plot(x=log(x=(cpm(AS_S_gene_filt[,i])+1), base=10), y=log(x=(cpm(AS_S_gene_filt[,grep(pattern=col_sense, x=colnames(AS_S_gene_filt))])+1), base=10), xlab="Log10_CPM_antisense", ylab="Log10_CPM_sense", main=col_antisense, col=color)
}
dev.off()

png(filename="CPM_sense_antisense.2.png", width=1000, height=1000, units="px")
par(mfrow=c(4,4))
for (i in 17:32) {
  col_antisense <- gsub(pattern="_AS", replacement="", x=colnames(AS_S_gene_filt)[i], perl=TRUE)
  col_sense <- gsub(pattern="$", replacement="_S", x=col_antisense, perl=TRUE)
  plot(x=log(x=(cpm(AS_S_gene_filt[,i])+1), base=10), y=log(x=(cpm(AS_S_gene_filt[,grep(pattern=col_sense, x=colnames(AS_S_gene_filt))])+1), base=10), xlab="Log10_CPM_antisense", ylab="Log10_CPM_sense", main=col_antisense, col=color)
}
dev.off()

png(filename="CPM_sense_antisense.3.png", width=1000, height=1000, units="px")
par(mfrow=c(4,4))
for (i in 33:48) {
  col_antisense <- gsub(pattern="_AS", replacement="", x=colnames(AS_S_gene_filt)[i], perl=TRUE)
  col_sense <- gsub(pattern="$", replacement="_S", x=col_antisense, perl=TRUE)
  plot(x=log(x=(cpm(AS_S_gene_filt[,i])+1), base=10), y=log(x=(cpm(AS_S_gene_filt[,grep(pattern=col_sense, x=colnames(AS_S_gene_filt))])+1), base=10), xlab="Log10_CPM_antisense", ylab="Log10_CPM_sense", main=col_antisense, col=color)
}
dev.off()

png(filename="CPM_sense_antisense.4.png", width=1000, height=1000, units="px")
par(mfrow=c(4,4))
for (i in 49:64) {
  col_antisense <- gsub(pattern="_AS", replacement="", x=colnames(AS_S_gene_filt)[i], perl=TRUE)
  col_sense <- gsub(pattern="$", replacement="_S", x=col_antisense, perl=TRUE)
  plot(x=log(x=(cpm(AS_S_gene_filt[,i])+1), base=10), y=log(x=(cpm(AS_S_gene_filt[,grep(pattern=col_sense, x=colnames(AS_S_gene_filt))])+1), base=10), xlab="Log10_CPM_antisense", ylab="Log10_CPM_sense", main=col_antisense, col=color)
}
dev.off()

png(filename="CPM_sense_antisense.5.png", width=1000, height=1000, units="px")
par(mfrow=c(4,4))
for (i in 65:78) {
  col_antisense <- gsub(pattern="_AS", replacement="", x=colnames(AS_S_gene_filt)[i], perl=TRUE)
  col_sense <- gsub(pattern="$", replacement="_S", x=col_antisense, perl=TRUE)
  plot(x=log(x=(cpm(AS_S_gene_filt[,i])+1), base=10), y=log(x=(cpm(AS_S_gene_filt[,grep(pattern=col_sense, x=colnames(AS_S_gene_filt))])+1), base=10), xlab="Log10_CPM_antisense", ylab="Log10_CPM_sense", main=col_antisense, col=color)
}
dev.off()

####################################################################################
# Filter out all antisense genes not identified by the differential ratio analysis #
####################################################################################

# Filter out antisense genes
Alv_filt <- Alv_dgelist[gsub(pattern="ENSBTAG", replacement="NATENSBTAG", x=rownames(DE_ratio_gene[DE_ratio_gene$adj_p_value < 0.05,])),]
dim(Alv_filt$counts)
head(Alv_filt$samples)
head(Alv_filt$counts)

# Recompute the library size
Alv_filt$samples$lib.size <- colSums(Alv_filt$counts)
head(Alv_filt$samples)
head(Alv_dgelist$samples)

# Output txt file of raw count
write.table(x=Alv_filt$counts, file="Alv_filt_rawcount.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

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
count_filt_log2 <- log2(Alv_norm$counts[,1:ncol(Alv_norm$counts)]+1)

# Plot density of count for all libraries
png(filename="Density_Alv_mac_post_filter.png", width=1366, height=768, units="px")
plot(density(count_filt_log2[,1]), main="Density plot of count per gene", lty=1, xlab="Log2 of count per gene", ylab="Density", col="black", ylim=c(0.0,0.25))
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
png(filename="MDS_Alv_mac(D1_vs_D2).png", width=1600, height=1000, units="px")
MDS <- plotMDS(x=Alv_norm, top=1000000, gene.selection="pairwise", xlab="Dimension 1", ylab="Dimension 4", col=color, labels=symbol, dim.plot=c(1,4), cex=2.0)
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

##################################################################
# Determine differential expression using negative binomial GLMs #
##################################################################

# Fit a negative binomial generalized linear model for each tag using the design matrix and calculated dispersion
Alv_fit <- glmFit(y=Alv_disp, design=design, )
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
  write.table(x=de$table[,c("sense_ensembl_gene_id", "external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  png(filename=paste(Smearfile, "png", sep="."), width=1366, height=768, units="px")
  plotSmear(object=lrt, de.tags=(rownames(lrt$table)[as.logical(decideTestsDGE(lrt, p.value = 0.05))]))
  abline(h=c(-1, 1), col="blue")
  dev.off()
  assign(x=LRTdata, value=lrt, envir=.GlobalEnv)
  assign(x=DEdata, value=de, envir=.GlobalEnv)
}

# Test for differential expression between the different time points/treatments
diff_expr_edgeR(treat1="MB.2H", treat2="CN.2H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.2H", DEdata="DE_MB.2H", DEfile="DE_MB_antisense.2H", Smearfile="Smear_FC_CPM_MB.2H")
diff_expr_edgeR(treat1="MB.6H", treat2="CN.6H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.6H", DEdata="DE_MB.6H", DEfile="DE_MB_antisense.6H", Smearfile="Smear_FC_CPM_MB.6H")
diff_expr_edgeR(treat1="MB.24H", treat2="CN.24H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.24H", DEdata="DE_MB.24H", DEfile="DE_MB_antisense.24H", Smearfile="Smear_FC_CPM_MB.24H")
diff_expr_edgeR(treat1="MB.48H", treat2="CN.48H", data=Alv_fit, design=design, group=group, adjpvalue=0.05, method="BH", LRTdata="Alv_lrt_MB.48H", DEdata="DE_MB.48H", DEfile="DE_MB_antisense.48H", Smearfile="Smear_FC_CPM_MB.48H")

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
colnames(Full_DE)[1] <- "bostaurus_gene_id"
head(Full_DE)

# Write into a table the full DE call data
write.matrix(x=Full_DE, file="Full_DE_antisense.txt", sep = "\t")

#######
# END #
#######