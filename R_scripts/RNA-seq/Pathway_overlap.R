#############################################################
# Pathway analyses of sense-novel gene counts (paired data) #
#############################################################

#####################################################
# Bovine gene annotation by orthology to H. sapiens #
#####################################################

# Read the input bovine gene ID and name
Bovine_to_human <- read.table(file="Annotated_count.txt", sep="\t", header=TRUE, fill=TRUE)
head(Bovine_to_human)
dim(Bovine_to_human)

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

# Read in the sense DE genes data
Full_DE <- read.table(file="Full_DE_sense.txt", sep="\t", header=TRUE, fill=TRUE)

# Merge the H. sapiens ortholog gene ID annotation with DE results
ortholog_full_DE <- merge(x=human_ortholog, y=Full_DE[,c(1,11:30)], by="ensembl_gene_id", all.y=TRUE)
colnames(ortholog_full_DE)[2] <- "Hsapiens_ensembl_gene_id"
ortholog_full_DE <- ortholog_full_DE[(!is.na(x=ortholog_full_DE$Hsapiens_ensembl_gene_id)),-1]
head(ortholog_full_DE)
dim(ortholog_full_DE)

# Read in the novel DE genes data
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

# Create a function to perform Sigora analysis based on user-defined parameters (this function can take the bovine ID and transform them in human ID)
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

# Create a function to perform Sigora analysis based on user-defined parameters (this function takes human ID only)
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

##############################################################################
# Compare the different KEGG pathways obtained with Sigora and Pathway-guide #
##############################################################################

# Load required library
library(gdata)

# Read in the lists of pathways obtained with Pathway-guide from KEGG
Path_MB2H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_2H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2, )
colnames(Path_MB2H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB2H_novel), perl=TRUE)
head(Path_MB2H_novel)
dim(Path_MB2H_novel)
Path_MB6H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_6H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB6H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB6H_novel), perl=TRUE)
head(Path_MB6H_novel)
dim(Path_MB6H_novel)
Path_MB24H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_24H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB24H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB24H_novel), perl=TRUE)
head(Path_MB24H_novel)
dim(Path_MB24H_novel)
Path_MB48H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/KEGG/Pathway_Sig_MB_sense_48H(novel)/Adipocytokine signaling pathway.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB48H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB48H_novel), perl=TRUE)
head(Path_MB48H_novel)
dim(Path_MB48H_novel)

# Create a function to read in the lists of pathways obtained with Sigora according to database
sigora_input <- function(file, database) {
  sigora_path <- read.table(file=file, header=FALSE, sep="\t", quote="", dec=".", na.strings="NA", col.names=c("row", "count",  "weightsum",  "pwyid",  "pwyname",  "p_value",  "adj.p_value",  "round_sum",	"m",	"n",	"k"), fill=TRUE, blank.lines.skip=TRUE, stringsAsFactors=FALSE)
  if (database == "KEGG") {
    sigora_path <- sigora_path[(grep(pattern="TRUE", x=(sigora_path[,2]=="Summary results (based on Kegg):")))+2:nrow(sigora_path),-1]
  }
  else if (database == "Reactome") {
    sigora_path <- sigora_path[((grep(pattern="TRUE", x=(sigora_path[,2]=="Summary results (based on Reactome):")))+2):((grep(pattern="TRUE", x=(sigora_path[,1]=='[1] "KEGG"')))-4),-1]
  }
  else {
    stop("Non authorised value for database! Please try again with either KEGG or Reactome values!")
  }
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

# Read in the lists of KEGG pathways obtained with Sigora
Sig_MB2H_novel <- sigora_input(file="./Sigora/Sigora_MB2H(novel).txt", database="KEGG")
head(Sig_MB2H_novel)
dim(Sig_MB2H_novel)
Sig_MB6H_novel <- sigora_input(file="./Sigora/Sigora_MB6H(novel).txt", database="KEGG")
head(Sig_MB6H_novel)
dim(Sig_MB6H_novel)
Sig_MB24H_novel <- sigora_input(file="./Sigora/Sigora_MB24H(novel).txt", database="KEGG")
head(Sig_MB24H_novel)
dim(Sig_MB24H_novel)
Sig_MB48H_novel <- sigora_input(file="./Sigora/Sigora_MB48H(novel).txt", database="KEGG")
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

# Identify the KEGG common significant pathways between Pathway-guide and Sigora
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

# Output the KEGG common pathways
write.table(x=Pathway_MB2H, file="Pathway_MB2H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB6H, file="Pathway_MB6H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB24H, file="Pathway_MB24H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB48H, file="Pathway_MB48H(kegg).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

##################################################################################
# Compare the different Reactome pathways obtained with Sigora and Pathway-guide #
##################################################################################

# Load required library
library(gdata)

# Read in the lists of pathways obtained with Pathway-guide from Reactome
Path_MB2H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_2H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2, )
colnames(Path_MB2H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB2H_novel), perl=TRUE)
head(Path_MB2H_novel)
dim(Path_MB2H_novel)
Path_MB6H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_6H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB6H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB6H_novel), perl=TRUE)
head(Path_MB6H_novel)
dim(Path_MB6H_novel)
Path_MB24H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_24H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB24H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB24H_novel), perl=TRUE)
head(Path_MB24H_novel)
dim(Path_MB24H_novel)
Path_MB48H_novel <- read.xls(xls="F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/Pathway-guide/Reactome/Pathway_Sig_MB_sense_48H(novel)/B Cell Activation.xls", sheet=1, header=TRUE, row.names=2)
colnames(Path_MB48H_novel) <- gsub(pattern="$", replacement="_Pathguide", x=colnames(Path_MB48H_novel), perl=TRUE)
head(Path_MB48H_novel)
dim(Path_MB48H_novel)

# Read in the lists of Reactome pathways obtained with Sigora (using previously defined function: sigora_input)
Sig_MB2H_novel <- sigora_input(file="./Sigora/Sigora_MB2H(novel).txt", database="Reactome")
head(Sig_MB2H_novel)
dim(Sig_MB2H_novel)
Sig_MB6H_novel <- sigora_input(file="./Sigora/Sigora_MB6H(novel).txt", database="Reactome")
head(Sig_MB6H_novel)
dim(Sig_MB6H_novel)
Sig_MB24H_novel <- sigora_input(file="./Sigora/Sigora_MB24H(novel).txt", database="Reactome")
head(Sig_MB24H_novel)
dim(Sig_MB24H_novel)
Sig_MB48H_novel <- sigora_input(file="./Sigora/Sigora_MB48H(novel).txt", database="Reactome")
head(Sig_MB48H_novel)
dim(Sig_MB48H_novel)

# Identify the Reactome common significant pathways between Pathway-guide and Sigora (using previously defined function: common_path)
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

# Output the Reactome common pathways
write.table(x=Pathway_MB2H, file="Pathway_MB2H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB6H, file="Pathway_MB6H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB24H, file="Pathway_MB24H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
write.table(x=Pathway_MB48H, file="Pathway_MB48H(reactome).txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#######
# END #
#######