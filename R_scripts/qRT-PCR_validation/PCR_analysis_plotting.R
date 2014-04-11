###################################################
# Analysis and plotting of real time qRT-PCR data #
###################################################

#############################
# List of required packages #
#############################

# Source the common functions used across this script
source(file="F:/nnalpas/Documents/PhD project/Bioinformatics/R/General_function.R")

# Load the required packages
loadpackage(package=ggplot2)
loadpackage(package=grid)

###############################
# Prepare working environment #
###############################

# Move to the appropriate folder
setwd("F:/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/PCR validation/Plotting in R(10_animals)")
getwd()
workDir <- getwd()
workDir

###############################
# Read in input file within R #
###############################

# Read in the data file
Data <- read.table(file="PCR_results.txt", header=TRUE, colClasses=c("character", "character", "character", "integer", "numeric", "numeric", "numeric", "numeric"), sep="\t")
head(Data)

# Split the data on a method basis
PCR_CNRQ_filt <- grep(pattern="PCR", x=Data$Methods, perl=TRUE)
PCR_CNRQ <- Data[PCR_CNRQ_filt,]
head(PCR_CNRQ)
RNAseq_filt <- grep(pattern="RNA-seq", x=Data$Methods)
RNAseq <- Data[RNAseq_filt,]
head(RNAseq)

###########################################
# Adjust p-value for multiple comparisons #
###########################################

# Adjust p-value of PCR results for multiple comparisons with Benjamini-Hochberg method
PCR_CNRQ[grep(pattern="0", x=PCR_CNRQ$Times..H., invert=TRUE),"Adj.p.value"] <- p.adjust(p=PCR_CNRQ[grep(pattern="0", x=PCR_CNRQ$Times..H., invert=TRUE),"p.value"], method="BH")
PCR_CNRQ$Adj.p.value

# Add significance label
PCR_CNRQ <- sig_label(arg1=PCR_CNRQ, arg2="Adj.p.value")
head(PCR_CNRQ)
RNAseq <- sig_label(arg1=RNAseq, arg2="Adj.p.value")
head(RNAseq)

# Bind together the PCR and RNA-seq data
Data_2_plot <- rbind(PCR_CNRQ, RNAseq)

#####################
# Plot the PCR data #
#####################

# Plot the data
expression_plot <- function (arg1) {
  gene <- c(arg1$Genes)
  gene <- unique(x=gene)
  for (i in 1:length(gene)) {
    file <- paste(gene[i], "png", sep=".")
    if (length(grep(pattern=gene[i], x=arg1$Genes)) == 0) {
      plot1 <- "Not performed"
    }
    else {
      plot1 <- qplot(data=arg1[grep(pattern=gene[i], x=arg1$Genes),], x=Time_in_Hours, y=Mean_log2_Fold.change, colour=Methods, main=gene[i])+geom_point(shape=19, size=15)+geom_line(size=5)+geom_errorbar(aes(x=Time_in_Hours, ymin=Mean_log2_Fold.change-Std._Error_Mean, ymax=Mean_log2_Fold.change+Std._Error_Mean), width=1, size=2.5)+geom_text(aes(x=Time_in_Hours+1.75, y=Mean_log2_Fold.change+Std._Error_Mean, label=Significance_label), size=22)+theme(panel.background=element_rect(fill='wheat'), legend.title=element_text(size=45, face="bold"), legend.text=element_text(size=45, face="bold"), axis.title.x=element_text(face="bold", size=45), axis.text.x=element_text(face="bold", size=45), axis.title.y=element_text(face="bold", size=45), axis.text.y=element_text(face="bold", size=45), plot.title=element_text(face="bold", size=55), legend.position=c(0.33, 0.08))+xlab("Time (hours)")+ylab("Log2 fold-change")
    }
    png(filename=file, width=1366, height=1366, units="px")
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(1, 1, heights = unit(c(4), "null"))))
    if (plot1 == "Not performed") {
      grid.text(label="Not performed", vp=viewport(layout.pos.row=1, layout.pos.col=1))
    }
    else {
      print(plot1, vp=viewport(layout.pos.row=1, layout.pos.col=1))
    }
    dev.off()
  }
}
colnames(Data_2_plot) <- c("Genes", "Methods", "Comparisons", "Time_in_Hours", "Mean_log2_Fold.change", "Std._Error_Mean", "p_value", "Adj_p_value", "Significance_label")
PCR_vs_RNAseq <- expression_plot(Data_2_plot)

#######################################
# Correlation between RNA-seq and PCR #
#######################################

# Filter out the log2 fold-change at 0H from both datasets
PCR_cor_filt <- grep(pattern="0", x=PCR_CNRQ$Times..H., invert=TRUE)
PCR_cor_no0 <- PCR_CNRQ[PCR_cor_filt,]
head(PCR_cor_no0)
RNAseq_cor_filt <- grep(pattern="0", x=RNAseq$Times..H., invert=TRUE)
RNAseq_cor_no0 <- RNAseq[RNAseq_cor_filt,]
head(RNAseq_cor_no0)

# Calculate correlation between fold-changes obtained with PCR and RNA-seq without time point 0H
cor.test(PCR_cor_no0$Mean, RNAseq_cor_no0$Mean, method="pearson")
cor.test(PCR_cor_no0$Mean, RNAseq_cor_no0$Mean, method="spearman")

# Calculate correlation between fold-changes obtained with PCR and RNA-seq for each gene
cor.test(PCR_cor_no0[grep(pattern="CCL4", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="CCL4", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="CCL4", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="CCL4", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")
cor.test(PCR_cor_no0[grep(pattern="FOS", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="FOS", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="FOS", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="FOS", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")
cor.test(PCR_cor_no0[grep(pattern="IL10", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="IL10", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="IL10", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="IL10", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")
cor.test(PCR_cor_no0[grep(pattern="IL1B", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="IL1B", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="IL1B", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="IL1B", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")
cor.test(PCR_cor_no0[grep(pattern="IL6", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="IL6", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="IL6", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="IL6", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")
cor.test(PCR_cor_no0[grep(pattern="PIK3IP1", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="PIK3IP1", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="PIK3IP1", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="PIK3IP1", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")
cor.test(PCR_cor_no0[grep(pattern="TLR2", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="TLR2", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="TLR2", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="TLR2", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")
cor.test(PCR_cor_no0[grep(pattern="TNF", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="TNF", x=RNAseq_cor_no0$Samples),"Mean"], method="pearson")
cor.test(PCR_cor_no0[grep(pattern="TNF", x=PCR_cor_no0$Samples),"Mean"], RNAseq_cor_no0[grep(pattern="TNF", x=RNAseq_cor_no0$Samples),"Mean"], method="spearman")

#######
# END #
#######