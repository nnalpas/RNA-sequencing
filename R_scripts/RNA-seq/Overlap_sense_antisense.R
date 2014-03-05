#######################################################################
# Compare the list of differentially expressed genes between datasets #
#######################################################################

# Set the working directory
setwd(dir="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Overlap")
getwd()

# Read the input sense DEG at the different time points
de_sense <- read.table(file="Full_DE_sense.txt", header=TRUE, sep="\t", quote="")
head(de_sense)
dim(de_sense)

# Read the input antisense DEG at the different time points
de_antisense <- read.table(file="Full_DE_antisense.txt", header=TRUE, sep="\t", quote="")
head(de_antisense)
dim(de_antisense)

# Create a function to find common gene between sense and antisense
sense_antisense_overlap <- function(time) {
  print(x="Up on sense versus Up on antisense:")
  print(x=summary(de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) > 0),1] %in% de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2]))
  print(x="Up on sense versus Down on antisense:")
  print(x=summary(de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) > 0),1] %in% de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2]))
  anticor <- intersect(x=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) > 0),1], y=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2])
  print(x=intersect(x=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) > 0),1], y=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2]))
  print(x="Down on sense versus Up on antisense:")
  print(x=summary(de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) < 0),1] %in% de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2]))
  anticor <- c(anticor, intersect(x=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) < 0),1], y=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2]))
  print(x=intersect(x=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) < 0),1], y=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2]))
  print(x="Down on sense versus Down on antisense:")
  print(x=summary(de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_sense$logFC", time, sep="_")), envir=globalenv()) < 0),1] %in% de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2]))
  anticor <- as.matrix(anticor)
  colnames(anticor) <- c("ensembl_gene_id")
  anticor <- merge(x=anticor, y=de_sense, by="ensembl_gene_id", all.x=TRUE)
  anticor <- merge(x=anticor, y=de_antisense, by.x="ensembl_gene_id", by.y="sense_ensembl_gene_id", all.x=TRUE)
  colnames(anticor) <- gsub(pattern=".x$", replacement="_S", x=colnames(anticor), perl=TRUE)
  colnames(anticor) <- gsub(pattern=".y$", replacement="_AS", x=colnames(anticor), perl=TRUE)
  return(anticor)
}

# Use the function to find common DEG at the different time points
anticor_2h <- sense_antisense_overlap(time="2H")
anticor_2h <- anticor_2h[,c(1:3,11,15,40,44)]
anticor_2h
anticor_6h <- sense_antisense_overlap(time="6H")
anticor_6h <- anticor_6h[,c(1:3,16,20,45,49)]
anticor_6h
anticor_24h <- sense_antisense_overlap(time="24H")
anticor_24h <- anticor_24h[,c(1:3,21,25,50,54)]
anticor_24h
anticor_48h <- sense_antisense_overlap(time="48H")
anticor_48h <- anticor_48h[,c(1:3,26,30,55,59)]
anticor_48h

# Output the results of anticorrelated sense and antisense gene expression
write.table(x=anticor_24h, file="Anticorrelated_sense_antisense.txt", sep="\t", append=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(x=anticor_48h, file="Anticorrelated_sense_antisense.txt", sep="\t", append=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE)

# Create a function to find unique antisense genes
antisense_unique <- function(time) {
  print(x="Up antisense and filtered out of sense:")
  print(x=summary(!(de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2] %in% de_sense[,1])))
  unique <- setdiff(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2], y=de_sense[,1])
  print(x=setdiff(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2], y=de_sense[,1]))
  print(x="Down antisense and filtered out of sense:")
  print(x=summary(!(de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2] %in% de_sense[,1])))
  unique <- c(unique, setdiff(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2], y=de_sense[,1]))
  print(x=setdiff(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2], y=de_sense[,1]))
  print(x="Up antisense and not DEG sense:")
  print(x=summary(de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2] %in% de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) > 0.05),1]))
  unique <- c(unique, intersect(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2], y=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) > 0.05),1]))
  print(x=intersect(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) > 0),2], y=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) > 0.05),1]))
  print(x="Down antisense and not DEG sense:")
  print(x=summary(de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2] %in% de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) > 0.05),1]))
  unique <- c(unique, intersect(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2], y=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) > 0.05),1]))
  print(x=intersect(x=de_antisense[(eval(expr=parse(text=paste("de_antisense$FDR", time, sep="_")), envir=globalenv()) < 0.05)&(eval(expr=parse(text=paste("de_antisense$logFC", time, sep="_")), envir=globalenv()) < 0),2], y=de_sense[(eval(expr=parse(text=paste("de_sense$FDR", time, sep="_")), envir=globalenv()) > 0.05),1]))
  unique <- as.matrix(unique)
  colnames(unique) <- c("ensembl_gene_id")
  unique <- merge(x=unique, y=de_antisense, by.x="ensembl_gene_id", by.y="sense_ensembl_gene_id", all.x=TRUE)
  unique <- merge(x=unique, y=de_sense, by="ensembl_gene_id", all.x=TRUE)
  colnames(unique) <- gsub(pattern=".x$", replacement="_AS", x=colnames(unique), perl=TRUE)
  colnames(unique) <- gsub(pattern=".y$", replacement="_S", x=colnames(unique), perl=TRUE)
  return(unique)
}

# Use the function to find common DEG at the different time points
unique_2h <- antisense_unique(time="2H")
unique_2h <- unique_2h[,c(1,3:4,11,15,40,44)]
unique_2h
unique_6h <- antisense_unique(time="6H")
unique_6h <- unique_6h[,c(1,3:4,16,20,45,49)]
unique_6h
unique_24h <- antisense_unique(time="24H")
unique_24h <- unique_24h[,c(1,3:4,21,25,50,54)]
unique_24h
unique_48h <- antisense_unique(time="48H")
unique_48h <- unique_48h[,c(1,3:4,26,30,55,59)]
unique_48h

# Output the results of unique antisense gene
write.table(x=unique_2h, file="Unique_antisense.txt", sep="\t", append=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(x=unique_6h, file="Unique_antisense.txt", sep="\t", append=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(x=unique_24h, file="Unique_antisense.txt", sep="\t", append=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE)
write.table(x=unique_48h, file="Unique_antisense.txt", sep="\t", append=TRUE, quote=FALSE, col.names=TRUE, row.names=FALSE)

# Create a matrix of antisense-sense genes comparison at each time point
comparison <- matrix(data=0, nrow=nrow(de_antisense), ncol=4)
rownames(comparison) <- de_antisense$bostaurus_gene_id
colnames(comparison) <- c("Comparison_2H", "Comparison_6H", "Comparison_24H", "Comparison_48H")
for (i in 1:nrow(de_antisense)) {
  if (de_antisense[i,"FDR_2H"] < 0.05) {
    if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == TRUE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_2H"] <- 2
      if ((de_antisense[i,"logFC_2H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_2H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_2H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_2H"] <- 3
      }
      else if ((de_antisense[i,"logFC_2H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_2H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_2H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_2H"] <- 4
      }
      else if ((de_antisense[i,"logFC_2H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_2H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_2H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_2H"] <- 5
      }
      else if ((de_antisense[i,"logFC_2H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_2H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_2H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_2H"] <- 6
      }
    }
    else if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == FALSE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_2H"] <- 1
    }
    else {
      stop("Value not correct at 2H!")
    }
  }
  if (de_antisense[i,"FDR_6H"] < 0.05) {
    if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == TRUE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_6H"] <- 2
      if ((de_antisense[i,"logFC_6H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_6H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_6H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_6H"] <- 3
      }
      else if ((de_antisense[i,"logFC_6H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_6H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_6H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_6H"] <- 4
      }
      else if ((de_antisense[i,"logFC_6H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_6H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_6H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_6H"] <- 5
      }
      else if ((de_antisense[i,"logFC_6H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_6H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_6H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_6H"] <- 6
      }
    }
    else if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == FALSE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_6H"] <- 1
    }
    else {
      stop("Value not correct at 6H!")
    }
  }
  if (de_antisense[i,"FDR_24H"] < 0.05) {
    if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == TRUE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_24H"] <- 2
      if ((de_antisense[i,"logFC_24H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_24H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_24H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_24H"] <- 3
      }
      else if ((de_antisense[i,"logFC_24H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_24H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_24H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_24H"] <- 4
      }
      else if ((de_antisense[i,"logFC_24H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_24H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_24H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_24H"] <- 5
      }
      else if ((de_antisense[i,"logFC_24H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_24H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_24H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_24H"] <- 6
      }
    }
    else if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == FALSE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_24H"] <- 1
    }
    else {
      stop("Value not correct at 24H!")
    }
  }
  if (de_antisense[i,"FDR_48H"] < 0.05) {
    if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == TRUE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_48H"] <- 2
      if ((de_antisense[i,"logFC_48H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_48H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_48H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_48H"] <- 3
      }
      else if ((de_antisense[i,"logFC_48H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_48H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_48H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_48H"] <- 4
      }
      else if ((de_antisense[i,"logFC_48H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_48H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_48H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_48H"] <- 5
      }
      else if ((de_antisense[i,"logFC_48H"] < 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"logFC_48H"] > 0) && (de_sense[grep(pattern=de_antisense[i,"sense_ensembl_gene_id"], x=de_sense[,"ensembl_gene_id"]),"FDR_48H"] < 0.05)) {
        comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_48H"] <- 6
      }
    }
    else if ((de_antisense[i,"sense_ensembl_gene_id"] %in% de_sense[,"ensembl_gene_id"]) == FALSE) {
      comparison[grep(pattern=de_antisense[i,"bostaurus_gene_id"], x=rownames(comparison)), "Comparison_48H"] <- 1
    }
    else {
      stop("Value not correct at 48H!")
    }
  }
}
head(comparison)
dim(comparison)

# Output the results of unique antisense gene
write.table(x=comparison, file="Comparison_AS_S.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=TRUE)

# Format the data for plotting
data <- merge(x=de_antisense[, c("sense_ensembl_gene_id", "logFC_2H", "FDR_2H", "logFC_6H", "FDR_6H", "logFC_24H", "FDR_24H", "logFC_48H", "FDR_48H")], y=de_sense[, c("ensembl_gene_id", "logFC_2H", "FDR_2H", "logFC_6H", "FDR_6H", "logFC_24H", "FDR_24H", "logFC_48H", "FDR_48H")], by.x="sense_ensembl_gene_id", by.y="ensembl_gene_id")
data <- cbind(data, c(rep(x="2", times=nrow(data))), c(rep(x="6", times=nrow(data))), c(rep(x="24", times=nrow(data))), c(rep(x="48", times=nrow(data))), c(rep(x="Antisense", times=nrow(data))), c(rep(x="Sense", times=nrow(data))))
head(data)
dim(data)
data_2_plot <- rbind(as.matrix(data[, c(1, 22, 18, 2:3)]), as.matrix(data[, c(1, 22, 19, 4:5)]), as.matrix(data[, c(1, 22, 20, 6:7)]), as.matrix(data[, c(1, 22, 21, 8:9)]), as.matrix(data[, c(1, 23, 18, 10:11)]), as.matrix(data[, c(1, 23, 19, 12:13)]), as.matrix(data[, c(1, 23, 20, 14:15)]), as.matrix(data[, c(1, 23, 21, 16:17)]))
colnames(data_2_plot) <- c("ensembl_gene_id", "Methods", "Time_in_Hours", "Log2FC", "FDR")
data_2_plot <- as.data.frame(data_2_plot)
data_2_plot[,3] <- as.numeric(as.character(x=data_2_plot[,3]))
data_2_plot[,4] <- as.numeric(as.character(x=data_2_plot[,4]))
data_2_plot[,5] <- as.numeric(as.character(x=data_2_plot[,5]))
head(data_2_plot)
dim(data_2_plot)

# Add significance label
sig_label <- function(arg1) {
  Significance_label <- vector()
  for (j in 1:length(arg1$FDR)) {
    if (arg1$FDR[j] < 0.001) {
      Significance_label <- c(Significance_label, "***")
    }
    else if (arg1$FDR[j] < 0.01) {
      Significance_label <- c(Significance_label, "**")
    }
    else if (arg1$FDR[j] < 0.05) {
      Significance_label <- c(Significance_label, "*")
    }
    else {
      Significance_label <- c(Significance_label, "")
    }
  }
  arg1 <- cbind(arg1, Significance_label)
  return(arg1)
}
data_2_plot <- sig_label(data_2_plot)
head(data_2_plot)

# Load required package for plotting
library(ggplot2)
library(grid)

# Plotting of the logFC for antisense and sense data for each gene
setwd("C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Overlap/plot_AS_S")
getwd()
not_plotted <- vector()
for (i in 1:nrow(data)) {
  if (((data[i, "FDR_2H.x"] < 0.05) == FALSE) && ((data[i, "FDR_6H.x"] < 0.05) == FALSE) && ((data[i, "FDR_24H.x"] < 0.05) == FALSE) && ((data[i, "FDR_48H.x"] < 0.05) == FALSE)) {
    not_plotted <- c(not_plotted, as.character(data[i, "sense_ensembl_gene_id"]))
  }
  else {
    plot1 <- qplot(data=data_2_plot[grep(pattern=as.character(data[i, "sense_ensembl_gene_id"]), x=data_2_plot$ensembl_gene_id),], x=Time_in_Hours, y=Log2FC, colour=Methods, main=as.character(data[i, "sense_ensembl_gene_id"]))+geom_point(shape=19, size=15)+geom_line(size=5)+geom_text(aes(x=Time_in_Hours+2, y=Log2FC, label=Significance_label), size=22)+theme(panel.background=element_rect(fill='wheat'), legend.title=element_text(size=45, face="bold"), legend.text=element_text(size=45, face="bold"), axis.title.x=element_text(face="bold", size=45), axis.text.x=element_text(face="bold", size=45), axis.title.y=element_text(face="bold", size=45), axis.text.y=element_text(face="bold", size=45), plot.title=element_text(face="bold", size=55), legend.position=c(0.33, 0.08))+xlab("Time (hours)")+ylab("Log2 fold-change")+scale_color_manual(values=c(rgb(red=102, green=146, blue=62, maxColorValue=255), rgb(red=138, green=116, blue=168, maxColorValue=255)))
    png(filename=paste(as.character(data[i, "sense_ensembl_gene_id"]), ".png", sep=""), width=1366, height=1366, units="px")
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(1, 1, heights = unit(c(4), "null"))))
    print(plot1, vp=viewport(layout.pos.row=1, layout.pos.col=1), )
    dev.off()
  }
}
not_plotted
length(not_plotted)

# Reset the working directory
setwd(dir="C:/Users/nnalpas/Documents/PhD project/Alveolar macrophages/RNA-seq analysis/Results/edgeR/Analysis 251113/Overlap")
getwd()

#######
# END #
#######