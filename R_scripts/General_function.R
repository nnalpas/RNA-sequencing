#################################################
# Define functions used accross several scripts #
#################################################

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
    source(file="http://bioconductor.org/biocLite.R", verbose=FALSE)
    biocLite(pkgs=deparse(substitute(package)), suppressUpdates=TRUE)
    if(require(package=deparse(substitute(package)), character.only=TRUE, quietly=TRUE)) {
      print(paste(deparse(substitute(package)), " is correctly installed and loaded!", sep=""))
    }
    else {
      stop(paste('"', "Could not install ", deparse(substitute(package)), '"', sep=""))
    }
  }
  print(paste(deparse(substitute(package)), " version: ", packageVersion(pkg=deparse(substitute(package))), sep=""))
}

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
  print("Summary of the number of edgeR DEG:")
  print(summary(decideTestsDGE(lrt, p.value=adjpvalue)))
  print("Names of the edgeR multiple correction test dataframe:")
  print(names(de))
  print("Heading of the edgeR multiple correction test dataframe:")
  print(head(de$table))
  if(length(grep(pattern="antisense", x=DEfile))==1) {
    write.table(x=de$table[,c("sense_ensembl_gene_id", "external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  }
  else if(length(grep(pattern="novel", x=DEfile))==1) {
    write.table(x=de$table[,c("Hsapiens_ensembl_gene_id", "external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  }
  else {
    write.table(x=de$table[,c("external_gene_id","description","logFC", "logCPM", "LR", "PValue", "FDR")], file=paste(DEfile, "txt", sep="."), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)
  }
  png(filename=paste(Smearfile, "png", sep="."), width=1366, height=768, units="px")
  plotSmear(object=lrt, de.tags=(rownames(lrt$table)[as.logical(decideTestsDGE(lrt, p.value = 0.05))]))
  abline(h=c(-1, 1), col="blue")
  dev.off()
  assign(x=LRTdata, value=lrt, envir=.GlobalEnv)
  assign(x=DEdata, value=de, envir=.GlobalEnv)
}

################################
# Significance label for graph #
################################

# Create a function to add significance label
sig_label <- function(arg1, arg2) {
  Significance_label <- vector()
  for (j in 1:length(arg1[,arg2])) {
    if (arg1[j,arg2] < 0.001) {
      Significance_label <- c(Significance_label, "***")
    }
    else if (arg1[j,arg2] < 0.01) {
      Significance_label <- c(Significance_label, "**")
    }
    else if (arg1[j,arg2] < 0.05) {
      Significance_label <- c(Significance_label, "*")
    }
    else {
      Significance_label <- c(Significance_label, "")
    }
  }
  arg1 <- cbind(arg1, Significance_label)
  return(arg1)
}

#######
# END #
#######