#################################################################
# Visualisation of KEGG pathway with overlaying gene expression #
#################################################################

# Once the pathway of interest have been determined, use KEGGprofile package to overlay gene expression across several time points on the pathway

#############################
# List of required packages #
#############################

# Source the common functions used across this script
source(file="C:/Users/nnalpas/Dropbox/Home_work_sync/Bioinformatics/R/General_function.R")

# Load the required packages
loadpackage(package=KEGGprofile)
loadpackage(package=gdata)
loadpackage(XML)
loadpackage(RCurl)

#####################################################################
# Adapt the KEGGprofile package functions fort better visualisation #
#####################################################################

# Define the plot_polygon function used within KEGGprofile package
plot_polygon<-function(ChosedGene,i,i_max,col,height,err_x_location,err_y_location,magnify,border=F,...) {
  polygon(c(ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*(i-1)*magnify+err_x_location,ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*i*magnify+err_x_location,ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*i*magnify+err_x_location,ChosedGene[,2]-ChosedGene[,4]/2*magnify+ChosedGene[,4]/i_max*(i-1)*magnify+err_x_location), c(height-ChosedGene[,3]-ChosedGene[,5]/2*magnify-err_y_location,height-ChosedGene[,3]-ChosedGene[,5]/2*magnify-err_y_location,height-ChosedGene[,3]+ChosedGene[,5]/2*magnify-err_y_location,height-ChosedGene[,3]+ChosedGene[,5]/2*magnify-err_y_location), col=col[as.character(ChosedGene[,1]),i],border=border,...)
}

# Edit the plot_profile function from KEGGprofile package to export result picture as high quality TIFF file
new_plot_profile <- function (gene_expr, pathway_name, KEGG_database, groups, bg_col = "white", text_col = "black", line_col, border_col = "grey", text_cex = 0.25, magnify = 1, type = c("lines", "bg"), pathway_min = 5, genes_kept = c("foldchange", "first", "random", "var", "abs"), specis = "hsa", database_dir = getwd(), max_dist, lwd = 1.2) {
  type <- if (missing(type)) 
    "lines"
  else match.arg(type)
  if (type == "lines" & ncol(gene_expr) <= 1) {
    print("When type=='lines', You should have more than one time points")
  }
  genes_kept <- if (missing(genes_kept)) 
    "foldchange"
  else match.arg(genes_kept)
  if (missing(groups)) 
    groups <- rep(1, ncol(gene_expr))
  groups <- factor(groups, levels = unique(groups))
  if (missing(line_col)) 
    line_col <- rainbow(length(unique(groups)))
  if (is.matrix(bg_col) | is.data.frame(bg_col)) {
  }
  else {
    bg_col <- matrix(rep(bg_col, nrow(gene_expr)), ncol = 1)
    row.names(bg_col) <- row.names(gene_expr)
  }
  if (is.matrix(border_col) | is.data.frame(border_col)) {
  }
  else {
    border_col <- matrix(rep(border_col, nrow(gene_expr)), ncol = 1)
    row.names(border_col) <- row.names(gene_expr)
  }
  if (is.matrix(text_col) | is.data.frame(text_col)) {
  }
  else {
    text_col <- matrix(rep(text_col, nrow(gene_expr)), ncol = 1)
    row.names(text_col) <- row.names(gene_expr)
  }
  return_expr <- NULL
  genes <- intersect(KEGG_database[, 1], row.names(gene_expr))
  if (length(genes) < pathway_min) {
    print(paste("The genes mapped in pathway ", pathway_name, " were less than ", pathway_min, ", skip this pathway.", sep = ""))
    return()
  }
  require(png)
  require(TeachingDemos)
  img <- readPNG(paste(database_dir, "/", pathway_name, ".png", sep = ""))
  width <- ncol(img)
  height <- nrow(img)
  err_x_location <- 1
  err_y_location <- 0
  tiff(paste(pathway_name, "_profile_", type, ".tiff", sep = ""), width = width * 2, height = height * 2, res = 600)
  par(yaxs = "i")
  par(xaxs = "i")
  par(mar = c(0, 0, 0, 0))
  plot(c(0, width), c(0, height), type = "n", xlab = "", ylab = "")
  rasterImage(img, 0, 0, width, height, interpolate = F)
  result_genes <- as.data.frame(KEGG_database[which(KEGG_database[, 1] %in% genes), ], stringsAsFactors = F)
  colnames(result_genes) <- c("genes", "x", "y", "width", "height", "name")
  result_genes <- transform(result_genes, x = as.numeric(x), y = as.numeric(y), width = as.numeric(width), height = as.numeric(height))
  if (missing(max_dist) & type == "lines") {
    max_dist <- max(apply(matrix(gene_expr[result_genes[, 1], ]), 1, function(x) range(x, na.rm = T)[2] - range(x, na.rm = T)[1]))
  }
  findUnique <- apply(result_genes[, 2:3], 1, function(x) paste(x, collapse = " "))
  temp <- split(as.data.frame(result_genes, stringsAsFactors = F), findUnique)
  for (xx in 1:length(temp)) {
    if (length(temp[[xx]][, 1]) > 1) {
      if (genes_kept == "foldchange") {
        ChosedGene <- temp[[xx]][which.max(apply(data.frame(gene_expr[temp[[xx]][, 1], ]), 1, function(x) range(x, na.rm = T)[2] - range(x, na.rm = T)[1])), 1]
      }
      else if (genes_kept == "first") {
        ChosedGene <- temp[[xx]][1, 1]
      }
      else if (genes_kept == "random") {
        ChosedGene <- sample(temp[[xx]][, 1], 1)
      }
      else if (genes_kept == "var") {
        ChosedGene <- temp[[xx]][which.max(apply(data.frame(gene_expr[temp[[xx]][, 1], ]), 1, var)), 1]
      }
      else if (genes_kept == "abs") {
        ChosedGene <- temp[[xx]][which.max(apply(data.frame(gene_expr[temp[[xx]][, 1], ]), 1, function(x) max(abs(x)))), 1]
      }
      ChosedGene <- temp[[xx]][which(temp[[xx]][, 1] == ChosedGene), ]
    }
    else {
      ChosedGene <- temp[[xx]]
    }
    i_max <- ncol(bg_col)
    for (i in 1:i_max) {
      plot_polygon(ChosedGene = ChosedGene, i = i, i_max = i_max, col = bg_col, height = height, err_x_location = err_x_location, err_y_location = err_y_location, magnify = magnify)
    }
    if (type == "lines") {
      ChosedGeneProfile <- as.matrix(gene_expr[as.character(ChosedGene[, 1]), ])
      ChosedGeneProfile <- sapply(split(ChosedGeneProfile, groups), function(x) x)
      gene_dist <- gene_expr[as.character(ChosedGene[, 1]), ]
      gene_dist <- range(gene_dist, na.rm = T)[2] - range(gene_dist, na.rm = T)[1]
      if (is.na(max_dist) | gene_dist == 0) {
        y_ratio <- 1
      }
      else {
        y_ratio <- gene_dist/max_dist
        if (y_ratio > 1) {
          y_ratio <- 1
        }
      }
      old_par <- par(no.readonly = TRUE)
      subplot(matplot(ChosedGeneProfile, main = "", type = "l", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", col = line_col, lty = 1, lwd = lwd), c(ChosedGene[, 2] - ChosedGene[, 4]/2 * magnify + err_x_location, ChosedGene[, 2] + ChosedGene[, 4] * magnify/2 + err_x_location), c(height - ChosedGene[, 3] - ChosedGene[, 5] * y_ratio/2 * magnify - err_y_location, height - ChosedGene[, 3] + ChosedGene[, 5] * y_ratio/2 * magnify - err_y_location))
      par(old_par)
    }
    plot_polygon(ChosedGene = ChosedGene, col = as.data.frame(NULL), i = 1, i_max = 1, height = height, err_x_location = err_x_location, err_y_location = err_y_location, magnify = magnify, border = border_col[as.character(ChosedGene[, 1]), 1], lty = 1, lwd = 0.5)
    text(ChosedGene[, 2], height - ChosedGene[, 3], labels = ChosedGene[, 6], cex = text_cex, col = text_col[as.character(ChosedGene[, 1]), 1])
    return_expr <- rbind(return_expr, c(ChosedGene[, 1], ChosedGene[, 6], gene_expr[as.character(ChosedGene[, 1]), ]))
  }
  if (type == "lines") {
    legend("topright", legend = unique(groups), lwd = 1, col = line_col, bty = "n", cex = 0.3)
    if (!(is.na(max_dist))) {
      polygon(c(width - 66, width - 20, width - 20, width - 66), c(height - 60, height - 60, height - 77, height - 77), border = "grey")
      text(width - 43, height - 68.5, round(max_dist, 2), cex = text_cex)
    }
  }
  dev.off()
  print(paste("The figure was generated in ", getwd(), "/", pathway_name, "_profile_", type, ".png", sep = ""))
  colnames(return_expr) <- c("Gene", "Name", colnames(gene_expr))
  return(unique(return_expr))
}

# Edit the plot_pathway function from KEGGprofile package to remove parsing of XML file (since that was already done) and also use the new adapted function
new_plot_pathway <- function (gene_expr, line_col, groups, pathway_id = "00010", specis = "hsa", pathway_min = 5, XML2data, database_dir = getwd(), ...) {
  if ((!file.exists(paste(database_dir, "/", specis, pathway_id, ".xml", sep = ""))) | (!file.exists(paste(database_dir, "/", specis, pathway_id, ".png", sep = "")))) {
    download_KEGGfile(pathway_id = pathway_id, specis = specis, target_dir = database_dir)
  }
  if (is.null(XML2data)) {
    return()
  }
  return_expr <- new_plot_profile(gene_expr = gene_expr, KEGG_database = XML2data, groups = groups, line_col = line_col, pathway_name = paste(specis, pathway_id, sep = ""), database_dir = database_dir, ...)
  return(return_expr)
}

###########################################################################
# Visualise the KEGG pathway of interest: RIG-I-like receptors signalling #
###########################################################################

# Download the pathway of interest (RIG-I-like receptors signalling)
download_KEGGfile(pathway_id="04622", specis="bta")

# Read in the table of gene expression used for RIG-I pathway (containing all genes within this pathway)
rig1 <- read.table(file="RIG1_gene.txt", sep="\t", header=TRUE, fill=TRUE)
rownames(rig1) <- rig1$entrezID
rig1 <- rig1[,-1]
head(rig1)

# Prepare the color code matrix for each expression value (without the NA value for gene expression)
temp <- apply(rig1, 1, function(x) length(which(is.na(x))))
tiff(filename="color.tiff", width=2000, height=2000, res=600)
col <- col_by_value(rig1[which(temp==0), 2:5], col=colorRampPalette(c("darkblue", "blue", "lightgrey", "red", "darkred"))(1025), range=c(-8, 8))
dev.off()

# Modify the color code matrix to include gene with NA expression value and replace their NA value by a white color code
colnames(col) <- c("logFC2h", "logFC6h", "logFC24h", "logFC48h")
col <- rbind(col, rig1[which(temp==4), -1])
col <- as.matrix(col)
col[is.na(col)] <- "#FFFFFF"
colnames(col) <- NULL

# Edit the XML file manually if gene are missing then parse it to obtain the plotting information and had official gene name to the resulting matrix
test_XML <- parse_XMLfile(pathway_id="04622", specis="bta")
for (x in 1:dim(test_XML)[1]){
  test_XML[x, 6] <- as.character(rig1[test_XML[x,1], "Name"])
}
rig1 <- rig1[,-1]

# Graph and obtain picture of RIG-I pathway after customisation of functions from KEGGprofile package
temp <- new_plot_pathway(rig1, XML2data=test_XML, type="bg", bg_col=col, text_col="black", magnify=1.0, text_cex=0.25, specis="bta", database_dir=system.file("extdata", package="KEGGprofile"), pathway_id="04622")

# Clean up variables
rm(temp, rig1, col, test_XML)

####################################################
# Visualise the KEGG pathway of interest: Lysosome #
####################################################

# Download the pathway of interest (Lysosome)
download_KEGGfile(pathway_id="04142", specis="bta")

# Read in the table of gene expression used for lysosome pathway (containing all genes within this pathway)
lysosome <- read.table(file="Lysosome_gene.txt", sep="\t", header=TRUE, fill=TRUE)
rownames(lysosome) <- lysosome$entrezID
lysosome <- lysosome[,-1]
head(lysosome)

# Prepare the color code matrix for each expression value (without the NA value for gene expression)
temp <- apply(lysosome, 1, function(x) length(which(is.na(x))))
col <- col_by_value(lysosome[which(temp==0), 2:5], col=colorRampPalette(c("darkblue", "blue", "lightgrey", "red", "darkred"))(1025), range=c(-8, 8))

# Modify the color code matrix to include gene with NA expression value and replace their NA value by a white color code
colnames(col) <- c("logFC2h", "logFC6h", "logFC24h", "logFC48h")
col <- rbind(col, lysosome[which(temp==4), -1])
col <- as.matrix(col)
col[is.na(col)] <- "#FFFFFF"
colnames(col) <- NULL

# Edit the XML file manually if gene are missing then parse it to obtain the plotting information and had official gene name to the resulting matrix
test_XML <- parse_XMLfile(pathway_id="04142", specis="bta")
for (x in 1:dim(test_XML)[1]){
  test_XML[x, 6] <- as.character(lysosome[test_XML[x,1], "Name"])
}
lysosome <- lysosome[,-1]

# Graph and obtain picture of Lysosome pathway after customisation of functions from KEGGprofile package
temp <- new_plot_pathway(lysosome, XML2data=test_XML, type="bg", bg_col=col, text_col="black", magnify=1.0, text_cex=0.25, specis="bta", database_dir=system.file("extdata", package="KEGGprofile"), pathway_id="04142")

# Clean up variables
rm(temp, lysosome, col, test_XML)

#######
# END #
#######