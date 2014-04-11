########################################################
# Analysis of alveolar macrophages flow cytometry data #
########################################################

#############################
# List of required packages #
#############################

# Source the common functions used across this script
source(file="F:/nnalpas/Documents/PhD project/Bioinformatics/R/General_function.R")

# Load the required packages
loadpackage(package=flowCore)
loadpackage(package=flowViz)

###############################
# Prepare working environment #
###############################

# Set the path to working directory
setwd(dir="F:/nnalpas/Documents/PhD project/Alveolar macrophages/MB-TB infection (Mar-Jun2012)/Flow cytometry/R analysis")
getwd()

####################################
# Read in all flow cytometry files #
####################################

# Create a vector of all fcs files
files <- list.files(path="F:/nnalpas/Documents/PhD project/Alveolar macrophages/MB-TB infection (Mar-Jun2012)/Flow cytometry/R analysis/files")
files <- grep(x=files, pattern="N.*((St)|(Se)|(Un)).*.fcs", perl=TRUE, value=TRUE)

# Read all the files into a flowframe
alv_mac <- read.flowSet(files=files, path="F:/nnalpas/Documents/PhD project/Alveolar macrophages/MB-TB infection (Mar-Jun2012)/Flow cytometry/R analysis/files")
pData(phenoData(alv_mac))

# Check data
alv_mac
summary(alv_mac)

#####################################
# Create gate1 FS log versus SS log #
#####################################

# Create gate 1 for FS log versus SS log
gate1 <- rectangleGate(filterId="gate1", `FS Log`=c(750,8000), `SS Log`=c(75,7800))

# Apply gate to the flowset
alv_mac_gate1 <- Subset(x=alv_mac, subset=gate1)

# Check the data
summary(filter(x=alv_mac, filter=gate1))
capture.output(summary(filter(x=alv_mac, filter=gate1)), file="Gate1.txt", append=TRUE)

######################################
# Create gate2 FS area versus FS lin #
######################################

# Create gate 2 for FS area versus FS lin
gate2 <- matrix(c(17000,60000,58000,17000,9000,33000,43000,16000),ncol=2,nrow=4)
colnames(gate2) <- c("FS Area","FS Lin")
gate2 <- polygonGate(filterId="gate2", .gate=gate2)

# Apply gate to the flowset
alv_mac_gate2 <- Subset(x=alv_mac_gate1, subset=gate2)

# Check the data
summary(filter(x=alv_mac_gate1, filter=gate2))
capture.output(summary(filter(x=alv_mac_gate1, filter=gate2)), file="Gate2.txt", append=TRUE)

###########################################
# Create gate3 SS area versus Pulse width #
###########################################

# Create gate 3 for SS area versus Pulse width
gate3 <- matrix(c(12000,59000,60000,25000,10000,40,70,110,100,80),ncol=2,nrow=5)
colnames(gate3) <- c("SS Area","Pulse Width")
gate3 <- polygonGate(filterId="gate3", .gate=gate3)

# Apply gate to the flowset
alv_mac_gate3 <- Subset(x=alv_mac_gate2, subset=gate3)

# Check the data
summary(filter(x=alv_mac_gate2, filter=gate3))
capture.output(summary(filter(x=alv_mac_gate2, filter=gate3)), file="Gate3.txt", append=TRUE)

###########################################
# Create gate4 FS area versus Pulse width #
###########################################

# Create gate 4 for FS area versus Pulse width
gate4 <- matrix(c(18000,61000,61000,31000,17000,40,70,125,115,90),ncol=2,nrow=5)
colnames(gate4) <- c("FS Area","Pulse Width")
gate4 <- polygonGate(filterId="gate4", .gate=gate4)

# Apply gate to the flowset
alv_mac_gate4 <- Subset(x=alv_mac_gate3, subset=gate4)

# Check the data
summary(filter(x=alv_mac_gate3, filter=gate4))
capture.output(summary(filter(x=alv_mac_gate3, filter=gate4)), file="Gate4.txt", append=TRUE)

#######################################
# Create gate5 FL 1 Log (alias CD14+) #
#######################################

# Create gate 5 for CD14 log versus Count
gate5 <- rectangleGate(filterId="gate5", `FL 1 Log`=c(17,10000))

# Obtain summary for CD14+ cells
summary(filter(x=alv_mac_gate4, filter=gate5))
capture.output(summary(filter(x=alv_mac_gate4, filter=gate5)), file="Gate5.txt", append=TRUE)

#####################################
# Create gate6 FL 2 Log (alias PE+) #
#####################################

# Create gate 6 for PE log versus Count
gate6 <- rectangleGate(filterId="gate6", `FL 2 Log`=c(17,10000))

# Obtain summary for PE+ cells
summary(filter(x=alv_mac_gate4, filter=gate6))
capture.output(summary(filter(x=alv_mac_gate4, filter=gate6)), file="Gate6.txt", append=TRUE)

##############################################
# Plot all gates into single plot per sample #
##############################################

# Plot the data for gate1,2,3,4,5,6
for (i in 1:36) {
  file_name <- gsub(pattern=".fcs$", replacement="_flow_cytometry.png", x=sampleNames(alv_mac)[i], perl=TRUE)
  png(filename=file_name, width=1000, height=1500, units="px")
  par(mfrow=c(3,2), cex=2.0)
  plot(transform(alv_mac[[i]], FS_log=(`FS Log`), SS_log=(`SS Log`)), c("FS_log", "SS_log"), log="xy", smooth=FALSE, col=rgb(240,34,34,150,maxColorValue=255))
  rect(750, 75, 8000, 7800)
  text(x=250, y=4900, labels="gate1")
  plot(transform(alv_mac_gate1[[i]], FS_area=(`FS Area`), FS_lin=(`FS Lin`)), c("FS_area", "FS_lin"), smooth=FALSE, col=rgb(240,34,34,150,maxColorValue=255))
  polygon(x=c(17000,60000,58000,17000), y=c(9000,33000,43000,16000))
  text(x=14000, y=21500, labels="gate2")
  plot(transform(alv_mac_gate2[[i]], SS_area=(`SS Area`), Pulse_width=(`Pulse Width`)), c("SS_area", "Pulse_width"), smooth=FALSE, col=rgb(240,34,34,150,maxColorValue=255))
  polygon(x=c(12000,59000,60000,25000,10000), y=c(40,70,110,100,80))
  text(x=9000, y=120, labels="gate3")
  plot(transform(alv_mac_gate3[[i]], FS_area=(`FS Area`), Pulse_width=(`Pulse Width`)), c("FS_area", "Pulse_width"), smooth=FALSE, col=rgb(240,34,34,150,maxColorValue=255))
  polygon(x=c(18000,61000,61000,31000,17000), y=c(40,70,125,115,90))
  text(x=14000, y=120, labels="gate4")
  plot(transform(alv_mac_gate4[[i]], CD14_log=(`FL 1 Log`), Event_Count=sqrt(`Event Count`)), c("CD14_log", "Event_Count"), log="x", smooth=FALSE, col=rgb(240,34,34,150,maxColorValue=255))
  rect(17, 0.01, 10000, 256)
  text(x=50, y=225, labels="gate5")
  plot(transform(alv_mac_gate4[[i]], PE_log=(`FL 2 Log`), Event_Count=sqrt(`Event Count`)), c("PE_log", "Event_Count"), log="x", smooth=FALSE, col=rgb(240,34,34,150,maxColorValue=255))
  rect(17, 0.01, 10000, 256)
  text(x=50, y=225, labels="gate6")
  mtext(text=file_name, side=3, line=-2, outer=TRUE, cex=2.5)
  dev.off()
}

#######
# END #
#######