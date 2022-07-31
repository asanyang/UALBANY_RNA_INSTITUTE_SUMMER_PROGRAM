#First step: downloading DESeq 2

#Install the DESeq2 Package (This will only need to be run one time)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#Only will have to download Bioconductor once

BiocManager::install("DESeq2")
#Only required to run one time

#Load the DESeq2 package 
library(DESeq2)

#Second Step: Count File generation

#Load your files in 
file.list <- list.files( path = "./", pattern = "*ReadsPerGene.out.tab$")

#Form your combined counts table
counts.files <- lapply(file.list, read.table, skip = 4)

#Select your counts column
counts <- as.data.frame( sapply( counts.files, function(x) x[ ,2] ) )

#Setting Column and Row Names for your Data Frame
colnames(counts) <- file.list
row.names(counts) <- counts.files[[1]]$V1

#Third Step: Conditions and Metadata

#Setting Conditions
condition <- c(rep("Control",5), rep("DM1",5))
 condition <- c(“Control”,”DM1”,”Control”... etc.) for all 10 samples in order to orient them correctly
 
#Meta Data Matrix Formation
sampleTable <- data.frame(sampleName = file.list, condition = condition)

#DESeq Dataset Formation
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)

