#DESeq2 Pipeline

#Install the DESeq2 Package (Only need to be run once)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
#only run once

#Step 1: Load the DESeq2package
library(DESeq2)

#Step 2: Count File Generation

#Load your files in 
file.list <- list.files( path = "./", pattern = "*ReadsPerGene.out.tab$")

#Form your combined counts table
counts.files <- lapply(file.list, read.table, skip = 4)

#Select your counts column
counts <- as.data.frame( sapply( counts.files, function(x) x[ ,2] ) )
View(counts)

#Setting Column and Row Names for your Data Frame
colnames(counts) <- file.list
row.names(counts) <- counts.files[[1]]$V1
View(counts)

#Step 3: Conditions and Metadata

#Setting Conditions
condition <- c(rep("Control",5), rep("DM1",5))

#Above code may not work sometimes, instaed use condition <- c(“Control”,”DM1”,”Control”... etc.) for all 10 samples in order to orient them correctly

#Meta Data Matrix Formation
sampleTable <- data.frame(sampleName = file.list, condition = condition)
View(sampleTable)

#DESeq Dataset Formation
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sampleTable, design = ~ condition)

#Step 4: Running and Writing Data

#Running DESeq
output <- DESeq(dds)

#Obtain your Results
results_Control_DM <- results(output, contrast=c("condition","DM1","Control"))

#Sort Results by Adjusted P-Value
results_Control_DM_PValue <- results_Control_DM[order(results_Control_DM$padj),]

#Head Your Results and Write Them to File
head(results_Control_DM_PValue)
write.csv(as.data.frame(results_Control_DM_PValue), file="Control_DM_DE.csv")
#This heads the results that you have just analyzed, and then saves the results to a CSV file that you can look at in a text editor or in Excel

#When you run differential gene expression analysis using DESeq2, the genes that are differentially expressed will be given as ensembl IDs instead of gene IDs
#It is much more useful to convert the ensemble IDs to gene IDs, which can be easily be done directly in Rstudio or online using the website biotools.fr

#Converting ENSEMBL IDs to Gene IDs
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library(AnnotationDbi)
library(org.Hs.eg.db)
DM1_Control <- read.table("Control_DM_DE.csv", header = TRUE, sep =',')
IDs <- c(DM1_Control$X)
DM1_Control$Symbol <- mapIds(org.Hs.eg.db, IDs, 'SYMBOL', 'ENSEMBL')
rownames(DM1_Control) <- DM1_Control$X

#This code can be used to write a new file with the GeneIDs added
write.csv(DM1_Control, file = “Control_DM_DE_GeneIDs.csv”)

#Done
