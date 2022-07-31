#Figure generation from differential gene expression analysis

#After running DESeq2, information on the differential expression of 60683 genes between DM1 and control has been obtained
#Next is to extract meaningful information from these large tables to be able to answer our research questions
#Data is viewed in different ways acc. to whatever research questions that is to be addressed (Data Visualization)

#Question 1: Do your samples cluster by condition? 
#A PCA Plot to plot each sample separately to check your control samples are all similar to each other but different to your DM1 samples
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("condition"))

#Question 2: Do you have a large number of genes that are differentially expressed? Are differentially expressed genes primarily upregulated, downregulated or evenly split between up and down regulation? 
#A Volcano plot or an MA plot can be used to plot each gene in your dataset to look at overall trends of differential gene expression
#Volcano Plot
#Install EnhancedVolcano package and load it
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
EnhancedVolcano(results_Control_DM,
                lab = DM1_Control$Symbol,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.05,
                FCcutoff = 1.5,
                pointSize = 1.5,
                labSize = 3.0,
               title = "DM1 vs Control",
               xlim = c(-10,10),
               ylim = c(0,10))
#Note: Adjust axis limits (xlim and ylim) accordingly to capture all points and not cut off any

#MA Plot
plotMA(results_Control_DM, ylim=c(-10,10))
#results_Control_DM is what the results from DESeq2 is named, this may need to be adjusted on the y coordinates for a better image

#Question 3: Do the top differentially expressed genes behave in the same way for all DM1 samples? 
#A Heat map can plot the most differentially expressed genes to look at changes between samples
#Install and load Pheatmap package
install.packages("pheatmap")
library("pheatmap")
rld <- rlog(dds, blind=FALSE)
topgenes <- head(rownames(results_Control_DM_PValue),20)
topgenes_symbols <- DM1_Control[topgenes,]
mat <- assay(rld)[topgenes,]

mat <- mat - rowMeans(mat)
rownames(mat) <- topgenes_symbols$Symbol
df <- as.data.frame(colData(dds)[,c("condition")])
colnames(df) <- "Condition"
row.names(df) <- sampleTable$sampleName
pheatmap(mat, annotation_col=df, cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=FALSE)

#Question 4: How different is your specific gene of interest between control and DM1? 
#Plotting normalized counts can plot specific genes of interest individually to show changes between control and DM1

plotCounts(dds, gene=which.min(results_Control_DM$padj), intgroup="condition")

#This will plot the normalized read counts for the most statistically significant gene based on adjusted p-value
#The gene can equal to any gene of interest in ENSG format (Example_gene=“ENSG00000286235”)
#Adding returnData = TRUE will give the counts for each individual library, useful for identifying outliers

#Tip to save plots in R
#tiff('MyGraph.tiff', width = 7, height = 7, units = 'in', res=300)
#*Your Graphing Code* i.e. plotMA(results_Control_DM, ylim=c(-10,10))
#dev.off()
