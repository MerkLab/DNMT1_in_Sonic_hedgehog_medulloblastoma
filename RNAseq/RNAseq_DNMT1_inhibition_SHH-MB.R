getwd()
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Foteini/R_DESeq2/Test/Data")
library(DESeq2)
# cannot load pheatmap although package installed apparently, load heatmap instead 
# read what the fehler says, execute it (delete some folder from documents bla bla) and try to install again, works!
install.packages("pheatmap")
library(pheatmap)
install.packages("ggplot2")
library(ggplot2)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("heatmaps")

# sep function separates columns by semicolon, comma is the default and recognizes as 1 variable instead of 13
data <-  read.csv("rawdata.csv", header = T, sep = ";")

# structure, chr is a character and int integer aka akeraios arithmos
str(data)

# aggregate sums up duplicate genes 
data.agg <- aggregate(data[-1], by = list(data$Gene.Symbol), FUN=sum)
# generate an additional column (Group1 with the gene names)  
row.names(data.agg) <- data.agg$Group.1
# kick out group1 [-1] and leave gene symbols alone, not part of a group
data.agg <- data.agg[, -1]

setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Foteini/R_DESeq2/DNMT1/Data")
# open metadata, row=1 keeps only 1 column 
meta <- read.csv("meta.csv", header=T, sep = ";", row=1)
# class checks data frames
class(data.agg)
class(meta)
# all checks whether data names match and names in order (true or false)
all(names(data.agg) %in% rownames(meta))
all(names(data.agg) == rownames(meta))
# dds is the Object, design refers to the column of the metadata, run deseq library, had an error before
dds <- DESeqDataSetFromMatrix(countData = data.agg, colData = meta, design = ~ sampletype)
# still not sure what sizefactor is, refers to gene length and seq depth
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
# nest within view function so we see it in script editor
View(counts(dds))

# filter out low/non expressed genes 
# dim gives the number of elements (rows, columns), 53.112 rows, 12 columns
dim(dds)

## genes out where less than 12 samples have norm counts greater/equal to 10, greater/equal to 12 samples
idx. <- rowSums(counts(dds, normalized=TRUE) >= 10) >= 12
dds. <- dds[idx.,]
# eventually 13,743 rows 
dim(dds.)
# retrieve normalized counts from dds 
normalized_counts. <- counts(dds., normalized=TRUE)
dim(normalized_counts.)
write.table(normalized_counts, file="normalized_counts_.csv", sep="\t", quote=F, col.names=NA)

## from now on, work with dds. and idx.
#QC,Regularized Log Transform of norm counts, intgroup is the column of metadata we're interested in  
rld <- rlog(dds., blind=TRUE)
plotPCA(rld, ntop=500, intgroup=("sampletype"))
rld.aza.only <- rld[,c(1,2,3,4,5,6,7,8,9)]
plotPCA(rld.aza.only, ntop=500, intgroup=("sampletype"))

rld_mat <- assay(rld) 
write.table(rld_mat, file="rlog_values_.txt", sep="\t", quote=F, col.names=NA)
head(colData())
# compute pairwise correlation values, correlation bar is missing 
rld_cor <- cor(rld_mat)
pheatmap(rld_cor)

# important function!
dds. <- DESeq(dds.)

# dispersion estimates reflect the variance in gene expression for a given mean value
getwd()
plotDispEsts(dds.)
# of raw reads per sample 
colSums(counts(dds.))
# of normalized counts per sample
colSums(counts(dds., normalized=T))

# Identify gene expression changes using pairwise Wald tests statistics (pairwise comparisons)
# coefficients 
resultsNames(dds.)

# specify which contrast you have for diff expression (both commands should be right), alpha refers to FDR
contrast_1. <- list( "sampletypeAZA_2h", "sampletypeDMSO")
contrast_1. <- c("sampletype", "AZA_2h", "DMSO")
res_table_1. <- results(dds., contrast=contrast_1., alpha = 0.05)
head(res_table_1.)
head(res_table_1)
contrast_2. <- c("sampletype", "AZA_48h", "DMSO")
res_table_2. <- results(dds., contrast=contrast_2., alpha = 0.05)

# I can remove tables and so on
rm(res_table_1_unshrunken)

# whats this table 
mcols(res_table_1, use.names=T)

# save in dataframe
res_1_df. <- data.frame(res_table_1.)
res_2_df. <- data.frame(res_table_2.)

# Summarizing results
summary(res_table_1., alpha=0.05)
summary(res_table_2., alpha=0.05)

# save tables as csvs
write.table(res_table_1., file="table_1..csv", sep="\t", quote=F, col.names=NA)
write.table(res_table_2., file="table_2..csv", sep="\t", quote=F, col.names=NA)

# define thresholds, TRUE genes meet criteria, FALSE not
padj.cutoff <- 0.05
lfc.cutoff <- 0.58
threshold_1 <- res_table_1$padj < padj.cutoff & abs(res_table_1$log2FoldChange) > lfc.cutoff
threshold_2 <- res_table_2$padj < padj.cutoff & abs(res_table_2$log2FoldChange) > lfc.cutoff

# of genes deferentially expressed, 546-1946-555, according to our criteria (comparisons)

# add a column to results table "res_table"
res_table_1$threshold <- threshold_1
res_table_2$threshold <- threshold_2

# only include those that are significant, CANNOT SAVE THESE TABLES
subset(res_table_1, threshold == TRUE)
subset(res_table_2, threshold == TRUE)

# visualize specific genes-1st way
plotCounts(dds, gene="Gli1", intgroup="sampletype")
plotCounts(dds, gene="Gli2", intgroup="sampletype")
plotCounts(dds, gene="Mycn", intgroup="sampletype")
plotCounts(dds, gene="Ptch1", intgroup="sampletype")
plotCounts(dds, gene="Hhip", intgroup="sampletype")
plotCounts(dds, gene="Smo", intgroup="sampletype")
# visualize specific genes-2nd way
Hhip<- plotCounts(dds, gene="Hhip", intgroup="sampletype", returnData=TRUE)
ggplot(Hhip, aes(x = sampletype, y = count, color = sampletype)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(Hhip))) + 
  theme_bw() +
  ggtitle("Hhip") +
  theme(plot.title = element_text(hjust = 0.5))

# volcano plot, create data frame for plotting
df_1 <- data.frame(res_table_1)
df_2 <- data.frame(res_table_2)

# plot it-1st way
ggplot(df_1) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold_1)) +
  xlim(c(-2,2))+
  ggtitle('Differential expression') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))

  ggplot(df_2) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold_2)) +
  xlim(c(-2,2))+
  ggtitle('Differential expression') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggplot(df_3) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour=threshold_3)) +
  xlim(c(-2,2)) +
  ggtitle('Differential expression') +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                           panel.background = element_blank(), axis.line = element_line(colour = "black"))


# plot it-2nd way
ggplot(data=df_1, aes(x=log2FoldChange, y=-log10(padj))) + geom_point() + theme_gray()

hypomethyl <-  read.csv("DEG.AZA48_Hypomethyl.AZA24.csv",  header = T, sep = ",")
ggplot(df_2) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), colour="Column1")) +
  xlim(c(-2,2))+ 
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)),
        axis.text = element_text(size = rel(1.25)))


# heatmap, sort the results tables, order the rows of res-tables according to column padj
res_table_1_sorted <- res_table_1[order(res_table_1$padj), ]
# get significant genes, creates a vector (row.names), for which threshold is TRUE! DOESNT WORK
sig1 <- row.names(res_table_1_sorted)[which(res_table_1_sorted$threshold_1)]
# Extract normalized counts for significant genes
norm_1sig <- normalized_counts[sig1,]



