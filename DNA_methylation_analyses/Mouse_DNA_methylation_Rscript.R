library(sesame)
library(readxl)
library(umap)
library(ComplexHeatmap)
library(ggplot2)
library(tidyverse)
library(SummarizedExperiment)
library(sesameData)
sesameDataCache()
library(pals)
install.packages("gprofiler2")
library(gprofiler2)
install.packages("data.table")
library(data.table)

getwd()
# get experiment information
experiment_anno = read_xlsx("experiment_anno.xlsx")

# make an anno matrix
anno_methylation = as.matrix(experiment_anno)
rownames(anno_methylation) = experiment_anno$Sample_Name

# read and preprocess idat files, open sesame also performs QC (e.g. detection p values)
# total number of probes equals 297,415 on CHIP; 268,966 after QC  
getwd()
idat_dir = "/Users/lab/Documents/Foteini/Mouse_methylation/data"

#get all betas and detP for GEO submission
pvalspOOBAH = openSesame(idat_dir, func = pOOBAH, return.pval = TRUE)
write.csv(pvalspOOBAH, file = "Mouse_methylation_all_detP.csv")
betas_all = openSesame(idat_dir, mask=FALSE)
betas_all_nomask = openSesame(idat_dir)
write.csv(betas_all, file = "Mouse_methylation_all_Betas.csv")

betas. = openSesame(idat_dir)
dim(betas.)
betas.NA.omit = na.omit(betas.)
dim(betas.NA.omit)
write.table(betas., file= "betavalues.csv", sep=",", row.names=FALSE)

# rename and reorder columns
#but the command for betas.df is missing!!!!
betas.2 <- betas.
colnames(betas.2) <- c("SSMB55_AZA_24h_A", "SMB55_AZA_24h_B", "SMB55_AZA_24h_C", "SMB55_DMSO_A", "SMB55_DMSO_B",  "SMB55_DMSO_C",  "SMB21_AZA_24h_A", "SMB21_AZA_24h_B",  "SMB21_AZA_24h_C", "SMB21_DMSO_A",  "SMB21_DMSO_B",  "SMB21_DMSO_C")
colnames(betas.2)
betas.2. <- betas.2[, c(1, 3, 5, 7, 9, 11, 2, 4, 6, 8, 10, 12, 13, 15, 17, 19, 21, 23, 14, 16, 18, 20, 22, 24)]

betas.MB55 <- betas.2. [,-(13:24)]
betas.MB21 <- betas.2. [,-(1:12)]

# differential methylation
#first, make summarizedExperiment
#for each cell line individually
## Start with SMB21 first
#needs to be dataframe, beta already dataframe, make annotation also dataframe
sampledata = as.data.frame(anno_methylation_21)
#make summarizedExperiment
library(SummarizedExperiment)
se.MB21 <- SummarizedExperiment(as.matrix(betas.MB21), colData = sampledata)
colData(se.MB21)
assay(se.MB21)
dim(se.MB21)
#get rid of NA probes (before 297,415, after 275,380)
se.ok = (checkLevels(assay(se.MB21), colData(se.MB21)$Condition))
sum(se.ok)
se.MB21 = se.MB21[se.ok,]
dim(se.MB21)
#turn discrete contrast variables to factor (before chr, after factor) with reference level equals DMSO
str(se.MB21$Condition)
colData(se.MB21)$Condition <- relevel(factor(colData(se.MB21)$Condition), "SMB21_DMSO")
str(se.MB21$Condition)
#model DNA methylation variation
smry.MB21 = DML(se.MB21, ~Condition, mc.cores = 6)
lm_results_MB21 = summaryExtractTest(smry.MB21)
head(lm_results_MB21)
str(lm_results_MB21)
lm_results_MB21$Pval_ConditionSMB21_AZA_24h
getwd()
setwd("/Users/lab/Documents/Foteini/Mouse_methylation/results")
write.table(lm_results_MB21, file ="lm_results_MB21.csv")
# adjust p values
padj.AZA24 = p.adjust(lm_results_MB21$Pval_ConditionSMB21_AZA_24h, method = "BH", n = length(lm_results_MB21$Pval_ConditionSMB21_AZA_24h))
lm_results_MB21_padj = cbind(lm_results_MB21, padj.AZA24)
padj.AZA2 = p.adjust(lm_results_MB21$Pval_ConditionSMB21_AZA_2h, method = "BH", n = length(lm_results_MB21$Pval_ConditionSMB21_AZA_2h))
lm_results_MB21_padj_ = cbind(lm_results_MB21_padj, padj.AZA2)
padj.SGI24 = p.adjust(lm_results_MB21$Pval_ConditionSMB21_SGI_24h, method = "BH", n = length(lm_results_MB21$Pval_ConditionSMB21_SGI_24h))
lm_results_MB21_padj. = cbind(lm_results_MB21_padj_,padj.SGI24)
write.table(lm_results_MB21_padj., file ="lm_results_MB21_padj.csv")
whatever = cbind(lm_results_MB21_padj_,padj.SGI24)

#pairwise comparison
library(ggplot2)
# p value-based figures
ggplot(lm_results_MB21) + geom_point(aes(Est_ConditionSMB21_AZA_24h, -log10(Pval_ConditionSMB21_AZA_24h)))

# p adjusted-based figures
# color=x & size=y
ggplot(lm_results_MB21_padj.) + geom_point(aes(x=Est_ConditionSMB21_AZA_24h, y=-log10(padj.AZA24), color=Est_ConditionSMB21_AZA_24h, size=-log10(padj.AZA24)))+scale_color_gradient2(low = "steelblue", mid= "steelblue",high = "darkred", midpoint = 0,space = "Lab") + 
  scale_size(range = c(0.3,3))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# color=y & size=x
getwd()
setwd("/Users/lab/Documents/Foteini/Mouse_methylation/results")

ggplot(lm_results_MB21_padj.) + geom_point(aes(x=Est_ConditionSMB21_AZA_24h, y=-log10(padj.AZA24), color=-log10(padj.AZA24), size=abs(Est_ConditionSMB21_AZA_24h)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.1, 0)) + 
  scale_size(range = c(0.3,3))+
  scale_y_continuous(breaks = seq(0,11,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

#DMR
dmContrasts(smry.MB21)
DMR.MB21.AZA24 = DMR(se.MB21, smry.MB21, "ConditionSMB21_AZA_24h")
DMR.MB21.AZA24 %>% dplyr::filter(Seg_Pval_adj < 0.01)

#DMR does not seem to a big help, continue with diff CpGs

#get probe IDs from modeling that have at least downregulation of 0.1 and padj < 0.05
# 73594 probes fit the criteria
intermediate = lm_results_MB21_padj. %>% dplyr::filter(Est_ConditionSMB21_AZA_24h < -0.1)
intermediate = intermediate %>% dplyr::filter(padj.AZA24 < 0.05)
sigdiffprobesdown.MB21.AZA24 = as.vector(intermediate$Est_Probe_ID)
sigdiffprobesdown.MB21.AZA24

#get genes associated with hypomethylated CpGs in AZA 24hrs, (n = 630) 
genes.down.MB21.AZA24 = testEnrichment(sigdiffprobesdown.MB21.AZA24, KYCG_buildGeneDBs(sigdiffprobesdown.MB21.AZA24))
genes.down.MB21.AZA24.sig = genes.down.MB21.AZA24 %>% dplyr::filter(FDR < 0.1)
genes.down.MB21.AZA24.sig = genes.down.MB21.AZA24.sig$gene_name
genes.down.MB21.AZA24.sig
write.csv(genes.down.MB21.AZA24.sig, file = "Genes_down_MB21_AZA24_sig.csv")

# get probe IDs with upregulation of min +0.1 and padj<0.05
# 0 genes for 24h Aza
intermediate.hyper.AZA.24 = lm_results_MB21_padj. %>% dplyr::filter(Est_ConditionSMB21_AZA_24h > 0.1)
intermediate.hyper.AZA.24 = intermediate.hyper.AZA.24 %>% dplyr::filter(padj.AZA24 < 0.05)

#GO analysis for GENES associated with hypomethylation; AZA 24hrs
getwd()
go.genes.MB21.AZA24.down <- gost(genes.down.MB21.AZA24.sig, organism = "mmusculus", sources = "GO")
go.genes.MB21.AZA24.down$result[order(go.genes.MB21.AZA24.down$result$p_value),]
go.genes.MB21.AZA24.down$result
results.GO.MB21.AZA24.down = go.genes.MB21.AZA24.down$result
fwrite(results.GO.MB21.AZA24.down, file = "GO_MB21_AZA24_down.csv")
gostplot(go.genes.MB21.AZA24.down, capped = FALSE, interactive = FALSE)


# differential methylation
## SMB55 ! 
anno_methylation_55 <- anno_methylation[-(13:24), ]
sampledata.55 = as.data.frame(anno_methylation_55)
#make summarizedExperiment
library(SummarizedExperiment)
se.MB55 <- SummarizedExperiment(as.matrix(betas.MB55), colData = sampledata.55)
colData(se.MB55)
assay(se.MB55)
dim(se.MB55)
#get rid of NA probes (before 297,415, after 272,886)
se.ok.55 = (checkLevels(assay(se.MB55), colData(se.MB55)$Condition))
sum(se.ok.55)
se.MB55 = se.MB55[se.ok.55,]
#turn discrete contrast variables to factor (before chr, after factor) with reference level equals DMSO
str(se.MB55$Condition)
colData(se.MB55)$Condition <- relevel(factor(colData(se.MB55)$Condition), "SMB55_DMSO")
str(se.MB55$Condition)
#model DNA methylation variation
smry.MB55 = DML(se.MB55, ~Condition, mc.cores = 6)
# cores=6 gives me an error
smry.MB55 = DML(se.MB55, ~Condition)
lm_results_MB55 = summaryExtractTest(smry.MB55)
setwd("/Users/lab/Documents/Foteini/Mouse_methylation/results")
write.table(lm_results_MB55, file ="lm_results_MB55.csv")
# adjust p values
padj.MB55.AZA24 = p.adjust(lm_results_MB55$Pval_ConditionSMB55_AZA_24h, method = "BH", n = length(lm_results_MB55$Pval_ConditionSMB55_AZA_24h))
lm_results_MB55_padj = cbind(lm_results_MB55, padj.MB55.AZA24)
padj.MB55.AZA2 = p.adjust(lm_results_MB55$Pval_ConditionSMB55_AZA_2h, method = "BH", n = length(lm_results_MB55$Pval_ConditionSMB55_AZA_2h))
lm_results_MB55_padj_ = cbind(lm_results_MB55_padj, padj.MB55.AZA2)
padj.MB55.SGI24 = p.adjust(lm_results_MB55$Pval_ConditionSMB55_SGI_24h, method = "BH", n = length(lm_results_MB55$Pval_ConditionSMB55_SGI_24h))
lm_results_MB55_padj. = cbind(lm_results_MB55_padj_,padj.MB55.SGI24)
write.table(lm_results_MB55_padj., file ="lm_results_MB55_padj.csv")

#pairwise comparison
# p value-based figures
ggplot(lm_results_MB55) + geom_point(aes(Est_ConditionSMB55_AZA_24h, -log10(Pval_ConditionSMB55_AZA_24h)))
ggplot(lm_results_MB55) + geom_point(aes(Est_ConditionSMB55_AZA_2h, -log10(Pval_ConditionSMB55_AZA_2h)))
ggplot(lm_results_MB55) + geom_point(aes(Est_ConditionSMB55_SGI_24h, -log10(Pval_ConditionSMB55_SGI_24h)))

# p adjusted-based figures
# color=y & size=x
ggplot(lm_results_MB55_padj.) + geom_point(aes(x=Est_ConditionSMB55_AZA_24h, y=-log10(padj.MB55.AZA24), color=-log10(padj.MB55.AZA24), size=abs(Est_ConditionSMB55_AZA_24h)))+
  scale_color_gradientn(colours = c("red", "steelblue", "darkblue"), values = c(1,0.1, 0)) + 
  scale_size(range = c(0.3,3))+
  scale_y_continuous(breaks = seq(0,11,2))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))

# identify hxpomethyl probes of 0.1 and padj < 0.05 
# 98794 probes fit the criteria
intermediate.MB55 = lm_results_MB55_padj. %>% dplyr::filter(Est_ConditionSMB55_AZA_24h < -0.1)
intermediate.MB55. = intermediate.MB55 %>% dplyr::filter(padj.MB55.AZA24 < 0.05)
sigdiffprobesdown.MB55.AZA24 = as.vector(intermediate.MB55.$Probe_ID)
sigdiffprobesdown.MB55.AZA24

#get genes associated with hypomethylated CpGs in AZA 24hrs (n = 607)
genes.down.MB55.AZA24 = testEnrichment(sigdiffprobesdown.MB55.AZA24, KYCG_buildGeneDBs(sigdiffprobesdown.MB55.AZA24))
genes.down.MB55.AZA24.sig = genes.down.MB55.AZA24 %>% dplyr::filter(FDR < 0.1)
genes.down.MB55.AZA24.sig = genes.down.MB55.AZA24.sig$gene_name
genes.down.MB55.AZA24.sig
write.csv(genes.down.MB55.AZA24.sig, file = "Genes_down_MB55_AZA24_sig.csv")

# identify hypermethyl probes of +0.1 and padj<0.05
# 1 probe for 24h Aza
intermediate.hyper.MB55.AZA.24 = lm_results_MB55_padj. %>% dplyr::filter(Est_ConditionSMB55_AZA_24h > 0.1)
intermediate.hyper.MB55.AZA.24. = intermediate.hyper.MB55.AZA.24 %>% dplyr::filter(padj.MB55.AZA24 < 0.05)
sigdiffprobesup.MB55.AZA24 = as.vector(intermediate.hyper.MB55.AZA.24.$Probe_ID)

#get genes associated with hypermethylated CpGs in AZA 24hrs 
genes.up.MB55.AZA24 = testEnrichment(sigdiffprobesup.MB55.AZA24, KYCG_buildGeneDBs(sigdiffprobesup.MB55.AZA24))
genes.up.MB55.AZA24.sig = genes.up.MB55.AZA24 %>% dplyr::filter(FDR < 0.1)
genes.up.MB55.AZA24.sig  = genes.up.MB55.AZA24.sig $gene_name
write.csv(genes.up.MB55.AZA24.sig, file = "Genes_up_MB55_AZA24_sig.csv")



