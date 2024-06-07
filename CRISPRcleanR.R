rm(list = ls())
library(openxlsx)
library(ggplot2)
library(pheatmap)
library(org.Mm.eg.db)
library(dplyr)
library(CRISPRcleanR)

setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/input")

reads_SMB21 = read.delim('counts_SMB21.txt')
reads_DAOY = read.delim("counts_DAOY.txt")
reads_DAOY = reads_DAOY[,-6]
colnames(reads_SMB21)[3] <- "SMB21_A"
colnames(reads_SMB21)[4] <- "SMB21_B"
colnames(reads_SMB21)[5] <- "SMB21_C"
colnames(reads_SMB21)[6] <- "Brie_reference"
colnames(reads_DAOY)[3] <- "DAOY_A"
colnames(reads_DAOY)[4] <- "DAOY_B"
colnames(reads_DAOY)[5] <- "DAOY_C"
colnames(reads_DAOY)[6] <- "Brunello_reference"
Brunello_reference = read.csv(file = "Brunello_reference.csv", header = T)
reads_DAOY = reads_DAOY[,-6]
reads_DAOY$Brunello_reference = Brunello_reference$Brunello_reference

#include gene information for Brie and Brunello, with correct code for crisprcleanR, which also leads to deletion of non-targeting controls
#Brunello first
code = read.delim('gRNA_to_Brunello_code.txt')
raw_reads_DAOY = cbind(code, reads_DAOY[,3:6])

#Brie, make library annotation first as it's not included in crisprcleanR, information taken from addgene and chr number generated below
x <- org.Mm.egCHR
mapped_genes <- mappedkeys(x)
xx = as.data.frame(as.numeric(x[mapped_genes]))
write.table(xx, '\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/input/Brie_chr_information.txt')

Brie.anno = read.csv(file="Brie_anno.csv")
rownames(Brie.anno) = Brie.anno[,1]
Brie.anno = Brie.anno[order(rownames(Brie.anno)),]
Brie.anno = Brie.anno[,-1]
str(Brie.anno)

Brie.anno.clean = read.csv(file="Brie_anno_clean.csv")
rownames(Brie.anno.clean) = Brie.anno.clean[,1]
Brie.anno.clean = Brie.anno.clean[order(rownames(Brie.anno.clean)),]
Brie.anno.clean = Brie.anno.clean[,-1]

raw_reads_SMB21 = cbind(Brie.anno, reads_SMB21[,3:6])
str(raw_reads_SMB21)

###perform crisprcleanR
#do normalization, LFC, and correction of LFCs and read counts individually for each replicate  
#make dataframe for individual replicates and plasmid first
#DAOY
#REP A
data("Brunello_Library")
counts_DAOY_A = cbind.data.frame(raw_reads_DAOY[,"Code"], raw_reads_DAOY[,"Gene"],raw_reads_DAOY[,"Brunello_reference"], raw_reads_DAOY[,"DAOY_A"])
colnames(counts_DAOY_A) = c("sgRNA", "gene", "plasmid", "DAOY_A")
str(counts_DAOY_A)
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_DAOY_A,
                                saveToFig = TRUE,
                                libraryAnnotation=Brunello_Library,
                                EXPname='DAOY_A',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="DAOY_A_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="DAOY_A_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brunello_Library)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "DAOY_A",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")

write.csv(correctedFCs$corrected_logFCs, file="DAOY_A_corrected_LFCs.csv")

#correct read counts
correctedCounts<-ccr.correctCounts("DAOY_A",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brunello_Library)
write.csv(correctedCounts, file="DAOY_A_corrected_counts.csv")

#REP B
data("Brunello_Library")
counts_DAOY_B = cbind.data.frame(raw_reads_DAOY[,"Code"], raw_reads_DAOY[,"Gene"],raw_reads_DAOY[,"Brunello_reference"], raw_reads_DAOY[,"DAOY_B"])
colnames(counts_DAOY_B) = c("sgRNA", "gene", "plasmid", "DAOY_B")
str(counts_DAOY_B)
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_DAOY_B,
                                saveToFig = TRUE,
                                libraryAnnotation=Brunello_Library,
                                EXPname='DAOY_B',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="DAOY_B_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="DAOY_B_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brunello_Library)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "DAOY_B",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")

write.csv(correctedFCs$corrected_logFCs, file="DAOY_B_corrected_LFCs.csv")

#correct read counts
correctedCounts<-ccr.correctCounts("DAOY_B",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brunello_Library)
write.csv(correctedCounts, file="DAOY_B_corrected_counts.csv")

#REP C
data("Brunello_Library")
counts_DAOY_C = cbind.data.frame(raw_reads_DAOY[,"Code"], raw_reads_DAOY[,"Gene"],raw_reads_DAOY[,"Brunello_reference"], raw_reads_DAOY[,"DAOY_C"])
colnames(counts_DAOY_C) = c("sgRNA", "gene", "plasmid", "DAOY_C")
str(counts_DAOY_C)
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_DAOY_C,
                                saveToFig = TRUE,
                                libraryAnnotation=Brunello_Library,
                                EXPname='DAOY_C',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="DAOY_C_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="DAOY_C_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brunello_Library)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "DAOY_C",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")

write.csv(correctedFCs$corrected_logFCs, file="DAOY_C_corrected_LFCs.csv")

#correct read counts
correctedCounts<-ccr.correctCounts("DAOY_C",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brunello_Library)
write.csv(correctedCounts, file="DAOY_C_corrected_counts.csv")

#all replicates
data("Brunello_Library")
counts_DAOY_all = cbind.data.frame(raw_reads_DAOY[,"Code"], raw_reads_DAOY[,"Gene"],raw_reads_DAOY[,"Brunello_reference"], raw_reads_DAOY[,"DAOY_A"],raw_reads_DAOY[,"DAOY_B"],raw_reads_DAOY[,"DAOY_C"] )
colnames(counts_DAOY_all) = c("sgRNA", "gene", "plasmid","DAOY_A","DAOY_B", "DAOY_C")
str(counts_DAOY_all)
getwd()
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/02 CRISPRcleanR/Output/crisprcleanR")
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_DAOY_all,
                                saveToFig = TRUE,
                                libraryAnnotation=Brunello_Library,
                                EXPname='DAOY_all',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/02 CRISPRcleanR/Output/crisprcleanR")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/02 CRISPRcleanR/Output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="DAOY_all_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="DAOY_all_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brunello_Library)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "DAOY_all",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/02 CRISPRcleanR/Output/crisprcleanR")

write.csv(correctedFCs$corrected_logFCs, file="DAOY_all_corrected_LFCs.csv")

#get gene level LFCs for DAOY
FCs_DAOY = correctedFCs$corrected_logFCs$avgFC
names(FCs_DAOY) = rownames(correctedFCs$corrected_logFCs)
geneFCs_DAOY = ccr.geneMeanFCs(FCs_DAOY, Brunello_Library)
getwd()
write.csv(geneFCs_DAOY, file="LFCs_DAOY_gene_level.csv")


#correct read counts
correctedCounts<-ccr.correctCounts("DAOY_all",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brunello_Library)
write.csv(correctedCounts, file="DAOY_all_corrected_counts.csv")



#SMB21
#REP A
counts_SMB21_A = cbind.data.frame(rownames(raw_reads_SMB21), raw_reads_SMB21[,"GENES"],raw_reads_SMB21[,"Brie_reference"], raw_reads_SMB21[,"SMB21_A"])
colnames(counts_SMB21_A) = c("sgRNA", "gene", "plasmid", "SMB21_A")
str(counts_SMB21_A)
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_SMB21_A,
                                saveToFig = TRUE,
                                libraryAnnotation=Brie.anno,min_reads = 0,
                                EXPname='SMB21_A',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="SMB21_A_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="SMB21_A_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brie.anno)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "SMB21_A",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")

write.csv(correctedFCs$corrected_logFCs, file="SMB21_A_corrected_LFCs.csv")

#correct read counts
correctedCounts<-ccr.correctCounts("SMB21_A",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brie.anno)
write.csv(correctedCounts, file="SMB21_A_corrected_counts.csv")

#REP B
counts_SMB21_B = cbind.data.frame(rownames(raw_reads_SMB21), raw_reads_SMB21[,"GENES"],raw_reads_SMB21[,"Brie_reference"], raw_reads_SMB21[,"SMB21_B"])
colnames(counts_SMB21_B) = c("sgRNA", "gene", "plasmid", "SMB21_B")
str(counts_SMB21_B)
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_SMB21_B,
                                saveToFig = TRUE,
                                libraryAnnotation=Brie.anno,min_reads = 0,
                                EXPname='SMB21_B',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="SMB21_B_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="SMB21_B_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brie.anno)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "SMB21_B",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")

write.csv(correctedFCs$corrected_logFCs, file="SMB21_B_corrected_LFCs.csv")

#correct read counts
correctedCounts<-ccr.correctCounts("SMB21_B",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brie.anno)
write.csv(correctedCounts, file="SMB21_B_corrected_counts.csv")

#REP C
counts_SMB21_C = cbind.data.frame(rownames(raw_reads_SMB21), raw_reads_SMB21[,"GENES"],raw_reads_SMB21[,"Brie_reference"], raw_reads_SMB21[,"SMB21_C"])
colnames(counts_SMB21_C) = c("sgRNA", "gene", "plasmid", "SMB21_C")
str(counts_SMB21_C)
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_SMB21_C,
                                saveToFig = TRUE,
                                libraryAnnotation=Brie.anno,min_reads = 0,
                                EXPname='SMB21_C',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="SMB21_C_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="SMB21_C_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brie.anno)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "SMB21_C",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/Dependency_screens/output/crisprcleanR/")

write.csv(correctedFCs$corrected_logFCs, file="SMB21_C_corrected_LFCs.csv")

#correct read counts
correctedCounts<-ccr.correctCounts("SMB21_C",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brie.anno)
write.csv(correctedCounts, file="SMB21_C_corrected_counts.csv")

#all replicates
counts_SMB21_all = cbind.data.frame(rownames(raw_reads_SMB21), raw_reads_SMB21[,"GENES"],raw_reads_SMB21[,"Brie_reference"], raw_reads_SMB21[,"SMB21_A"],raw_reads_SMB21[,"SMB21_B"],raw_reads_SMB21[,"SMB21_C"])
colnames(counts_SMB21_all) = c("sgRNA", "gene", "plasmid","SMB21_A","SMB21_B", "SMB21_C")
str(counts_SMB21_all)
normANDfcs<-ccr.NormfoldChanges(Dframe=counts_SMB21_all,
                                saveToFig = TRUE,
                                libraryAnnotation=Brie.anno,min_reads = 0,
                                EXPname='SMB21_all',
                                outdir = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/02 CRISPRcleanR/Output/crisprcleanR")
setwd("\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/02 CRISPRcleanR/Output/crisprcleanR")
write.csv(normANDfcs$norm_counts, file="SMB21_all_norm_counts.csv")
write.csv(normANDfcs$logFCs, file="SMB21_all_logFCs.csv")

#correct LFCs
gwSortedFCs <- ccr.logFCs2chromPos(normANDfcs$logFCs,Brie.anno)
correctedFCs <- ccr.GWclean(gwSortedFCs, display = TRUE, label = "SMB21_all",min.ngenes=5,
                            saveTO = "\\\\fileserv1/groups/Neuro-Oncology/Clinical Neuro-Oncology/Daniel/SHH MB Screens/02 CRISPRcleanR/Output/crisprcleanR")

write.csv(correctedFCs$corrected_logFCs, file="SMB21_all_corrected_LFCs.csv")

#get gene level LFCs for SMB21
FCs_SMB21 = correctedFCs$corrected_logFCs$avgFC
names(FCs_SMB21) = rownames(correctedFCs$corrected_logFCs)
geneFCs_SMB21 = ccr.geneMeanFCs(FCs_SMB21, Brie.anno)
getwd()
write.csv(geneFCs_SMB21, file="LFCs_SMB21_gene_level.csv")



#correct read counts
correctedCounts<-ccr.correctCounts("SMB21_all",
                                   normANDfcs$norm_counts,
                                   correctedFCs,
                                   Brie.anno)
write.csv(correctedCounts, file="SMB21_all_corrected_counts.csv")



