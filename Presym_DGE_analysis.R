###############################################################################
# Author: Gabriel Ferreras Garrucho                                           #
# Script template: Tximport + DESeq2 RNA-Seq analysis                         #
# Experiment:     Pre-symbiotic Network                                       #
###############################################################################

# LOADING PACKAGES ############################################################

library(tximport)
library(DESeq2)
library(rhdf5)
library(genefilter)
library(plotly)
library(RNAseqQC)
library(apeglm)
library(BiocParallel)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(scales)
library(ggtext)
library(ggh4x)
library(ggpmisc)

dir.create("plots/")
dir.create("plots/QC/")
dir.create("plots/volcano_plots/")
dir.create("gene_count_tables/")
dir.create("DEG_tables/")
dir.create("comparison_tables/")


register(MulticoreParam(4))

treatment_colours<-c("gray30",
                     "palegreen3","springgreen3", "springgreen4", "olivedrab3", "olivedrab4", 
                     "goldenrod1","darkorchid1", "deeppink")

genotype_colours<-c("gray30","darkorchid1","deepskyblue","dodgerblue","goldenrod1",
                    "palegreen2", "springgreen3", "chartreuse2")


treatment_shapes<-c(19,0,7,12,14,22,15,18,17)
genotype_shapes<-c(19,25,7,8,1,13,11,14)

sampleTable<-readxl::read_xlsx("info_sheets/sample_info.xlsx", sheet=1)
sampleTable$condition<-factor(sampleTable$condition, levels=unique(sampleTable$condition))
sampleTable$genotype<-factor(sampleTable$genotype, levels=unique(sampleTable$genotype))
sampleTable$genotype_name<-factor(sampleTable$genotype_name, levels=unique(sampleTable$genotype_name))
sampleTable$treatment_name<-factor(sampleTable$treatment_name, levels=unique(sampleTable$treatment_name))
sample_list<-sampleTable$name

names(treatment_colours)<-unique(sampleTable$treatment_name)
names(genotype_colours)<-unique(sampleTable$genotype_name)

publication_theme <- function() 
{
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(face="bold", color="black", size=10),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size=10, color="black", face="italic"),
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=10, color="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(angle=0, vjust=0.5, size=10, color="black"),
    axis.title.y = element_text(face="bold", color="black", size=10),
    axis.title.x = element_text(face="bold", color="black", size=10),
    plot.title = element_text (face="bold", color="black", hjust = 0.5),
    strip.text = element_text(size=10, color="white", face="italic"),
    strip.background = element_rect(color="black", fill="black", linewidth=1.5, linetype="solid")
  )
}

# TXIMPORT ##############################################################
# Importing data from Kallisto counting results
# Not required to run this script as the gene count matrix is provided.

files_list<-paste("../samples","/",sample_list,"/","abundance.h5",sep="")
names(files_list)<-sample_list
tx2gene.table_RAP<-read_csv("info_sheets/tx2gene_table_RAP.csv")
nrow(tx2gene.table_RAP)
txi.kallisto <- tximport(files_list, type = "kallisto", txOut = F, tx2gene = tx2gene.table_RAP)
head(txi.kallisto$counts)
nrow(txi.kallisto$counts)
write.csv(txi.kallisto$counts, file="gene_count_tables/gene_count_matrix.csv")
write.csv(txi.kallisto$length, file="gene_count_tables/gene_length_matrix.csv")
write.csv(txi.kallisto$abundance, file="gene_count_tables/gene_abundance_matrix.csv")

counts <- gene_count_matrix
counts_df <- as.data.frame(counts)
counts_long <- counts_df %>% pivot_longer(cols = everything(), names_to = "name", values_to = "count")
counts_long <- merge(counts_long, sampleTable, by = "name")
ggplot(counts_long, aes(x = name, y = count, fill = genotype_name)) +
  geom_violin(alpha = 0.6, linewidth=0.25, color="black") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA, linewidth=0.2) +
  labs(x = "Sample", y = "Counts", title = "Kallisto Count Distribution") +
  publication_theme() +
  scale_fill_manual(name = "Genotype", values=genotype_colours) +
  scale_y_log10(labels = scales::comma_format(scale = 10))
ggsave(filename = "plots/QC/kallisto_counts_boxplot.pdf", width = 15, height = 4)

# DESEQ2 #####################################################################

# Reading raw gene count matrix, naming the columns and creating DESeq DataSet.

gene_count_matrix <- read.csv("gene_count_tables/gene_count_matrix.csv", row.names = 1)
gene_count_matrix <- round(as.matrix(gene_count_matrix))

#dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~ condition)
dds <- DESeqDataSetFromMatrix(gene_count_matrix, sampleTable, ~ condition)
colData(dds)
head(counts(dds))
nrow(counts(dds))

# Quality control

plot_total_counts(dds)
ggsave(filename = "plots/QC/total_counts.pdf", width = 5, height = 5)
plot_library_complexity(dds)
ggsave(filename = "plots/QC/library_complexity.pdf", width = 5, height = 5)
plot_gene_detection(dds)
ggsave(filename = "plots/QC/gene_detection.pdf", width = 5, height = 5)

counts <- counts(dds)
counts_df <- as.data.frame(counts)
counts_long <- counts_df %>% pivot_longer(cols = everything(), names_to = "name", values_to = "count")
counts_long <- merge(counts_long, sampleTable, by = "name")
ggplot(counts_long, aes(x = name, y = count, fill = genotype_name)) +
  geom_violin(alpha = 0.6, linewidth=0.25, color="black") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA, linewidth=0.2) +
  labs(x = "Sample", y = "Counts", title = "DESeq2 Raw Count Distribution") +
  publication_theme() +
  scale_fill_manual(name = "Genotype", values=genotype_colours) +
  scale_y_log10(labels = scales::comma_format(scale = 10))
ggsave(filename = "plots/QC/DESEq2_counts_boxplot.pdf", width = 15, height = 4)

# Prefiltering: eliminating all genes with less than 10 reads between all samples

keep <- rowSums(counts(dds)) >= 10
dds<-dds[keep,]
head(counts(dds))
nrow(counts(dds))

# Factoring the DESeq DataSet
dds$condition <- factor(dds$condition, levels = levels(sampleTable$condition))
dds$condition <- relevel(dds$condition, ref = "A01")

# Obtaining DESeq results
dds <- DESeq(dds, parallel=TRUE)

normalized_counts <- counts(dds, normalized=TRUE)
write.csv(counts(dds, normalized=F), file="gene_count_tables/DESeq2_counts.csv")
write.csv(normalized_counts, file="gene_count_tables/DESeq2_normalized_counts.csv")
counts_df <- as.data.frame(normalized_counts)
counts_long <- counts_df %>% pivot_longer(cols = everything(), names_to = "name", values_to = "count")
counts_long <- merge(counts_long, sampleTable, by = "name")
ggplot(counts_long, aes(x = name, y = count, fill = genotype_name)) +
  geom_violin(alpha = 0.6, linewidth=0.25, color="black") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA, linewidth=0.2) +
  labs(x = "Sample", y = "Counts", title = "DESeq2 Normalised Count Distribution") +
  publication_theme() +
  scale_fill_manual(name = "Genotype", values=genotype_colours) +
  scale_y_log10(labels = scales::comma_format(scale = 10))
ggsave(filename = "plots/QC/DESEq2_normalised_counts_boxplot.pdf", width = 15, height = 4)

# Variance-stabilising transformation
vsd <- vst(dds, blind=FALSE)
mean_sd_plot(vsd)
ggsave(filename = "plots/QC/mean_sd_plot.pdf", width = 5, height = 5)
vsd_counts <- assay(vsd)
write.csv(vsd_counts, file="gene_count_tables/DESeq2_VSD_counts.csv")

counts <- assay(vsd)
counts_df <- as.data.frame(counts)
counts_long <- counts_df %>% pivot_longer(cols = everything(), names_to = "name", values_to = "count")
counts_long <- merge(counts_long, sampleTable, by = "name")
ggplot(counts_long, aes(x = name, y = count, fill = genotype_name)) +
  geom_violin(alpha = 0.6, linewidth=0.25, color="black") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.shape = NA, linewidth=0.2) +
  labs(x = "Sample", y = "Counts", title = "DESeq2 VST Count Distribution") +
  publication_theme() +
  scale_fill_manual(name = "Genotype", values=genotype_colours) +
  scale_y_log10(labels = scales::comma_format(scale = 10))
ggsave(filename = "plots/QC/VSD_normalised_counts_boxplot.pdf", width = 15, height = 4)

# Other QC

set.seed(1)
pdf(file = "plots/QC/correlation_map.pdf", width = 30, height = 30)
plot_sample_clustering(vsd, anno_vars = c("treatment_name", "genotype_name","replicate"), distance = "euclidean")
dev.off()

ma_plots <- plot_sample_MAs(vsd, group = "condition")
cowplot::plot_grid(plotlist = ma_plots[1:12], ncol = 3)
ggsave(filename = "plots/QC/replicate_variability.pdf", width = 20, height = 30)

# DEG RESULTS ##################################################################
rice_annotation<-read_tsv("info_sheets/Osativa_gene_annotation_RAPcollapsed.tsv")
tx2gene.table_RAP<-read_csv("info_sheets/tx2gene_table_RAP.csv")

get_DEGs<-function(genotypes, treatments, condition, ref)
{
  name<-paste(condition, "_vs_", ref, sep="")
  
  # Getting data, importing and tximport
  sampleTable_filtered<-filter(sampleTable, genotype %in% genotypes, treatment %in% treatments)
  sample_list_filtered<-sampleTable_filtered$name
  
  files_list_filtered<-paste("samples","/",sample_list_filtered,"/","abundance.h5",sep="")
  names(files_list_filtered)<-sample_list_filtered
  txi.kallisto <- tximport(files_list_filtered, type = "kallisto", txOut = F, tx2gene = tx2gene.table_RAP)
  
  dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable_filtered, ~condition)
  keep <- rowSums(counts(dds)) >= 10
  dds<-dds[keep,]
  
  # Factoring the DESeq DataSet
  dds$condition <- factor(dds$condition, levels = unique(sampleTable_filtered$condition))
  
  # Stablishing reference level
  dds$condition <- relevel(dds$condition, ref = ref)
  
  # Obtaining DESeq results.
  dds <- DESeq(dds, parallel=TRUE)
  
  # Getting the comparison table.
  comparison<- results(dds, name=paste("condition_",name, sep=""), alpha=0.05)
  summary(comparison)
  comparison <- comparison[order(comparison$padj),]
  
  ### MA plots.
  pdf(file = paste("plots/volcano_plots/", name ,"_MA.pdf", sep=""), width = 6, height = 4)
  par(mar = c(5, 5, 3, 3))
  plotMA(comparison, main=name)
  dev.off()
  
  ### Extracting DEGs according to criteria.
  fold.change.comparison <- comparison$log2FoldChange
  p.value.comparison<-comparison$padj
  genes.ids.comparison <- rownames(comparison)
  log.p.value.comparison <- -log10(p.value.comparison)
  names(fold.change.comparison)<- genes.ids.comparison
  names(p.value.comparison)<-genes.ids.comparison
  names(log.p.value.comparison)<- genes.ids.comparison
  
  activated.genes.comparison <- genes.ids.comparison[fold.change.comparison > (log2(1.5)) & p.value.comparison < 0.05]
  activated.genes.comparison<-activated.genes.comparison[!is.na(activated.genes.comparison)]
  repressed.genes.comparison <- genes.ids.comparison[fold.change.comparison < (-log2(1.5)) & p.value.comparison < 0.05]
  repressed.genes.comparison<-repressed.genes.comparison[!is.na(repressed.genes.comparison)]
  
  ### Annotating DEGs.

  annotate_gene_list<-function(gene.list)
  {
    table<-tibble("RAP"=gene.list, 
                  "logFC"=fold.change.comparison[gene.list],
                  "adj.p.value"=p.value.comparison[gene.list])
    table <- table %>% left_join(., rice_annotation, by = "RAP") 
    return(table)
  }

  comparison.table<-as.data.frame(comparison)
  comparison.table<-mutate(comparison.table, RAP=row.names(comparison.table)) %>% dplyr::select(RAP, baseMean, log2FoldChange, lfcSE, pvalue, padj) 
  comparison.table<-left_join(comparison.table, rice_annotation, by = "RAP")
  
  activated.genes.comparison.table<-annotate_gene_list(activated.genes.comparison)
  repressed.genes.comparison.table<-annotate_gene_list(repressed.genes.comparison)
  
  write_tsv(file=paste("DEG_tables/activated_genes_table_",name,".tsv",sep=""), activated.genes.comparison.table)
  write_tsv(file=paste("DEG_tables/repressed_genes_table_", name,".tsv",sep=""), repressed.genes.comparison.table)
  write_tsv(file=paste("comparison_tables/all_genes_table_", name,".tsv",sep=""), comparison.table)
  
  ### Volcano plot
  pdf(file = paste("plots/volcano_plots/", name, "_volcano.pdf", sep=""), width = 8, height = 6)
  par(mar = c(5, 5, 3, 3))
  plot(fold.change.comparison, log.p.value.comparison,
       pch=19,cex=0.3,col="grey",ylab=expression(-log[10](p-value)), xlab=expression(log[2](FC)), 
       main=name)
  
  points(fold.change.comparison[activated.genes.comparison],
         log.p.value.comparison[activated.genes.comparison],
         pch=19,cex=0.3,col="blue")
  
  points(fold.change.comparison[repressed.genes.comparison],
         log.p.value.comparison[repressed.genes.comparison],
         pch=19,cex=0.3,col="red")
  
  abline(v=log2(1), lty = 3, col="grey")
  abline(v=-log2(1), lty = 3, col="grey")
  abline(h=-log10(0.1), lty = 3, col="grey")
  
  legend("topleft", legend=c("Not significant", "Activated genes", "Repressed genes"), 
         col=c("grey", "blue", "red"), pch=19, bty = "n")
  dev.off()
}

get_DEGs(c("A"), c("06", "01"), condition="A06", ref="A01")
get_DEGs(c("A"), c("07", "01"), condition="A07", ref="A01")
get_DEGs(c("A"), c("08", "01"), condition="A08", ref="A01")
get_DEGs(c("A"), c("09", "01"), condition="A09", ref="A01")
get_DEGs(c("A"), c("10", "01"), condition="A10", ref="A01")
get_DEGs(c("A"), c("13", "01"), condition="A13", ref="A01")

get_DEGs(c("A"), c("06", "07"), condition="A07", ref="A06")
get_DEGs(c("A"), c("06", "08"), condition="A08", ref="A06")
get_DEGs(c("A"), c("06", "09"), condition="A09", ref="A06")
get_DEGs(c("A"), c("06", "10"), condition="A10", ref="A06")

get_DEGs(c("A"), c("07", "08"), condition="A08", ref="A07")
get_DEGs(c("A"), c("09", "10"), condition="A10", ref="A09")

for(i in c("01", "11", "12"))
{
  for(j in c("B", "C", "F", "G", "H", "I", "J"))
  {
    get_DEGs(c("A", j), c(i), condition=paste(j,i,sep=""), ref=paste("A",i,sep=""))
  }
}

for(i in c("11", "12"))
{
  for(j in c("A","B", "C", "F", "G", "H", "I", "J"))
  {
    get_DEGs(c(j), c(i, "01"), condition=paste(j,i,sep=""), ref=paste(j,"01",sep=""))
  }
}
