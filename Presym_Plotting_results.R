###############################################################################
# Author: Gabriel Ferreras Garrucho                                           #
# Script template: RNA-Seq analysis plotting results                          #
# Experiment:     Pre-symbiotic Network                                       #
###############################################################################

# LOADING PACKAGES ############################################################

library(plotly)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(multcompView)
library(scales)
library(ggtext)
library(ggh4x)
library(ggpmisc)
library(eulerr)
library(graphics)
library(gplots)
library(ggvenn)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(FSA)
library(ggnewscale)
library(NbClust)
library(cluster)
library(hopkins)
library(factoextra)
library(fpc)
library(clValid)
library(MASS)
library(reshape2)
library(reshape)
library(matrixStats)

dir.create("plots/Average_expression_plots/")
dir.create("plots/Heatmaps/")
dir.create("plots/Venn/")
dir.create("plots/Correlation/")
dir.create("plots/Scatterplots/")
dir.create("plots/Gene_count_plots/")
dir.create("plots/DEG_number/")
dir.create("plots/PCA/")
dir.create("plots/GO/")
dir.create("plots/GO/GO_tables/")
dir.create("plots/GO/GO_plots/")

treatment_colours<-c("gray30",
                     "springgreen3", "springgreen4", "olivedrab3", "olivedrab4", "palegreen3",
                     "goldenrod1", "deeppink","darkorchid1")

genotype_colours<-c("gray30","darkorchid1","deepskyblue","dodgerblue","goldenrod1",
                    "palegreen2", "springgreen3", "chartreuse2")


treatment_shapes<-c(19,0,7,12,14,13,15,18,17)
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

publication_theme_PCA <- function() 
{
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(face="bold", color="black", size=10),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size=10, color="black"),
    axis.text.x = element_text(vjust=0.3, size=10, color="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(angle=0, vjust=0.3, size=10, color="black"),
    axis.title.y = element_text(face="bold", color="black", size=10),
    axis.title.x = element_text(face="bold", color="black", size=10),
    plot.title = element_text (face="bold", color="black", hjust = 0.5)
  )
}

publication_theme_average<- function() 
{
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size=10, color="white"),
    strip.background = element_rect(color="black", fill="black", linewidth=1.5, linetype="solid"),
    legend.position = "none",
    axis.text.x = element_text(angle=90, vjust=0.3, hjust=1, size=10, color="black", face="italic"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(angle=0, vjust=0.3, size=10, color="black"),
    axis.title.y = element_text(face="bold", color="black", size=10),
    axis.title.x = element_text(face="bold", color="black", size=10),
    plot.title = element_text (face="bold", color="black", hjust = 0.5)
  )
}

publication_theme_counts<- function() 
{
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2),
    strip.text = element_text(size=10, color="white"),
    strip.background = element_rect(color="black", fill="black", linewidth=1.5, linetype="solid"),
    legend.position="none",
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(face="bold", color="black", size=10),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size=10, color="black",face="italic"),
    axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=10, color="black", face="italic"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(angle=0, vjust=0.3, size=10, color="black"),
    axis.title.y = element_text(face="bold", color="black", size=10),
    axis.title.x = element_blank(),
    plot.title = element_text (face="bold", color="black", hjust = 0.5)
  )
}

publication_theme_scatterplot <- function() 
{
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(face="bold", color="black", size=10),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size=10, color="black"),
    axis.text.x = element_text(vjust=0.3, size=9, color="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(angle=0, vjust=0.3, size=10, color="black"),
    axis.title.y = element_text(face="bold", color="black", size=10),
    axis.title.x = element_text(face="bold", color="black", size=10),
    plot.title = element_text (face="bold", color="black", hjust = 0.5),
    strip.text = element_text(size=10, color="white", face="italic"),
    strip.background = element_rect(color="black", fill="black", linewidth=1.5, linetype="solid")
  )
}

publication_theme_GO <- function() 
{
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill = "white"),
    legend.title = element_text(face="bold", color="black", size=10),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size=10, color="black"),
    axis.text.x = element_text(vjust=0.3, size=9, color="black"),
    axis.ticks = element_line(color = "black"),
    axis.text.y = element_text(angle=0, vjust=0.3, size=10, color="black"),
    axis.title.y = element_text(face="bold", color="black", size=10),
    axis.title.x = element_text(face="bold", color="black", size=10),
    plot.title = element_text (face="bold", color="black", hjust = 0.5),
    strip.text = element_text(size=10, color="white"),
    strip.background = element_rect(color="black", fill="black", linewidth=1.5, linetype="solid")
  )
}

# LOADING DATA #############################################################
# Compiled rice gene annotation
rice_annotation<-read_tsv("info_sheets/Osativa_gene_annotation_RAPcollapsed.tsv")

# Gene lists from the literature
AMgenes<-readxl::read_xlsx("info_sheets/AM_genes.xlsx", sheet=1)
AMfaves<-readxl::read_xlsx("info_sheets/AM_genes.xlsx", sheet=2)
Yang2021.degs_flg22_overlap_up<-readxl::read_xlsx("info_sheets/Yang2021.degs_flg22_overlap_up.xlsx")
Yang2021.degs<-readxl::read_xlsx("info_sheets/Yang2021_genes.xlsx", sheet = 2)
Gutjahr2015_gse_deg<-read_tsv("info_sheets/Gutjahr2015_gse_deg.tsv")
GSRgenes.rice<-read_tsv("info_sheets/Bjornson2021_GSRgenes_rice.tsv")
Tang2021.flg22<-readxl::read_xlsx("info_sheets/Tang2021_genes.xlsx", sheet = 2)
Yang2021.degs<-readxl::read_xlsx("info_sheets/Yang2021_genes.xlsx", sheet = 2)
Tian2018_magna_up<-read_tsv("info_sheets/Tian2018_Mo_UP.tsv")
Tian2018_magna_up<-Tian2018_magna_up %>% filter(!is.na(RAP))
Das2022_WT_Myc_LP_vs_WT_Mock_LP<-read_tsv("info_sheets/Das2022_WT_Myc_LP_vs_WT_Mock_LP.tsv")

# Gene count tables
normalized_counts<-read.csv("gene_count_tables/DESeq2_normalized_counts.csv", row.names = 1)
vsd_counts<-read.csv("gene_count_tables/DESeq2_VSD_counts.csv", row.names = 1)
vsd_counts_for_plotting<- vsd_counts %>% rownames_to_column("RAP") %>% gather(name, value, as.vector(unique(sampleTable$name))) %>% as_tibble()
vsd_counts_for_plotting<-left_join(vsd_counts_for_plotting, sampleTable, by="name") 
normalized_counts_for_plotting<- normalized_counts %>% rownames_to_column("RAP") %>% gather(name, value, as.vector(unique(sampleTable$name))) %>% as_tibble()
normalized_counts_for_plotting<-left_join(normalized_counts_for_plotting, sampleTable, by="name") 

annotate_gene_list<-function(gene.list)
{
  table<-tibble("RAP"=gene.list)
  table <- table %>% left_join(., rice_annotation, by = "RAP") 
  return(table)
}

# DEG lists

DEG.table.names<-substr(list.files("DEG_tables"),1,nchar(list.files("DEG_tables"))-4)
DEG.list<-list()
for(i in 1:length(DEG.table.names))
{
  DEG.table<-read_tsv(paste("DEG_tables/",DEG.table.names[i],".tsv", sep=""), show_col_types = FALSE)
  DEG.list[[i]]<-c(DEG.table$RAP)
}
names(DEG.list)<-DEG.table.names

# Code to obtain a FC matrix

FC.table.names<-substr(list.files("comparison_tables"),17,nchar(list.files("comparison_tables"))-4)
FC.table<-tibble(RAP=rownames(normalized_counts))
for(i in 1:length(FC.table.names))
{
  comparison.table<-read_tsv(paste("comparison_tables/all_genes_table_",FC.table.names[i],".tsv", sep=""), show_col_types = FALSE)
  FC.table.i<-comparison.table %>% dplyr::select(RAP, log2FoldChange)
  colnames(FC.table.i) <- c("RAP", FC.table.names[i])
  FC.table<-left_join(FC.table, FC.table.i, by="RAP")
}
write_tsv(x = FC.table, file = "gene_count_tables/all_FC_table.tsv")

# Reading the FC matrix already provided
FC.table<-read_tsv("gene_count_tables/all_FC_table.tsv")
FC.df<-as.data.frame(FC.table)
FC.df<-column_to_rownames(FC.table, var = "RAP")
FC.df<-FC.df[rowSums(is.na(FC.df)) != ncol(FC.df), ]
FC_for_plotting<- FC.table %>% gather(comparison, FC, colnames(FC.df))
FC_for_plotting[is.na(FC_for_plotting)] <- 0

# Code to obtain a FDR matrix

padj.names<-substr(list.files("comparison_tables"),17,nchar(list.files("comparison_tables"))-4)
padj.table<-tibble(RAP=rownames(normalized_counts))
for(i in 1:length(padj.names))
{
  comparison.table<-read_tsv(paste("comparison_tables/all_genes_table_",padj.names[i],".tsv", sep=""), show_col_types = FALSE)
  padj.table.i<-comparison.table %>% dplyr::select(RAP, padj)
  colnames(padj.table.i) <- c("RAP", padj.names[i])
  padj.table<-left_join(padj.table, padj.table.i, by="RAP")
}
write_tsv(x = padj.table, file = "gene_count_tables/all_padj_table.tsv")

# Reading the FDR matrix already provided and generating significance matrix

padj.table<-read_tsv("gene_count_tables/all_padj_table.tsv")
padj.df<-as.data.frame(padj.table)
padj.df<-column_to_rownames(padj.table, var = "RAP")
padj.df<-padj.df[rowSums(is.na(padj.df)) != ncol(padj.df), ]
padj_for_plotting<- padj.table %>% gather(comparison, padj, colnames(padj.df))
padj_for_plotting[is.na(padj_for_plotting)] <- 1
padj.matrix<-as.matrix(padj.df)
padj_sign_matrix <- ifelse(
  !is.na(padj.matrix) & padj.matrix <= 0.001, "***",
  ifelse(
    !is.na(padj.matrix) & padj.matrix <= 0.01, "**",
    ifelse(
      !is.na(padj.matrix) & padj.matrix <= 0.05, "*", ""
    )
  )
)

# Count to create VST count matrix averaged by replicates:

vsd_counts_averages<-data.frame(RAP = rownames(vsd_counts))
for (i in unique(sampleTable$condition))
{
  print(i)
  samples <- filter(sampleTable, condition == i) %>% pull(name)
  print(samples)
  new_column <- vsd_counts %>% dplyr::select(all_of(samples)) %>% rowMeans(na.rm = T)
  vsd_counts_averages<-cbind(vsd_counts_averages, new_column)
}
colnames(vsd_counts_averages)<-c("RAP", as.vector(unique(sampleTable$condition)))
vsd_counts_averages<-as_tibble(vsd_counts_averages)
write.csv(vsd_counts_averages, "gene_count_tables/DESeq2_VSD_counts_averages.csv")

# Reading the VST count matrix already provided:
vsd_counts_averages<-read.csv("gene_count_tables/DESeq2_VSD_counts_averages.csv", row.names = 1)
vsd_counts_averages_for_plotting<-gather(vsd_counts_averages, condition, value, as.vector(unique(sampleTable$condition)))
sampleTable_reduced<-unique(dplyr::select(sampleTable, condition, genotype, treatment, genotype_name, treatment_name))
vsd_counts_averages_for_plotting<-left_join(vsd_counts_averages_for_plotting, sampleTable_reduced, by = "condition")
vsd_counts_averages<-column_to_rownames(as.data.frame(vsd_counts_averages), "RAP")

# PCA -------------------------------------------------------------------------
# Figure 2G
ntop <- 5000
rv <- rowVars(as.matrix(vsd_counts))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(vsd_counts[select,])
pca <- prcomp(mat)
summary.pca<-summary(pca)
pca.loadings<-unlist(as.vector(as.data.frame(summary.pca$importance)[2,]))

scores <- as.data.frame(pca$x)
scores <- as_tibble(rownames_to_column(scores, var = "name"))
scores <- left_join(scores[,1:5], sampleTable, by="name")

ggplot(filter(scores, genotype %in% c("A"), treatment %in% c("01","06", "07", "08", "09", "10","11","12","13")), 
       aes(PC1, PC2, color=treatment_name, shape=treatment_name)) +
  xlab(paste0("PC1: ",pca.loadings[1]*100,"% variance")) +
  ylab(paste0("PC2: ",pca.loadings[2]*100,"% variance")) + 
  scale_y_continuous(limits = c(-40,25))+
  scale_x_continuous(limits = c(-60,90))+
  ggforce::geom_mark_ellipse(aes(color=treatment_name, fill=treatment_name), alpha=0.1)+
  geom_point(alpha=1, size=3) +
  publication_theme_PCA()+
  scale_color_manual(name = "Treatment", values=treatment_colours)+
  scale_fill_manual(name = "Treatment", values=treatment_colours)+
  scale_shape_manual(name = "Treatment", values=treatment_shapes, guide = guide_legend(label.theme = element_text(angle = 0, size=10,face = "italic")))+
  geom_vline(xintercept=c(0), linetype="dashed")+
  geom_hline(yintercept=c(0), linetype="dashed")
ggsave(filename = "plots/PCA/figure_2G.pdf", width = 5, height = 4)

# Figure 4B
ntop <- 5000
rv <- genefilter::rowVars(vsd_counts_averages)
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t(vsd_counts_averages[select,])
pca <- prcomp(mat)
summary.pca<-summary(pca)
pca.loadings<-unlist(as.vector(as.data.frame(summary.pca$importance)[2,]))

scores <- as.data.frame(pca$x)
scores <- as_tibble(rownames_to_column(scores, var = "condition"))
scores <- left_join(scores, sampleTable_reduced, by="condition")

ggplot(filter(scores, genotype %in% c("A","B","C","F","G","H","I","J"), treatment %in% c("01","11","12","13")), 
       aes(PC1, PC2, color=treatment_name, shape=genotype_name)) +
  geom_vline(xintercept=c(0), linetype="dashed")+
  geom_hline(yintercept=c(0), linetype="dashed")+
  geom_point(alpha=1, size=3) +
  xlab(paste0("PC1: ",pca.loadings[1]*100,"% variance")) +
  ylab(paste0("PC2: ",pca.loadings[2]*100,"% variance")) + 
  geom_point(alpha=1, size=3) +
  scale_shape_manual(name = "Genotype", values=genotype_shapes, guide = guide_legend(label.theme = element_text(angle = 0, size=10,face = "italic")))+
  publication_theme_PCA()+
  scale_color_manual(name = "Treatment", values=treatment_colours)

ggsave(filename = "plots/PCA/figure_4B.pdf", width = 4.5, height = 3.5)


# Heatmaps ---------------------------------------------------------------------
# Figure 1C

filtered.module<-filter(AMfaves, RAP %in% row.names(FC.df))
significance_matrix<-padj_sign_matrix[filtered.module$RAP,
                                      c("A13_vs_A01","A12_vs_A01")]
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A13_vs_A01","A12_vs_A01")])
fc.das <- setNames(Das2022_WT_Myc_LP_vs_WT_Mock_LP$LFC, Das2022_WT_Myc_LP_vs_WT_Mock_LP$RAP)[rownames(FC.matrix)]
FC.matrix <- cbind(FC.matrix, WT_Myc_LP_vs_WT_Mock_LP = fc.das)
padj.das <- setNames(Das2022_WT_Myc_LP_vs_WT_Mock_LP$FDR, Das2022_WT_Myc_LP_vs_WT_Mock_LP$RAP)[rownames(FC.matrix)]
padj_sign.das <- ifelse(
  !is.na(padj.das) & padj.das <= 0.001, "***",
  ifelse(
    !is.na(padj.das) & padj.das <= 0.01, "**",
    ifelse(
      !is.na(padj.das) & padj.das <= 0.05, "*", ""
    )
  )
)
significance_matrix <- cbind(significance_matrix, WT_Myc_LP_vs_WT_Mock_LP = padj_sign.das)
FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("flg22","Free AMF, 4 wpi (this study)","Free AMF, 6 wpi (Das et al., 2022)")
col_heatmap = colorRamp2(c(floor(min(c(FC.matrix))),-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
pdf(file = "plots/heatmaps/figure_1C.pdf", width = 11, height = 2.5)
set.seed(123)
ComplexHeatmap::Heatmap(t(FC.matrix), 
                        name="log2(FC)", rect_gp = gpar(col = "white", lwd = 0.2), na_col = "grey",
                        
                        col=col_heatmap, 
                        
                        cluster_row = F, show_row_dend = FALSE, row_km = 0, cluster_row_slices = F,
                        row_split = c("A","B","C"),
                        
                        cluster_column = T, show_column_dend = T, column_km = 6, cluster_column_slices = T,

                        column_labels = paste(filtered.module$Name), 
                        column_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 1, fontface="bold"),
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 9), row_title_gp = gpar(fontsize = 0, fontface="bold"),
                        
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(significance_matrix[j, i],x = x+unit(0.8, "mm"), y = y,
                                    just = "center", rot = 90, gp = gpar(fontsize = 8.5, col = "red"))}
)
dev.off()

# Figure 2D
filtered.module<-annotate_gene_list(c(DEG.list$activated_genes_table_A06_vs_A01, DEG.list$activated_genes_table_A07_vs_A01,
                                      DEG.list$activated_genes_table_A08_vs_A01, DEG.list$activated_genes_table_A09_vs_A01,
                                      DEG.list$activated_genes_table_A10_vs_A01, DEG.list$activated_genes_table_A11_vs_A01,
                                      DEG.list$repressed_genes_table_A06_vs_A01, DEG.list$repressed_genes_table_A07_vs_A01,
                                      DEG.list$repressed_genes_table_A08_vs_A01, DEG.list$repressed_genes_table_A09_vs_A01,
                                      DEG.list$repressed_genes_table_A10_vs_A01, DEG.list$repressed_genes_table_A11_vs_A01))
filtered.module<-filter(filtered.module, RAP %in% row.names(FC.df))
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A06_vs_A01","A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A11_vs_A01","A12_vs_A01","A13_vs_A01")])
FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("Plant exudates","Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Encased AMF","Free AMF","flg22")
col_heatmap = colorRamp2(c(floor(min(c(FC.matrix))),-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
png(file = "plots/heatmaps/figure_2D.png", width = 10, height = 2, res = 600, units = "in")
pdf(file = "plots/heatmaps/figure_2D.pdf", width = 10, height = 2)
set.seed(123)
ComplexHeatmap::Heatmap(t(FC.matrix), 
                        name="log2(FC)",
                        col=col_heatmap, 
                        cluster_columns = T, cluster_rows = T,
                        cluster_row_slices = T, cluster_column_slices = T,
                        column_km = 7, row_km = 3,
                        show_row_names = T, show_column_names = F,
                        show_column_dend = F, show_row_dend = T, 
                        row_names_side = "left", row_dend_side = "right",
                        column_labels = paste(filtered.module$RAP),
                        row_names_gp = gpar(fontsize = 9),
                        na_col = "grey",
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        column_names_gp = gpar(fontsize = 8),
                        column_title_gp = gpar(fontsize = 1, fontface="bold"),
                        row_title_gp = gpar(fontsize = 0.1, fontface="bold"),
)
dev.off()

# Figure 3A

filtered.module<-filter(AMfaves, RAP %in% row.names(FC.df))
significance_matrix<-padj_sign_matrix[filtered.module$RAP,
                                      c("A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A06_vs_A01","A11_vs_A01","A13_vs_A01","A12_vs_A01")]
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A06_vs_A01","A11_vs_A01","A13_vs_A01","A12_vs_A01")])

FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Plant exudates","Encased AMF","flg22","Free AMF")
col_heatmap = colorRamp2(c(floor(min(c(FC.matrix))),-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
pdf(file = "plots/heatmaps/figure_3A.pdf", width = 11, height = 3.2)
set.seed(1)
ComplexHeatmap::Heatmap(t(FC.matrix), 
                        name="log2(FC)", rect_gp = gpar(col = "white", lwd = 0.2), na_col = "grey",
                        
                        col=col_heatmap, 
                        
                        cluster_row = F, show_row_dend = FALSE, row_km = 0, cluster_row_slices = F,
                        row_split = c(rep("A",6),"B","C"),
                        
                        cluster_column = T, show_column_dend = T, column_km = 5, cluster_column_slices = T,
                        
                        column_labels = paste(filtered.module$Name), 
                        column_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 1, fontface="bold"),
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 9), row_title_gp = gpar(fontsize = 0, fontface="bold"),
                        
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(significance_matrix[j, i],x = x+unit(0.8, "mm"), y = y,
                                    just = "center", rot = 90, gp = gpar(fontsize = 8.5, col = "red"))}
)
dev.off()

# Figure 4D
scaled_counts<-scale(t(as.matrix(vsd_counts_averages[c(DEG.list$activated_genes_table_A11_vs_A01, 
                                                       DEG.list$activated_genes_table_A12_vs_A01),
                                                     c("A01","A13","A11","A12",
                                                       "B01","I01","C01","J01","F01","G01","H01",
                                                       "B11","I11","C11","J11","F11","G11","H11",
                                                       "B12","I12","C12","J12","F12","G12","H12")])), center=TRUE, scale=TRUE)

png(file = "plots/heatmaps/figure_3D.png", width = 7, height = 3.5, res = 600, units = "in")
pdf(file = "plots/heatmaps/figure_3D.pdf", width = 7, height = 3.5)
set.seed(123)
ComplexHeatmap::Heatmap(scaled_counts, use_raster = F, name="Z-score",
                        col=colorRamp2(c(-2,0,2), c("blue", "black", "yellow")),
                        cluster_columns = T, cluster_rows = F, column_km =3,
                        row_split = c("A", "A","A","A", rep("B",7), rep("C",7),rep("D",7)),
                        show_column_dend = FALSE, show_row_dend = FALSE, show_column_names = F,
                        column_names_gp = gpar(fontsize = 8),
                        row_names_gp = gpar(fontsize = 8),
                        column_title_gp = gpar(fontsize = 0.1),
                        row_title_gp = gpar(fontsize = 0.1))
dev.off()

# Figure 4E
filtered.module<-filter(AMfaves, RAP %in% row.names(FC.df))
significance_matrix<-padj_sign_matrix[filtered.module$RAP,
                                      c("A11_vs_A01","B11_vs_B01","I11_vs_I01","C11_vs_C01","J11_vs_J01","F11_vs_F01","G11_vs_G01","H11_vs_H01",
                                        "A12_vs_A01","B12_vs_B01","I12_vs_I01","C12_vs_C01","J12_vs_J01","F12_vs_F01","G12_vs_G01","H12_vs_H01",
                                        "A13_vs_A01")]
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A11_vs_A01","B11_vs_B01","I11_vs_I01","C11_vs_C01","J11_vs_J01","F11_vs_F01","G11_vs_G01","H11_vs_H01",
                             "A12_vs_A01","B12_vs_B01","I12_vs_I01","C12_vs_C01","J12_vs_J01","F12_vs_F01","G12_vs_G01","H12_vs_H01",
                             "A13_vs_A01")])
FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("WT, Encased AMF","d14l, Encased AMF","cerk1, Encased AMF","cerk1/2, Encased AMF","nfr5, Encased AMF",
                       "pollux, Encased AMF","ccamk, Encased AMF","cyclops, Encased AMF",
                       "WT, Free AMF","d14l, Free AMF","cerk1, Free AMF","cerk1/2, Free AMF","nfr5, Free AMF",
                       "pollux, Free AMF","ccamk, Free AMF","cyclops, Free AMF",
                       "WT, flg22")
col_heatmap = colorRamp2(c(floor(min(c(FC.matrix))),-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
pdf(file = "plots/heatmaps/figure_3E.pdf", width = 11, height = 4.6)
set.seed(123)
ComplexHeatmap::Heatmap(t(FC.matrix), 
                        name="log2(FC)", rect_gp = gpar(col = "white", lwd = 0.2), na_col = "grey",
                        
                        col=col_heatmap, 
                        
                        cluster_row = F, show_row_dend = FALSE, row_km = 0, cluster_row_slices = F,
                        row_split = c(rep("Encased AMF",8),rep("Free AMF",8),"Z"),
                        
                        cluster_column = T, show_column_dend = T, column_km = 5, cluster_column_slices = T,
                        
                        column_labels = paste(filtered.module$Name), 
                        column_names_gp = gpar(fontsize = 8), column_title_gp = gpar(fontsize = 1, fontface="bold"),
                        
                        row_labels = c("WT","d14l","cerk1","cerk1/2","nfr5",
                                       "pollux","ccamk","cyclops",
                                       "WT","d14l","cerk1","cerk1/2","nfr5",
                                       "pollux","ccamk","cyclops",
                                       "WT, flg22"),
                        row_names_side = "left",
                        row_names_gp = gpar(fontsize = 9, fontface="italic"), row_title_gp = gpar(fontsize = 9),
                        
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(significance_matrix[j, i],x = x+unit(0.8, "mm"), y = y,
                                    just = "center", rot = 90, gp = gpar(fontsize = 8.5, col = "red"))}
)
dev.off()

# Figure S2B
filtered.module<-filter(AMfaves, RAP %in% row.names(FC.df))
significance_matrix<-padj_sign_matrix[filtered.module$RAP,
                                      c("A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A06_vs_A01","A11_vs_A01","A12_vs_A01","A13_vs_A01")]
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A06_vs_A01","A11_vs_A01","A12_vs_A01","A13_vs_A01")])
fc.das <- setNames(Das2022_WT_Myc_LP_vs_WT_Mock_LP$LFC, Das2022_WT_Myc_LP_vs_WT_Mock_LP$RAP)[rownames(FC.matrix)]
FC.matrix <- cbind(FC.matrix, WT_Myc_LP_vs_WT_Mock_LP = fc.das)
padj.das <- setNames(Das2022_WT_Myc_LP_vs_WT_Mock_LP$FDR, Das2022_WT_Myc_LP_vs_WT_Mock_LP$RAP)[rownames(FC.matrix)]
padj_sign.das <- ifelse(
  !is.na(padj.das) & padj.das <= 0.001, "***",
  ifelse(
    !is.na(padj.das) & padj.das <= 0.01, "**",
    ifelse(
      !is.na(padj.das) & padj.das <= 0.05, "*", ""
    )
  )
)
significance_matrix <- cbind(significance_matrix, WT_Myc_LP_vs_WT_Mock_LP = padj_sign.das)
FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Plant exudates","Encased AMF","Free AMF","flg22","Free AMF (Das et al. 2022)")
col_heatmap = colorRamp2(c(floor(min(c(FC.matrix))),-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
pdf(file = "plots/heatmaps/figure_S2B.pdf", width = 5, height = 11)
set.seed(123)
ComplexHeatmap::Heatmap(FC.matrix, 
                        name="log2(FC)", rect_gp = gpar(col = "white", lwd = 0.2),
                        col=col_heatmap, cluster_columns = F, row_km = 6,
                        cluster_rows = T, show_column_dend = FALSE, show_row_dend = T, 
                        row_labels = paste(filtered.module$RAP,"  ",filtered.module$Name),
                        column_split = c(rep("Exudates",6),"Free AMF","flg22","Free AMF (Das et al. 2022)"),
                        column_names_gp = gpar(fontsize = 9),
                        na_col = "grey",
                        row_names_gp = gpar(fontsize = 8),
                        row_title_gp = gpar(fontsize = 1, fontface="bold"),
                        column_title_gp = gpar(fontsize = 0.1, fontface="bold"),
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(significance_matrix[i, j],x = x,y = y - unit(0.6, "mm"),just = "center",gp = gpar(fontsize = 10, col = "red"))}
)
dev.off()

# Figure S3B
filtered.module<-filter(Yang2021.degs_flg22_overlap_up, RAP %in% row.names(FC.df), Name != "NA")
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A06_vs_A01","A11_vs_A01","A12_vs_A01","A13_vs_A01")])
significance_matrix<-padj_sign_matrix[filtered.module$RAP,
                                      c("A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A06_vs_A01","A11_vs_A01","A12_vs_A01","A13_vs_A01")]
new_column <- setNames(Yang2021.degs$fc, Yang2021.degs$RAP)
FC.matrix <- cbind(FC.matrix, Yang2021 = new_column[rownames(FC.matrix)])
significance_matrix <- cbind(significance_matrix, Yang2021 = rep("*", nrow(significance_matrix)))
FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Plant exudates","Encased AMF","Free AMF","flg22","M. oryzae (Yang et al. 2021)")
colnames(significance_matrix)<-c("Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Plant exudates","Encased AMF","Free AMF","flg22","M. oryzae (Yang et al. 2021)")
col_heatmap = colorRamp2(c(-3,-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
pdf(file = "plots/heatmaps/figure_S3B.pdf", width = 4.5, height = 12)
set.seed(122)
ComplexHeatmap::Heatmap(FC.matrix, 
                        name="log2(FC)",rect_gp = gpar(col = "white", lwd = 0.2),
                        col=col_heatmap, cluster_columns = F, row_km = 4,
                        cluster_rows = T, show_column_dend = FALSE, show_row_dend = T, 
                        row_labels = paste(filtered.module$RAP,"  ",filtered.module$Name),
                        column_split = c(rep("Exudates",6),"Free AMF","flg22","Mo"),
                        column_names_gp = gpar(fontsize = 9),
                        na_col = "grey",
                        row_names_gp = gpar(fontsize = 7),
                        row_title_gp = gpar(fontsize = 0.1, fontface="bold"),
                        column_title_gp = gpar(fontsize = 0.1, fontface="bold"),
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(significance_matrix[i, j],x = x,y = y - unit(0.6, "mm"),just = "center",gp = gpar(fontsize = 10, col = "red"))}
                        )
dev.off()

# Figure S4C
filtered.module<-filter(Gutjahr2015_gse_deg, RAP %in% row.names(FC.df))
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01","A06_vs_A01","A11_vs_A01","A12_vs_A01","A13_vs_A01")])
new_column <- setNames(Gutjahr2015_gse_deg$Log2FC_WT.M_vs_WT.G, Gutjahr2015_gse_deg$RAP)
FC.matrix <- cbind(FC.matrix, RiGSE_Gutjahr2015 = new_column[rownames(FC.matrix)])
FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Plant exudates","Encased AMF","Free AMF","flg22","Ri GSE (Gutjahr et al., 2015)")
col_heatmap = colorRamp2(c(floor(min(c(FC.matrix))),-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
pdf(file = "plots/heatmaps/figure_S4C.pdf", width = 2.5, height = 6)
set.seed(123)
ComplexHeatmap::Heatmap(FC.matrix, 
                        name="log2(FC)", rect_gp = gpar(col = "white", lwd = 0.2),
                        col=col_heatmap, cluster_columns = F, show_row_names = F,
                        cluster_rows = T, show_column_dend = FALSE, show_row_dend = T, 
                        row_labels = paste(filtered.module$RAP),
                        row_split = filtered.module$`Up(1)_Down(-1)_WT.M_vs_WT.G`,
                        column_split = c(rep("Exudates",6),"Free AMF","flg22","2015"),
                        column_names_gp = gpar(fontsize = 9),
                        na_col = "grey",
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        row_names_gp = gpar(fontsize = 7),
                        row_title_gp = gpar(fontsize = 0.1, fontface="bold"),
                        column_title_gp = gpar(fontsize = 0.1, fontface="bold"))
dev.off()

# Figure S7
filtered.module<-filter(Yang2021.degs_flg22_overlap_up, RAP %in% row.names(FC.df), Name != "NA")
significance_matrix<-padj_sign_matrix[filtered.module$RAP,
                                      c("A11_vs_A01","B11_vs_B01","I11_vs_I01","C11_vs_C01","J11_vs_J01","F11_vs_F01","G11_vs_G01","H11_vs_H01",
                                        "A12_vs_A01","B12_vs_B01","I12_vs_I01","C12_vs_C01","J12_vs_J01","F12_vs_F01","G12_vs_G01","H12_vs_H01",
                                        "A13_vs_A01")]
FC.matrix<-as.matrix(FC.df[filtered.module$RAP,
                           c("A11_vs_A01","B11_vs_B01","I11_vs_I01","C11_vs_C01","J11_vs_J01","F11_vs_F01","G11_vs_G01","H11_vs_H01",
                             "A12_vs_A01","B12_vs_B01","I12_vs_I01","C12_vs_C01","J12_vs_J01","F12_vs_F01","G12_vs_G01","H12_vs_H01",
                             "A13_vs_A01")])
FC.matrix[is.na(FC.matrix)] <- 0
colnames(FC.matrix)<-c("WT, Encased AMF","d14l, Encased AMF","cerk1, Encased AMF","cerk1/2, Encased AMF","nfr5, Encased AMF",
                       "pollux, Encased AMF","ccamk, Encased AMF","cyclops, Encased AMF",
                       "WT, Free AMF","d14l, Free AMF","cerk1, Free AMF","cerk1/2, Free AMF","nfr5, Free AMF",
                       "pollux, Free AMF","ccamk, Free AMF","cyclops, Free AMF",
                       "WT, flg22")
col_heatmap = colorRamp2(c(floor(min(c(FC.matrix))),-2,-1,0,1,2,ceiling(max(FC.matrix))), 
                         c('#B6B6FF',"#0000FF",'#000048',"black",'#484800',"#FFFF00",'#FFFFB6'))
pdf(file = "plots/heatmaps/figure_S7.pdf", width = 6.5, height = 13)
set.seed(123)
ComplexHeatmap::Heatmap(FC.matrix, 
                        name="log2(FC)", rect_gp = gpar(col = "white", lwd = 0.2),
                        col=col_heatmap, cluster_columns = F,row_km=3,
                        cluster_rows = T, show_column_dend = FALSE, show_row_dend = T, 
                        row_labels = paste(filtered.module$RAP,"  ",filtered.module$Name),
                        column_split = c(rep("1",8),rep("2",8),"3"),
                        column_names_gp = gpar(fontsize = 9, fontface="italic"),
                        na_col = "grey",
                        row_names_gp = gpar(fontsize = 8),
                        row_title_gp = gpar(fontsize = 0.1, fontface="bold"),
                        column_title_gp = gpar(fontsize = 0.1, fontface="bold"),
                        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface="bold"),labels_gp = gpar(fontsize = 8),
                                                    legend_height = unit(2, "cm"),legend_width = unit(1, "cm")),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(significance_matrix[i, j],x = x,y = y - unit(0.6, "mm"),just = "center",gp = gpar(fontsize = 10, col = "red"))}
                        )
dev.off()

# Average expression plot ------------------------------------------------------
# Figure 3C

vsd_counts_for_plot <- filter(vsd_counts_for_plotting, RAP %in% GSRgenes.rice$RAP, 
                              genotype %in% c("A"), 
                              treatment %in% c("01","06","07","08","09","10","11","12","13"))

# Calculate means for each group to position the asterisks
mean_values <- vsd_counts_for_plot %>%
  group_by(treatment_name) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

# Perform pairwise t-tests against H2O reference
t_test_results <- compare_means(
  value ~ treatment_name, 
  data = vsd_counts_for_plot,
  ref.group = "H2O",  # Make sure this matches your reference group name
  method = "t.test"
) %>%
  # Add mean values for positioning
  left_join(mean_values, by = c("group2" = "treatment_name")) %>%
  # Create significance labels
  mutate(signif = case_when(
    p < 0.001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  ))

ggplot(vsd_counts_for_plot, aes(x = treatment_name, y = value)) +
  stat_summary(fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95), 
               geom = "ribbon", alpha = 0.5, aes(group = 1), color = "white", fill = "grey80") +
  stat_summary(fun.data = mean_se, geom = "pointrange", size = 0.2) +
  # Add ANOVA p-value
  stat_compare_means(label.y = 7.8, label.x = 1.5, size = 3, method = "anova") + 
  # Add t-test results as text - positioned 0.01 above the mean
  geom_text(data = t_test_results, 
            aes(x = group2, y = mean_value + 0.05, label = signif),
            vjust = 0, size = 4) +
  geom_hline(yintercept = mean(vsd_counts_for_plot$value[vsd_counts_for_plot$condition == "A01"]), 
             linetype = "dashed", colour = "black", alpha = 0.4, size = 0.4) +
  scale_y_continuous(name = "VST counts for GSR genes") +
  scale_x_discrete(name = "") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 10, color = "white"),
        strip.background = element_rect(color = "black", fill = "black", linewidth = 1.5, linetype = "solid"),
        legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.3, hjust = 1, size = 10, color = "black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, color = "black"),
        axis.title.y = element_text(face = "bold", color = "black", size = 10),
        axis.title.x = element_text(face = "bold", color = "black", size = 10),
        plot.title = element_text(face = "bold", color = "black", hjust = 0.5))

ggsave(filename = paste("plots/Average_expression_plots/figure_3C.pdf", sep = ""), width = 3, height = 4)

# Figure S8A
vsd_counts_for_plot<-filter(vsd_counts_for_plotting, RAP %in% GSRgenes.rice$RAP)
vsd_counts_for_plot<-filter(vsd_counts_for_plot, genotype %in% c("A","B","C","J","F","G","H","I"),treatment %in% c("01","11","12"))
vsd_counts_for_plot$condition<-factor(vsd_counts_for_plot$condition,levels=c("A01","B01","I01","C01","J01","F01","G01","H01",
                                                                             "A11","B11","I11","C11","J11","F11","G11","H11",
                                                                             "A12","B12","I12","C12","J12","F12","G12","H12"))

# Calculate means for each group to position the asterisks
mean_values <- vsd_counts_for_plot %>%
  group_by(condition) %>%
  summarise(mean_value = mean(value, na.rm = TRUE))

# Perform pairwise t-tests against H2O reference
t_test_results <- compare_means(
  value ~ condition, 
  data = vsd_counts_for_plot,
  ref.group = "A01",  # Make sure this matches your reference group name
  method = "t.test"
) %>%
  # Add mean values for positioning
  left_join(mean_values, by = c("group2" = "condition")) %>%
  # Create significance labels
  mutate(signif = case_when(
    p < 0.001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  ))

ggplot(vsd_counts_for_plot, aes(x=condition, y=value))+
  stat_summary(fun.data = mean_cl_normal, fun.args=list(conf.int=0.95), geom="ribbon", alpha=0.5, aes(group=1), color="white", fill="grey80")+
  stat_summary(fun.data=mean_se, geom="pointrange", size=0.2, color="black")+
  stat_compare_means(label.y = 7.75, label.x = 0.7, hjust=0, size=3, method="anova")+ 
  geom_text(data = t_test_results, 
            aes(x = group2, y = mean_value + 0.05, label = signif),
            vjust = 0, size = 4) +
  scale_y_continuous(name="VST counts for GSR genes")+
  scale_x_discrete(name="", labels=c("WT","d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops",
                                     "WT","d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops",
                                     "WT","d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops"))+
  geom_vline(xintercept = 8.5, color = "black", linewidth=0.4)+
  geom_vline(xintercept = 16.5, color = "black", linewidth=0.4)+
  annotate(geom="text", x = 4.5, y = 7, label="H2O", fontface="bold", size=3)+
  annotate(geom="text", x = 12.5, y = 7, label="Encased AMF", fontface="bold", size=3)+
  annotate(geom="text", x = 20.5, y = 7, label="Free AMF", fontface="bold", size=3)+
  publication_theme_average()+
  geom_hline(yintercept=mean(vsd_counts_for_plot $ value[vsd_counts_for_plot$condition == "A01"]), linetype="dashed", colour = "black", alpha=0.4, size=0.4)

ggsave(filename = paste("plots/Average_expression_plots/figure_S8A.pdf",sep=""), width = 7, height = 3)

# Venn -------------------------------------------------------------------------
# Figure 1A
colors.venn<-c("skyblue", "dodgerblue", "palegreen")
venn.list<-list("Free AMF (this study)"=c(DEG.list$activated_genes_table_A12_vs_A01),
                "Free AMF (Das et al. 2022)"=filter(Das2022_WT_Myc_LP_vs_WT_Mock_LP, FDR < 0.05, LFC > log2(1.5))$RAP,
                "AM gene list"=unique(AMgenes$RAP))

ggvenn(venn.list,
       fill_color = colors.venn, 
       stroke_size = 0.5, stroke_color = "white",
       set_name_color = colors.venn, 
       text_size = 6, show_percentage = F,
       set_name_size = 5)+scale_x_continuous(expand = expansion(mult = .2))
ggsave("plots/Venn/figure_1A.pdf", height=4, width=4)

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
view(filter(result_tibble,  ! is.na(AM_gene_description)))
write_tsv(result_tibble, "plots/Venn/figure_1A.tsv")

# Figure 1B
colors.venn<-c("magenta", "skyblue", "palegreen")
venn.list<-list("flg22 DOWN"=c(DEG.list$repressed_genes_table_A13_vs_A01),
                "Free AMF UP"=c(DEG.list$activated_genes_table_A12_vs_A01),
                "AM gene list"=unique(AMgenes$RAP))
ggvenn(venn.list,
       fill_color = colors.venn, 
       stroke_size = 0.5, stroke_color = "white",
       set_name_color = colors.venn, 
       text_size = 6, show_percentage = F,
       set_name_size = 5)+scale_x_continuous(expand = expansion(mult = .2))
ggsave("plots/Venn/figure_1B.pdf", height=4, width=4)

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
view(filter(result_tibble,  ! is.na(AM_gene_description)))
write_tsv(result_tibble, "plots/Venn/figure_1B.tsv")

# Figure 2A
colors.venn<-c("palegreen", "seagreen3", "dodgerblue")
venn.list<-list("Ri naïve GSE"=c(DEG.list$repressed_genes_table_A07_vs_A01,DEG.list$activated_genes_table_A07_vs_A01),
                "Mo naïve GSE"=c(DEG.list$repressed_genes_table_A09_vs_A01,DEG.list$activated_genes_table_A09_vs_A01),
                "Free AMF"=c(DEG.list$repressed_genes_table_A12_vs_A01,DEG.list$activated_genes_table_A12_vs_A01))
venn_fit <- euler(venn.list)
pdf("plots/Venn/figure_2A.pdf", height=4, width=4)
plot(venn_fit,
     fills = list(fill = colors.venn, alpha = 0.6),
     edges = list(col = "white"),
     labels = list(col = NA),
     #labels = list(col = colors.venn, font = 2, cex = 1.2),
     quantities = TRUE, main = NULL, legend = T)
dev.off()

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
write_tsv(result_tibble, "plots/Venn/figure_2A.tsv")

# Figure 2B
colors.venn<-c("palegreen", "seagreen3", "magenta")
venn.list<-list("Ri naïve GSE"=c(DEG.list$repressed_genes_table_A07_vs_A01,DEG.list$activated_genes_table_A07_vs_A01),
                "Mo naïve GSE"=c(DEG.list$repressed_genes_table_A09_vs_A01,DEG.list$activated_genes_table_A09_vs_A01),
                "flg22"=c(DEG.list$repressed_genes_table_A13_vs_A01,DEG.list$activated_genes_table_A13_vs_A01))
venn_fit <- euler(venn.list)
pdf("plots/Venn/figure_2B.pdf", height=4, width=4)
plot(venn_fit,
     fills = list(fill = colors.venn, alpha = 0.6),
     edges = list(col = "white"),
     labels = list(col = NA),
     #labels = list(col = colors.venn, font = 2, cex = 1.2),
     quantities = TRUE, main = NULL, legend = T)
dev.off()

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
write_tsv(result_tibble, "plots/Venn/figure_2B.tsv")

# Figure 2E
colors.venn<-c("goldenrod1", "palegreen", "seagreen3", "dodgerblue")
venn.list<-list("Encased AMF"=c(DEG.list$repressed_genes_table_A11_vs_A01,DEG.list$activated_genes_table_A11_vs_A01),
                "Plant exudates"=c(DEG.list$repressed_genes_table_A06_vs_A01,DEG.list$activated_genes_table_A06_vs_A01),
                "Ri primed GSE"=c(DEG.list$repressed_genes_table_A08_vs_A01,DEG.list$activated_genes_table_A08_vs_A01),
                "Free AMF"=c(DEG.list$repressed_genes_table_A12_vs_A01,DEG.list$activated_genes_table_A12_vs_A01))
ggvenn(venn.list,
       fill_color = colors.venn, 
       stroke_size = 0.5, stroke_color = "white",
       set_name_color = colors.venn, 
       text_size = 6, show_percentage = F,
       set_name_size = 5)+scale_x_continuous(expand = expansion(mult = .2))
ggsave("plots/Venn/figure_2E.pdf", height=5, width=5)

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
write_tsv(result_tibble, "plots/Venn/figure_2E.tsv")

# Figure S3A
colors.venn<-c("magenta","darkorchid1", "yellow", "goldenrod1")
venn.list<-list("flg22, 6h, roots (this study)"=c(DEG.list$activated_genes_table_A13_vs_A01),
                "flg22, 3h, seedlings (Tang 2021)"=filter(Tang2021.flg22, log2FoldChange.flg22 > 0)$RAP,
                "M. oryzae infection, roots (Tian 2019)"=c(Tian2018_magna_up$RAP),
                "M. oryzae infection, leaves (Yang 2021)"=filter(Yang2021.degs, fc > 0)$RAP)
ggvenn(venn.list,
       fill_color = colors.venn, 
       stroke_size = 0.5, stroke_color = "white",
       set_name_color = "white", 
       text_size = 6, show_percentage = F,
       set_name_size = 5)+scale_x_continuous(expand = expansion(mult = .2))
ggsave("plots/Venn/figure_S3A.pdf", height=5, width=5)

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
write_tsv(result_tibble, "plots/Venn/figure_S3A.tsv")

# Figure S3D
colors.venn<-c("magenta", "dodgerblue", "skyblue")
venn.list<-list("flg22 DOWN"=c(DEG.list$repressed_genes_table_A13_vs_A01),
                "Free AMF (Das et al. 2022)"=filter(Das2022_WT_Myc_LP_vs_WT_Mock_LP, FDR < 0.05, LFC > log2(1.5))$RAP,
                "AM gene list"=unique(AMgenes$RAP))

ggvenn(venn.list,
       fill_color = colors.venn, 
       stroke_size = 0.5, stroke_color = "white",
       set_name_color = colors.venn, 
       text_size = 6, show_percentage = F,
       set_name_size = 5)+scale_x_continuous(expand = expansion(mult = .2))
ggsave("plots/Venn/figure_S3D.pdf", height=4, width=4)

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
view(filter(result_tibble,  ! is.na(AM_gene_description)))
write_tsv(result_tibble, "plots/Venn/figure_S3D.tsv")

# Figure S5A
colors.venn<-c("dodgerblue", "green", "goldenrod1")
venn.list<-list("Ri naïve GSE (Gutjahr 2015)"=Gutjahr2015_gse_deg$RAP,
                "Ri naïve GSE"=c(DEG.list$repressed_genes_table_A07_vs_A01,DEG.list$activated_genes_table_A07_vs_A01),
                "Mo naïve GSE"=c(DEG.list$repressed_genes_table_A09_vs_A01,DEG.list$activated_genes_table_A09_vs_A01)
)
ggvenn(venn.list,
       fill_color = colors.venn, 
       stroke_size = 0.5, stroke_color = "white",
       set_name_color = colors.venn, 
       text_size = 6, show_percentage = F,
       set_name_size = 5)+scale_x_continuous(expand = expansion(mult = .2))
ggsave("plots/Venn/figure_S5A.pdf", height=4, width=4)

result_tibble <- tibble(RAP = unique(unlist(venn.list)))
for (set_name in names(venn.list)) {result_tibble[[set_name]] <- ifelse(result_tibble$RAP %in% venn.list[[set_name]], 1, 0)}
result_tibble<-result_tibble %>% left_join(dplyr::select(rice_annotation, RAP, MSU, Name, AM_gene_description), by = "RAP")
write_tsv(result_tibble, "plots/Venn/figure_S5A.tsv")

# DEG number plot --------------------------------------------------------------
DEG.length.table<-tibble("name"=names(DEG.list), "length"=as.integer(summary(DEG.list)[,1]))
DEG.length.table<-mutate(DEG.length.table, "comparison"=substr(name, 23, nchar(name)))
DEG.length.table<-mutate(DEG.length.table, "up.or.down"=substr(name, 1, 9))
DEG.length.table<-mutate(DEG.length.table, "length.adjusted"=length)
for(i in 1:nrow(DEG.length.table))
{
  if(DEG.length.table$up.or.down[i] == "repressed")
  {
    DEG.length.table$length.adjusted[i]<-(-DEG.length.table$length.adjusted[i])
  }
}
# Figure 1C

DEG.length.table.1 <- filter(DEG.length.table, comparison %in% c("A07_vs_A01",
                                                                 "A08_vs_A01",
                                                                 "A09_vs_A01",
                                                                 "A10_vs_A01",
                                                                 "A06_vs_A01",
                                                                 "A11_vs_A01",
                                                                 "A12_vs_A01",
                                                                 "A13_vs_A01",
                                                                 "A07_vs_A06",
                                                                 "A08_vs_A06",
                                                                 "A09_vs_A06",
                                                                 "A10_vs_A06",
                                                                 "A08_vs_A07",
                                                                 "A10_vs_A09"))

DEG.length.table.1$comparison<-factor(DEG.length.table.1$comparison, levels=c("A07_vs_A01",
                                                                              "A08_vs_A01",
                                                                              "A09_vs_A01",
                                                                              "A10_vs_A01",
                                                                              "A06_vs_A01",
                                                                              "A11_vs_A01",
                                                                              "A12_vs_A01",
                                                                              "A13_vs_A01",
                                                                              "A07_vs_A06",
                                                                              "A08_vs_A06",
                                                                              "A09_vs_A06",
                                                                              "A10_vs_A06",
                                                                              "A08_vs_A07",
                                                                              "A10_vs_A09"))

ggplot(DEG.length.table.1, aes(x=comparison, y=length.adjusted, fill=up.or.down))+
  geom_bar(stat="identity", color="black", width=0.7)+
  geom_text(aes(y =(length.adjusted+300*sign(length.adjusted)), label = length), size = 3)+
  theme(    panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle=90, vjust=0.3, hjust=1, size=10, color="black"),
            axis.ticks.x = element_line(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.text.y = element_text(angle=0, vjust=0.3, size=10, color="black"),
            axis.title.y = element_text(face="bold", color="black", size=10),
            axis.title.x = element_text(face="bold", color="black", size=10),
            plot.title = element_text (face="bold", color="black", hjust = 0.5))+
  geom_vline(xintercept = 8.5, color = "black", linewidth=0.5)+
  geom_vline(xintercept = 12.5, color = "black", linewidth=0.5)+
  annotate(geom="text", x = 4.5, y = -3000, label="vs H2O", fontface="bold", size=3)+
  annotate(geom="text", x = 10.5, y = -3000, label="vs Plant exudates", fontface="bold", size=3)+
  annotate(geom="text", x = 13.5, y = -3000, label="vs naïve", fontface="bold", size=3)+
  annotate(geom="text", x = 0.5, y = 3000, label="UP", fontface="bold", size=3, color="yellow4", hjust=0)+
  annotate(geom="text", x = 0.5, y = -3000, label="DOWN", fontface="bold", size=3, color="blue",hjust=0)+
  scale_fill_manual(values=c("yellow", "blue"), name="")+
  scale_x_discrete(name="", labels=c("Ri naïve GSE", "Ri primed GSE","Mo naïve GSE", "Mo primed GSE","Plant exudates", "Encased AMF","Free AMF","flg22",
                                     "Ri naïve GSE", "Ri primed GSE","Mo naïve GSE", "Mo primed GSE","Ri primed GSE", "Mo primed GSE"))+
  scale_y_continuous(name="Number of DEGs", breaks = c(-3000,-2000,-1000,0,1000,2000,3000), limits=c(-3000,3000))+ 
  geom_hline(yintercept=0, color = "black", linewidth=0.4)
ggsave("plots/DEG_number/figure_1C.pdf", height=3.5, width=5)


# Figure S6

DEG.length.table.1 <- filter(DEG.length.table, comparison %in% c("A12_vs_A01",
                                                                 "B12_vs_B01",
                                                                 "I12_vs_I01",
                                                                 "C12_vs_C01",
                                                                 "J12_vs_J01",
                                                                 "F12_vs_F01",
                                                                 "G12_vs_G01",
                                                                 "H12_vs_H01",
                                                                 
                                                                 "A11_vs_A01",
                                                                 "B11_vs_B01",
                                                                 "I11_vs_I01",
                                                                 "C11_vs_C01",
                                                                 "J11_vs_J01",
                                                                 "F11_vs_F01",
                                                                 "G11_vs_G01",
                                                                 "H11_vs_H01",
                                                                 
                                                                 "B12_vs_A12",
                                                                 "I12_vs_A12",
                                                                 "C12_vs_A12",
                                                                 "J12_vs_A12",
                                                                 "F12_vs_A12",
                                                                 "G12_vs_A12",
                                                                 "H12_vs_A12",
                                                                 
                                                                 "B11_vs_A11",
                                                                 "I11_vs_A11",
                                                                 "C11_vs_A11",
                                                                 "J11_vs_A11",
                                                                 "F11_vs_A11",
                                                                 "G11_vs_A11",
                                                                 "H11_vs_A11"
                                                                 ))

DEG.length.table.1$comparison<-factor(DEG.length.table.1$comparison, levels=c("A12_vs_A01",
                                                                              "B12_vs_B01",
                                                                              "I12_vs_I01",
                                                                              "C12_vs_C01",
                                                                              "J12_vs_J01",
                                                                              "F12_vs_F01",
                                                                              "G12_vs_G01",
                                                                              "H12_vs_H01",
                                                                              
                                                                              "A11_vs_A01",
                                                                              "B11_vs_B01",
                                                                              "I11_vs_I01",
                                                                              "C11_vs_C01",
                                                                              "J11_vs_J01",
                                                                              "F11_vs_F01",
                                                                              "G11_vs_G01",
                                                                              "H11_vs_H01",
                                                                              
                                                                              "B12_vs_A12",
                                                                              "I12_vs_A12",
                                                                              "C12_vs_A12",
                                                                              "J12_vs_A12",
                                                                              "F12_vs_A12",
                                                                              "G12_vs_A12",
                                                                              "H12_vs_A12",
                                                                              
                                                                              "B11_vs_A11",
                                                                              "I11_vs_A11",
                                                                              "C11_vs_A11",
                                                                              "J11_vs_A11",
                                                                              "F11_vs_A11",
                                                                              "G11_vs_A11",
                                                                              "H11_vs_A11"
))


ggplot(DEG.length.table.1, aes(x=comparison, y=length.adjusted, fill=up.or.down))+
  geom_bar(stat="identity", color="black", width=0.7)+
  geom_text(aes(y =(length.adjusted+200*sign(length.adjusted)), label = length), size = 3)+
  theme(    panel.background = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.text.x = element_text(angle=90, vjust=0.3, hjust=1, size=10, color="black",face="italic"),
            axis.ticks.x = element_line(color = "black"),
            axis.ticks.y = element_line(color = "black"),
            axis.text.y = element_text(angle=0, vjust=0.3, size=10, color="black"),
            axis.title.y = element_text(face="bold", color="black", size=10),
            axis.title.x = element_text(face="bold", color="black", size=10),
            plot.title = element_text (face="bold", color="black", hjust = 0.5))+
  geom_vline(xintercept = 8.5, color = "black", linewidth=0.5)+
  geom_vline(xintercept = 16.5, color = "black", linewidth=0.5)+
  geom_vline(xintercept = 23.5, color = "black", linewidth=0.5)+
  geom_vline(xintercept = 32.5, color = "black", linewidth=0.5)+
  annotate(geom="text", x = 5, y = -3000, label="Free AMF vs H2O", fontface="bold", size=3)+
  annotate(geom="text", x = 12.5, y = -3000, label="Encased AMF vs H2O", fontface="bold", size=3)+
  annotate(geom="text", x = 20, y = -3000, label="vs WT (Free AMF)", fontface="bold", size=3)+
  annotate(geom="text", x = 27, y = -3000, label="vs WT (Encased AMF)", fontface="bold", size=3)+
  annotate(geom="text", x = 0.5, y = 2500, label="UP", fontface="bold", size=3, hjust=0, color="yellow4")+
  annotate(geom="text", x = 0.5, y = -3000, label="DOWN", fontface="bold", size=3, hjust=0, color="blue")+
  scale_fill_manual(values=c("yellow", "blue"), name="")+
  scale_x_discrete(name="", labels=c("WT","d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops",
                                     "WT","d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops",
                                     "d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops",
                                     "d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops"))+
  scale_y_continuous(name="Number of DEGs", breaks = c(-3000,-2000,-1000,0,1000,2000), limits=c(-3000,2500))+ 
  geom_hline(yintercept=0, color = "black", linewidth=0.4)

ggsave("plots/DEG_number/figure_S6.pdf", height=3.5, width=9)

# Normalised gene count plots --------------------------------------------------
# Figure S9B
genotype_colours_dark<-c("gray40","purple","dodgerblue3","blue3","orange4",
                         "darkgreen", "springgreen4", "chartreuse4")

genotype_colours_light<-c("grey","lavender","skyblue","dodgerblue","goldenrod1",
                          "palegreen", "springgreen1", "chartreuse1")


d<-filter(normalized_counts_for_plotting, RAP == "Os01g0657100", genotype %in% c("A","B","I","C","J","E","F","G","H"), treatment %in% c("01","12","13"))
d$genotype<-factor(d$genotype, levels=c("A","B","I","C","J","E","F","G","H"))
d<-mutate(d,value = ifelse(value == 0, 1, value))
ggplot(d, aes(x=genotype_name, y=value, color=genotype_name))+
  geom_point(mapping = aes(fill=genotype_name, colour=genotype_name), size = 1.2, shape=21)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="crossbar", width = 0.3)+
  publication_theme_counts()+
  scale_fill_manual(values=genotype_colours_light, name="Genotype")+
  scale_colour_manual(values=genotype_colours_dark, name="Genotype") +
  scale_x_discrete(name="Genotype")+
  scale_y_continuous(
                     trans=log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x),
                     labels = trans_format("log10", math_format(10^.x)),
                     name="Normalised counts")+
  facet_grid(~ treatment_name, scales = "free_x", labeller=as_labeller(c("H2O"="H2O","Free AMF"="Free AMF","flg22"="flg22")))+
  ggtitle("PT11")+
  stat_compare_means(label.y =3.1, label.x = 0.7, hjust=0, size=2.5, method="anova")+ 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "WT", label.y = 2.9, hide.ns = TRUE, size=3)+
  force_panelsizes(cols = unit(c(4.2,4.2,1), "cm"), TRUE)
ggsave(filename = "plots/gene_count_plots/figure_S9B_1.pdf", width = 4.5, height = 3.7)

d<-filter(normalized_counts_for_plotting, RAP == "Os01g0783000", genotype %in% c("A","B","I","C","J","E","F","G","H"), treatment %in% c("01","12","13"))
d$genotype<-factor(d$genotype, levels=c("A","B","I","C","J","E","F","G","H"))
d<-mutate(d,value = ifelse(value == 0, 1, value))
ggplot(d, aes(x=genotype_name, y=value, color=genotype_name))+
  geom_point(mapping = aes(fill=genotype_name, colour=genotype_name), size = 1.2, shape=21)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="crossbar", width = 0.3)+
  publication_theme_counts()+
  scale_fill_manual(values=genotype_colours_light, name="Genotype")+
  scale_colour_manual(values=genotype_colours_dark, name="Genotype") +
  scale_x_discrete(name="Genotype")+
  scale_y_continuous(
    trans=log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    name="Normalised counts")+
  facet_grid(~ treatment_name, scales = "free_x", labeller=as_labeller(c("H2O"="H2O","Free AMF"="Free AMF","flg22"="flg22")))+
  ggtitle("AM3")+
  stat_compare_means(label.y =3.1,label.x = 0.7, hjust=0, size=2.5, method="anova")+ 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "WT", label.y = 2.9, hide.ns = TRUE, size=3)+
  force_panelsizes(cols = unit(c(4.2,4.2,1), "cm"), TRUE)
ggsave(filename = "plots/gene_count_plots/figure_S9B_2.pdf", width = 4.5, height = 3.7)

d<-filter(normalized_counts_for_plotting, RAP == "Os04g0134800", genotype %in% c("A","B","I","C","J","E","F","G","H"), treatment %in% c("01","12","13"))
d$genotype<-factor(d$genotype, levels=c("A","B","I","C","J","E","F","G","H"))
d<-mutate(d,value = ifelse(value == 0, 1, value))
ggplot(d, aes(x=genotype_name, y=value, color=genotype_name))+
  geom_point(mapping = aes(fill=genotype_name, colour=genotype_name), size = 1.2, shape=21)+
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  stat_summary(fun=mean, geom="crossbar", width = 0.3)+
  publication_theme_counts()+
  scale_fill_manual(values=genotype_colours_light, name="Genotype")+
  scale_colour_manual(values=genotype_colours_dark, name="Genotype") +
  scale_x_discrete(name="Genotype")+
  scale_y_continuous(
    trans=log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)),
    name="Normalised counts")+
  facet_grid(~ treatment_name, scales = "free_x", labeller=as_labeller(c("H2O"="H2O","Free AMF"="Free AMF","flg22"="flg22")))+
  ggtitle("AM1")+
  stat_compare_means(label.y =2.5, label.x = 0.7, hjust=0, size=2.5, method="anova")+ 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "WT", label.y = 2.3, hide.ns = TRUE, size=3)+
  force_panelsizes(cols = unit(c(4.2,4.2,1), "cm"), TRUE)
ggsave(filename = "plots/gene_count_plots/figure_S9B_3.pdf", width = 4.5, height = 3.7)


# Scatterplots ----------------------------------------------------------------------
# Figure S5D
scatter_data_gse<-left_join(Gutjahr2015_gse_deg, filter(FC_for_plotting, comparison == "A07_vs_A01", RAP %in% Gutjahr2015_gse_deg$RAP), by="RAP")
ggplot(scatter_data_gse, aes(x=Log2FC_WT.M_vs_WT.G, y=FC))+
  geom_point(color = "black", size = 2, alpha=0.5) +
  geom_smooth(method = "lm", color = "dodgerblue", fill="dodgerblue", se = TRUE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, color="dodgerblue") +
  publication_theme_scatterplot() +
  labs(
    x = "Ri naïve GSE vs H2O (Gutjahr et al., 2015)",
    y = "Ri naïve GSE vs H2O (this study)"
  )
ggsave("plots/Scatterplots/figure_S5D_1.pdf", width = 3.5, height = 2.5)

scatter_data_gse<-left_join(Gutjahr2015_gse_deg, filter(FC_for_plotting, comparison == "A09_vs_A01", RAP %in% Gutjahr2015_gse_deg$RAP), by="RAP")
ggplot(scatter_data_gse, aes(x=Log2FC_WT.M_vs_WT.G, y=FC))+
  geom_point(color = "black", size = 2, alpha=0.5) +
  geom_smooth(method = "lm", color = "dodgerblue", fill="dodgerblue", se = TRUE) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
               formula = y ~ x, 
               parse = TRUE, color="dodgerblue") +
  publication_theme_scatterplot() +
  labs(
    x = "Ri naïve GSE vs H2O (Gutjahr et al., 2015)",
    y = "Mo naïve GSE vs H2O (this study)"
  )

ggsave("plots/Scatterplots/figure_S5D_2.pdf", width = 3.5, height = 2.5)



# Correlations -----------------------------------------------------------------
# Figure 2F
# Calculate correlation matrix
cor_matrix <- FC.table[,c("RAP","A06_vs_A01","A07_vs_A01","A08_vs_A01","A09_vs_A01","A10_vs_A01",
                          "A11_vs_A01","A12_vs_A01","A13_vs_A01")] %>%
  mutate(across(-RAP, ~ replace_na(., 0))) %>%
  dplyr::select(-RAP) %>% # Remove the "RAP" column
  cor(method = "pearson") # Compute Pearson correlation coefficients
colnames(cor_matrix)<-c("Plant exudates","Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Encased AMF","Free AMF","flg22")
rownames(cor_matrix)<-c("Plant exudates","Ri naïve GSE","Ri primed GSE","Mo naïve GSE","Mo primed GSE","Encased AMF","Free AMF","flg22")

# Convert correlation matrix to long format for heatmap plotting
cor_long <- cor_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Comparison1") %>%
  pivot_longer(cols = -Comparison1, names_to = "Comparison2", values_to = "Correlation")

# Plot the heatmap
pdf(file = "plots/Correlation/figure_2F.pdf", width = 4.2, height = 3.8)
Heatmap(
  cor_matrix,
  name = "Correlation",
  col = colorRamp2(c(0.3, 1), c("blue", "yellow")),
  cluster_columns = TRUE, row_km = 3, column_km =3,
  row_names_side = "left",
  cluster_rows = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE, show_column_dend = T, show_row_dend = F,
  column_names_rot = 90,
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  row_title_gp = gpar(fontsize = 0.1, fontface="bold"),
  column_title_gp = gpar(fontsize = 0.1, fontface="bold")
)
dev.off()

# Figure 4C
# Calculate correlation matrix
cor_matrix <- FC.table[,c("RAP","A11_vs_A01","A12_vs_A01","A13_vs_A01",
                          "B11_vs_B01",
                          "I11_vs_I01",
                          "C11_vs_C01",
                          "J11_vs_J01",
                          "F11_vs_F01",
                          "G11_vs_G01",
                          "H11_vs_H01"
                          )] %>%
  dplyr::mutate(across(-RAP, ~ replace_na(., 0))) %>%
  dplyr::select(-RAP) %>% # Remove the "RAP" column
  cor(method = "pearson") # Compute Pearson correlation coefficients
colnames(cor_matrix)<-c("WT","Free AMF, WT","flg22, WT","d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops")
rownames(cor_matrix)<-c("WT","Free AMF, WT","flg22, WT","d14l","cerk1","cerk1/2","nfr5","pollux","ccamk","cyclops")

# Plot the heatmap
pdf(file = "plots/Correlation/figure_4C.pdf", width = 4.4, height = 4)
Heatmap(
  cor_matrix,
  name = "Correlation",
  col = colorRamp2(c(0, 1), c("blue", "yellow")),
  cluster_columns = T,row_km = 2, column_km =2,
  cluster_rows = T,row_names_side = "left",
  show_row_names = TRUE,
  show_column_names = TRUE, show_column_dend = T, show_row_dend = F,
  column_names_rot = 90,
  column_names_side = "bottom",
  row_names_gp = gpar(fontsize = 10,fontface="italic"),
  column_names_gp = gpar(fontsize = 10,fontface="italic"),
  row_title_gp = gpar(fontsize = 0.1, fontface="bold"),
  column_title_gp = gpar(fontsize = 0.1, fontface="bold"),
  na_col = "white"
)
dev.off()


# GO enrichment ----------------------------------------------------------------
# GO enrichment analysis for each DEG set
rice_TERM2GENE<-read_tsv(file = "info_sheets/Osativa_GO_TERM2GENE.tsv", col_names = T, show_col_types = FALSE)
rice_TERM2NAME<-read_tsv(file = "info_sheets/Osativa_GO_TERM2NAME.tsv", col_names = T, show_col_types = FALSE)
rice_universe_GO<-filter(rice_TERM2GENE, RAP %in% rownames(normalized_counts)) %>% pull(RAP) %>% unique()

GO_enrichment_Osativa<-function(gene.list, name, folder_tables="plots/GO/GO_tables/")
{
  ego <- enricher(gene         = gene.list,
                  universe     = rice_universe_GO, 
                  TERM2GENE = rice_TERM2GENE, 
                  TERM2NAME = rice_TERM2NAME, 
                  pvalueCutoff = 0.05)
  write_tsv(file=paste(folder_tables,name,"_GO_table.tsv",sep=""), ego@result)
}

for(i in 1:length(DEG.list))
{
  try(GO_enrichment_Osativa(DEG.list[[i]], names(DEG.list)[i]))
}

GO_enrichment_Osativa(GSRgenes.rice$RAP, "GSRgenes.rice")

# Reading back the GO enrichment tables

GO.table.names<-substr(list.files("plots/GO/GO_tables"),1,nchar(list.files("plots/GO/GO_tables"))-4)
GO.table.list<-list()
for(i in 1:length(GO.table.names))
{
  GO.table<-read_tsv(paste("plots/GO/GO_tables/",GO.table.names[i],".tsv", sep=""), show_col_types = FALSE)
  GO.table.list[[i]]<-GO.table
}
names(GO.table.list)<-GO.table.names
GO.table.complete <- bind_rows(GO.table.list, .id = "comparison")

# Generating figures
# Figure 3B
GO.table<-filter(GO.table.complete, 
                 ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term,
                 comparison %in% c("activated_genes_table_A11_vs_A01_GO_table",
                                   "GSRgenes.rice_GO_table",
                                   "activated_genes_table_A13_vs_A01_GO_table",
                                   "activated_genes_table_A12_vs_A01_GO_table"),
                 pvalue <= 0.05)

GeneRatios<-c()
for(i in 1:length(GO.table$GeneRatio))
{
  GeneRatios[i]<-eval(parse(text = GO.table$GeneRatio[i]))
}
GO.table$GeneRatio<-GeneRatios

GO.table.selection<-rbind(
  filter(GO.table, comparison == "activated_genes_table_A12_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A11_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "GSRgenes.rice_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A13_vs_A01_GO_table")[1:5,]
)
GO.table.filtered<-filter(GO.table, ID %in% unique(GO.table.selection$ID))
idx<-order(GO.table.filtered$GeneRatio, decreasing = T)
GO.table<-filter(GO.table.filtered[idx,], ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term)
GO.table.filtered$Description<-factor(GO.table.filtered$Description, levels=unique(GO.table.selection$Description))
GO.table.filtered$comparison<-factor(GO.table.filtered$comparison, levels=c("activated_genes_table_A11_vs_A01_GO_table",
                                                                            "GSRgenes.rice_GO_table",
                                                                            "activated_genes_table_A13_vs_A01_GO_table",
                                                                            "activated_genes_table_A12_vs_A01_GO_table"))


ggplot(GO.table.filtered, aes(x=GeneRatio, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  geom_point(shape = 1, colour = "black") +
  scale_color_gradientn(colours=c("yellow","blue","midnightblue"),values=c(0,0.5,1), name = "FDR") +
  ylab("GO terms") + 
  xlab("Gene Ratio") + 
  publication_theme_GO() +
  scale_size(name="Count", range=c(1,10), breaks=c(10,50,100,250)) +
  facet_wrap(vars(comparison), nrow = 1, ncol=8, labeller = as_labeller(c("activated_genes_table_A11_vs_A01_GO_table"="Encased AMF",
                                                                          "GSRgenes.rice_GO_table"="GSR genes",
                                                                          "activated_genes_table_A12_vs_A01_GO_table"="Free AMF",
                                                                          "activated_genes_table_A13_vs_A01_GO_table"="flg22"
  )))
ggsave("plots/GO/GO_plots/figure_3B.pdf", height=3.7, width=8.5)

# Figure S2A
GO.table<-filter(GO.table.complete, 
                 ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term,
                 comparison %in% c("activated_genes_table_A12_vs_A01_GO_table",
                                   "repressed_genes_table_A12_vs_A01_GO_table"),
                 pvalue <= 0.05)

GeneRatios<-c()
for(i in 1:length(GO.table$GeneRatio))
{
  GeneRatios[i]<-eval(parse(text = GO.table$GeneRatio[i]))
}
GO.table$GeneRatio<-GeneRatios
idx<-order(GO.table$GeneRatio, decreasing = T)
GO.table<-filter(GO.table[idx,], ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term)
GO.table$Description<-factor(GO.table$Description, levels=unique(GO.table$Description))

GO.table.selection<-rbind(
  filter(GO.table, comparison == "activated_genes_table_A12_vs_A01_GO_table")[1:20,],
  filter(GO.table, comparison == "repressed_genes_table_A12_vs_A01_GO_table")[1:20,])

ggplot(GO.table.selection, aes(x=GeneRatio, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  geom_point(shape = 1, colour = "black") +
  scale_color_gradientn(colours=c("yellow","blue","midnightblue"),values=c(0,0.5,1), name = "FDR") +
  ylab("GO terms") + 
  xlab("Gene Ratio") + 
  publication_theme_GO() +
  scale_size(name="Count") +
  facet_wrap(vars(comparison), nrow = 2, scales = "free", labeller = as_labeller(c("activated_genes_table_A12_vs_A01_GO_table"="Free AMF  up-regulated",
                                                                                   "repressed_genes_table_A12_vs_A01_GO_table"="Free AMF  down-regulated")))
ggsave("plots/GO/GO_plots/figure_S2C.pdf", height=7, width=6)

# Figure S3C
GO.table<-filter(GO.table.complete, 
                 ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term,
                 comparison %in% c("activated_genes_table_A13_vs_A01_GO_table",
                                   "repressed_genes_table_A13_vs_A01_GO_table"),
                 pvalue <= 0.05)

GeneRatios<-c()
for(i in 1:length(GO.table$GeneRatio))
{
  GeneRatios[i]<-eval(parse(text = GO.table$GeneRatio[i]))
}
GO.table$GeneRatio<-GeneRatios
idx<-order(GO.table$GeneRatio, decreasing = T)
GO.table<-filter(GO.table[idx,], ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term)
GO.table$Description<-factor(GO.table$Description, levels=unique(GO.table$Description))

GO.table.selection<-rbind(
  filter(GO.table, comparison == "activated_genes_table_A13_vs_A01_GO_table")[1:20,],
  filter(GO.table, comparison == "repressed_genes_table_A13_vs_A01_GO_table")[1:20,])

ggplot(GO.table.selection, aes(x=GeneRatio, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  geom_point(shape = 1, colour = "black") +
  scale_color_gradientn(colours=c("yellow","blue","midnightblue"),values=c(0,0.5,1), name = "FDR") +
  ylab("GO terms") + 
  xlab("Gene Ratio") + 
  publication_theme_GO() +
  scale_size(name="Count", range=c(1,10), breaks=c(50,100,150,200,250)) +
  facet_wrap(vars(comparison), nrow = 2, scales = "free", labeller = as_labeller(c("activated_genes_table_A13_vs_A01_GO_table"="flg22  up-regulated",
                                                                                   "repressed_genes_table_A13_vs_A01_GO_table"="flg22  down-regulated")))
ggsave("plots/GO/GO_plots/figure_S3C.pdf", height=7, width=6)

# Figure S4
GO.table<-filter(GO.table.complete, 
                 ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term,
                 comparison %in% c(
                                   "activated_genes_table_A07_vs_A01_GO_table",
                                   "activated_genes_table_A08_vs_A01_GO_table",
                                   "activated_genes_table_A09_vs_A01_GO_table",
                                   "activated_genes_table_A10_vs_A01_GO_table",
                                   "activated_genes_table_A06_vs_A01_GO_table",
                                   "activated_genes_table_A11_vs_A01_GO_table",
                                   "activated_genes_table_A13_vs_A01_GO_table",
                                   "activated_genes_table_A12_vs_A01_GO_table"
                                   ),
                 pvalue <= 0.05)

GeneRatios<-c()
for(i in 1:length(GO.table$GeneRatio))
{
  GeneRatios[i]<-eval(parse(text = GO.table$GeneRatio[i]))
}
GO.table$GeneRatio<-GeneRatios

GO.table.selection<-rbind(
  
  filter(GO.table, comparison == "activated_genes_table_A07_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A08_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A09_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A10_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A06_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A11_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A13_vs_A01_GO_table")[1:5,],
  filter(GO.table, comparison == "activated_genes_table_A12_vs_A01_GO_table")[1:5,]
)
GO.table.filtered<-filter(GO.table, ID %in% unique(GO.table.selection$ID))
idx<-order(GO.table.filtered$GeneRatio, decreasing = T)
GO.table<-filter(GO.table.filtered[idx,], ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term)
GO.table.filtered$Description<-factor(GO.table.filtered$Description, levels=unique(GO.table.selection$Description))
GO.table.filtered$comparison<-factor(GO.table.filtered$comparison, levels=c(
  "activated_genes_table_A07_vs_A01_GO_table",
  "activated_genes_table_A08_vs_A01_GO_table",
  "activated_genes_table_A09_vs_A01_GO_table",
  "activated_genes_table_A10_vs_A01_GO_table",
  "activated_genes_table_A06_vs_A01_GO_table",
  "activated_genes_table_A11_vs_A01_GO_table",
  "activated_genes_table_A13_vs_A01_GO_table",
  "activated_genes_table_A12_vs_A01_GO_table"
))


ggplot(GO.table.filtered, aes(x=GeneRatio, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  geom_point(shape = 1, colour = "black") +
  scale_color_gradientn(colours=c("yellow","blue","midnightblue"),values=c(0,0.5,1), name = "FDR") +
  ylab("GO terms") + 
  xlab("Gene Ratio") + 
  publication_theme_GO() +
  scale_size(name="Count", range=c(1,15), breaks=c(10,50,100,250)) +
  facet_wrap(vars(comparison), nrow = 2, ncol=4, labeller = as_labeller(c("activated_genes_table_A07_vs_A01_GO_table"="Ri naïve GSE",
                                                                          "activated_genes_table_A08_vs_A01_GO_table"="Ri primed GSE",
                                                                          "activated_genes_table_A09_vs_A01_GO_table"="Mo naïve GSE",
                                                                          "activated_genes_table_A10_vs_A01_GO_table"="Mo primed GSE",
                                                                          "activated_genes_table_A06_vs_A01_GO_table"="Plant exudates",
                                                                          "activated_genes_table_A11_vs_A01_GO_table"="Encased AMF",
                                                                          "activated_genes_table_A12_vs_A01_GO_table"="Free AMF",
                                                                          "activated_genes_table_A13_vs_A01_GO_table"="flg22")))
ggsave("plots/GO/GO_plots/figure_S4.pdf", height=7, width=10)

# Figure S8B

GO.table<-filter(GO.table.complete, 
                 ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term,
                 comparison %in% c("activated_genes_table_A11_vs_A01_GO_table",
                                   "activated_genes_table_B11_vs_B01_GO_table",
                                   "activated_genes_table_I11_vs_I01_GO_table",
                                   "activated_genes_table_C11_vs_C01_GO_table",
                                   "activated_genes_table_J11_vs_J01_GO_table",
                                   "activated_genes_table_F11_vs_F01_GO_table",
                                   "activated_genes_table_G11_vs_G01_GO_table",
                                   "activated_genes_table_H11_vs_H01_GO_table",
                                   "GSRgenes.rice_GO_table"),
                 pvalue <= 0.05)

GeneRatios<-c()
for(i in 1:length(GO.table$GeneRatio))
{
  GeneRatios[i]<-eval(parse(text = GO.table$GeneRatio[i]))
}
GO.table$GeneRatio<-GeneRatios

GO.table.selection<-rbind(
  filter(GO.table, comparison == "activated_genes_table_A11_vs_A01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_B11_vs_B01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_I11_vs_I01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_C11_vs_C01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_J11_vs_J01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_F11_vs_F01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_G11_vs_G01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_H11_vs_H01_GO_table")[1:3,],
  filter(GO.table, comparison == "GSRgenes.rice_GO_table")[1:3,])
GO.table.filtered<-filter(GO.table, ID %in% unique(GO.table.selection$ID))
idx<-order(GO.table.filtered$GeneRatio, decreasing = T)
GO.table.filtered<-filter(GO.table.filtered[idx,], ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term)
GO.table.filtered$Description<-factor(GO.table.filtered$Description, levels=unique(GO.table.selection$Description))
GO.table.filtered$comparison <- factor(GO.table.filtered$comparison,levels=c("activated_genes_table_A11_vs_A01_GO_table",
                                                                             "activated_genes_table_B11_vs_B01_GO_table",
                                                                             "activated_genes_table_I11_vs_I01_GO_table",
                                                                             "activated_genes_table_C11_vs_C01_GO_table",
                                                                             "activated_genes_table_J11_vs_J01_GO_table",
                                                                             "activated_genes_table_F11_vs_F01_GO_table",
                                                                             "activated_genes_table_G11_vs_G01_GO_table",
                                                                             "activated_genes_table_H11_vs_H01_GO_table",
                                                                             "GSRgenes.rice_GO_table"))


ggplot(GO.table.filtered, aes(x=GeneRatio, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  geom_point(shape = 1, colour = "black") +
  scale_color_gradientn(colours=c("yellow","blue","midnightblue"),values=c(0,0.5,1), name = "FDR") +
  ylab("GO terms") + 
  xlab("Gene Ratio") + 
  publication_theme_GO() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=9, color="black"),
        strip.text = element_text(size=10, color="white", face="italic"))+
  theme(axis.text.x = element_text(angle=90, vjust=0.3, hjust=1, size=9, color="black"))+
  scale_size(name="Count") +
  facet_wrap(vars(comparison), nrow = 1, ncol=9, labeller = as_labeller(c("activated_genes_table_A11_vs_A01_GO_table"="WT",
                                                                          "activated_genes_table_B11_vs_B01_GO_table"="d14l",
                                                                          "activated_genes_table_I11_vs_I01_GO_table"="cerk1",
                                                                          "activated_genes_table_C11_vs_C01_GO_table"="cerk1/2",
                                                                          "activated_genes_table_J11_vs_J01_GO_table"="nfr5",
                                                                          "activated_genes_table_F11_vs_F01_GO_table"="pollux",
                                                                          "activated_genes_table_G11_vs_G01_GO_table"="ccamk",
                                                                          "activated_genes_table_H11_vs_H01_GO_table"="cyclops",
                                                                          "GSRgenes.rice_GO_table"="GSR")))
ggsave("plots/GO/GO_plots/figure_S8B.pdf", height=3.5, width=11)

# Figure S8C

GO.table<-filter(GO.table.complete, 
                 ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term,
                 comparison %in% c("activated_genes_table_A12_vs_A01_GO_table",
                                   "activated_genes_table_B12_vs_B01_GO_table",
                                   "activated_genes_table_I12_vs_I01_GO_table",
                                   "activated_genes_table_C12_vs_C01_GO_table",
                                   "activated_genes_table_J12_vs_J01_GO_table",
                                   "activated_genes_table_F12_vs_F01_GO_table",
                                   "activated_genes_table_G12_vs_G01_GO_table",
                                   "activated_genes_table_H12_vs_H01_GO_table"),
                 pvalue <= 0.05)

GeneRatios<-c()
for(i in 1:length(GO.table$GeneRatio))
{
  GeneRatios[i]<-eval(parse(text = GO.table$GeneRatio[i]))
}
GO.table$GeneRatio<-GeneRatios

GO.table.selection<-rbind(
  filter(GO.table, comparison == "activated_genes_table_A12_vs_A01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_B12_vs_B01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_I12_vs_I01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_C12_vs_C01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_J12_vs_J01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_F12_vs_F01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_G12_vs_G01_GO_table")[1:3,],
  filter(GO.table, comparison == "activated_genes_table_H12_vs_H01_GO_table")[1:3,])
GO.table.filtered<-filter(GO.table, ID %in% unique(GO.table.selection$ID))
idx<-order(GO.table.filtered$GeneRatio, decreasing = T)
GO.table.filtered<-filter(GO.table.filtered[idx,], ID %in% filter(rice_TERM2NAME, Ontology == "BP")$GO.term)
GO.table.filtered$Description<-factor(GO.table.filtered$Description, levels=unique(GO.table.selection$Description))
GO.table.filtered$comparison <- factor(GO.table.filtered$comparison,levels=c("activated_genes_table_A12_vs_A01_GO_table",
                                                           "activated_genes_table_B12_vs_B01_GO_table",
                                                           "activated_genes_table_I12_vs_I01_GO_table",
                                                           "activated_genes_table_C12_vs_C01_GO_table",
                                                           "activated_genes_table_J12_vs_J01_GO_table",
                                                           "activated_genes_table_F12_vs_F01_GO_table",
                                                           "activated_genes_table_G12_vs_G01_GO_table",
                                                           "activated_genes_table_H12_vs_H01_GO_table"))


ggplot(GO.table.filtered, aes(x=GeneRatio, y=Description, size=Count, color=p.adjust)) +
  geom_point() +
  geom_point(shape = 1, colour = "black") +
  scale_color_gradientn(colours=c("yellow","blue","midnightblue"),values=c(0,0.5,1), name = "FDR") +
  ylab("GO terms") + 
  xlab("Gene Ratio") + 
  publication_theme_GO() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1, size=9, color="black"),
        strip.text = element_text(size=10, color="white", face="italic"))+
  theme(axis.text.x = element_text(angle=90, vjust=0.3, hjust=1, size=9, color="black"))+
  scale_size(name="Count") +
  facet_wrap(vars(comparison), nrow = 1, ncol=9, labeller = as_labeller(c("activated_genes_table_A12_vs_A01_GO_table"="WT",
                                                                          "activated_genes_table_B12_vs_B01_GO_table"="d14l",
                                                                          "activated_genes_table_I12_vs_I01_GO_table"="cerk1",
                                                                          "activated_genes_table_C12_vs_C01_GO_table"="cerk1/2",
                                                                          "activated_genes_table_J12_vs_J01_GO_table"="nfr5",
                                                                          "activated_genes_table_F12_vs_F01_GO_table"="pollux",
                                                                          "activated_genes_table_G12_vs_G01_GO_table"="ccamk",
                                                                          "activated_genes_table_H12_vs_H01_GO_table"="cyclops")))
ggsave("plots/GO/GO_plots/figure_S8C.pdf", height=4, width=11)
