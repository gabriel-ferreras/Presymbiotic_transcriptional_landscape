# Presymbiotic transcriptional landscape
Source code used to analyse the data and generate graphics for the publication:

### Defining the pre-symbiotic landscape of rice roots

Gabriel Ferreras Garrucho<sup>1</sup>, Uta Paszkowski<sup>1*</sup><sup>§</sup>, Chai Hao Chiu<sup>1*</sup><sup>§</sup>

<sup>1</sup> Crop Science Centre, Department of Plant Sciences, University of Cambridge, 93 Lawrence Weaver Rd, Cambridge, United Kingdom, CB3 0LE.  
<sup>*</sup> up220@cam.ac.uk, chc59@cam.ac.uk  
<sup>§</sup> shared last author

---

### Repository Contents

- `Presym_DGE_analysis.R`: Script for processing transcriptomic data and conducting differential gene expression analysis.
- `Presym_Plotting_results.R`: Script for generating publication-quality figures from processed DGE results.
- `gene_count_tables/`: Contains pre-processed gene count matrices including `gene_count_matrix.csv`.
- `info_sheets/`: Contains tables with gene annotations and related metadata.

---

### Usage

1. **Clone or download this repository locally**

   ```bash
   git clone https://github.com/your-username/presymbiotic-landscape.git
   cd presymbiotic-landscape

2. **Obtain the raw sequencing data**

   Raw reads are available at NCBI GEO under accession GSE28377:
   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28377
   
   To reproduce the results:
   
   1. Download the raw reads (fastq files) from GEO.
   2. Index the Oryza sativa Nipponbare reference transcriptome (IRGSP-1.0.52) using Kallisto v0.46.1.
   3. Run pseudo-alignment using the following Kallisto command:
   
            kallisto quant -i rice_index.idx -o samples/sample1 sample1_R1.fastq sample1_R2.fastq
   
   Repeat for all samples. The output folders must be placed in a directory called 'samples/' inside the main repo folder.

3. **Run the analysis pipeline**

      A. Run `Presym_DGE_analysis.R` to process data and perform DGE analysis.
      
      Sections include:
      
      1. LOADING PACKAGES  
         - Install and load required R packages.  
         - Also loads general settings for plotting themes.
      
      2. TXIMPORT  
         - Imports transcript-level abundance estimates from Kallisto results in 'samples/'.  
         - Generates a gene-level count matrix using the tximport package.  
         - Alternatively, loads the precomputed count matrix from 'gene_count_tables/gene_count_matrix.csv'.
      
      3. DESEQ2  
         - Constructs DESeqDataSet from the count matrix.  
         - Performs sample QC, pre-filtering, factor design, and differential expression analysis.  
         - Applies variance-stabilising transformation (VST).  
         - Saves key output matrices for downstream plotting.
      
      4. DEG RESULTS  
         - Extracts statistically significant DEGs.  
         - Outputs results tables to the working directory.


      B. Run `Presym_Plotting_results.R` to generate figures for the publication.
      
      Sections include:
      
      1. LOADING PACKAGES  
         - Loads libraries and custom ggplot themes.
      
      2. LOADING DATA  
         - Loads gene annotations, gene lists, DEGs, and count matrices from local files.  
         - Includes FDR and log2 fold change matrices as well as VST-normalized counts.
      
      3. PCA  
         - Generates PCA plots from VST-transformed counts.  
         - Includes figures corresponding to published results (e.g., Figure 2G, 4B).

      4. HEATMAPS  
         - Plots gene expression heatmaps for specific gene groups.  
         - Includes main and supplemental figure panels.
      
      5. AVERAGE EXPRESSION PLOTS  
         - Plots average expression of gene sets or individual markers.  
         - Includes pairwise t-tests for significance and publication-quality bar/line plots.
      
      6. VENN  
         - Generates Venn diagrams for DEG overlap between conditions.
      
      7. DEG NUMBER PLOTS  
         - Plots number of DEGs per comparison for figure panels.
      
      8. NORMALISED GENE COUNT PLOTS  
         - Individualised plots for selected genes from the DEG sets.
      
      9. SCATTERPLOTS and CORRELATIONS  
         - Generates scatterplots and heatmaps showing fold-change and correlation matrices.
      
      10. GO ENRICHMENT  
         - Runs GO enrichment analysis for each DEG set and produces summary plots.

---

### Citation

If you use this code or dataset, please cite the associated publication (DOI and citation to be added upon acceptance).


