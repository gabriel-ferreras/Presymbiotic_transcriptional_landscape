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

    Raw sequences are available in the NCBI GEO database under accession GSE28377:
    https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28377

    To replicate the results from this study, pseudo-align these reads to the
    Oryza sativa Nipponbare reference transcriptome (IRGSP-1.0.52) using Kallisto v0.46.1
    with default parameters.

    Example command for one sample:

       kallisto quant -i rice_index.idx -o samples/sample1 sample1_R1.fastq sample1_R2.fastq

    Repeat for all samples. All output folders should be saved inside a directory
    named 'samples' within the main folder of this repository.

4. **Run the analysis pipeline**

        Presym_DGE_analysis.R:
        
          1. LOADING PACKAGES
             Install and load required packages.
        
          2. TXIMPORT
             Import Kallisto abundance estimates using the tximport package.
             Alternatively, use the pre-supplied gene count matrix in:
             gene_count_tables/gene_count_matrix.csv
        
          3. DESEQ2
             Normalize count data and perform differential gene expression analysis.
        
          4. CONTRASTS
             Define contrasts between sample groups to extract DEGs.
        
          5. OUTPUT RESULTS
             Export DEG tables to the working directory.


        Presym_Plotting_results.R:
        
          1. LOADING PACKAGES
             Load required libraries and general plot themes.
        
          2. PCA PLOT
             Generate a principal component analysis plot of the samples.
        
          3. SAMPLE DISTANCE HEATMAP
             Create a heatmap of Euclidean distances between all samples.
        
          4. GENE EXPRESSION HEATMAPS
             Use curated gene lists from info_sheets/ to generate expression heatmaps.
        
          5. SELECTED GENE PLOTS
             Create publication-quality bar and line plots for key marker genes.

---

### Citation

If you use this code or dataset, please cite the associated publication (DOI and citation to be added upon acceptance).


