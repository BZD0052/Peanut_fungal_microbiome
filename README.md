
# Peanut Fungal ITS Amplicon Sequencing

[![DOI](https://zenodo.org/badge/924459211.svg)](https://doi.org/10.5281/zenodo.15278087)

Author: Bibek Dabargainya

Date: Feb. 9, 2024

Please read and understand all the following steps
## Bioinformatic pipeline for peanut ITS amplicon sequencing 

Welcome to the fungal ITS amplicon sequencing pipeline — a comprehensive workflow designed to explore how fungal
communities shift, survive, and interact with peanut plants across environmental stress conditions. This pipeline 
was developed to process high-throughput ITS sequence data from plant-associated environments using the ASAX High
Performance Computing (HPC) at the Alabama Supercomputer Authority. This specific project focuses on understanding
the fungal community diversity in peanut cropping systems under contrasting moisture regimes. The analysis uses 
ITS-1 targeted single-end Illumina MiSeq reads and was applied to 384 demultiplexed samples collected across three
growing seasons (2022, 2023, and 2024). These samples were taken from three biologically significant peanut tissue
types (root, peg, and soil); representing both drought-tolerant and a drought-sensitive peanut variety. Field experiments
was conducted at E.V. Smith Plant Breeding Unit, Auburn University, Wiregrass, Alabama.The core aim of this project 
is to identify and quantify shifts in fungal community composition in response to water stress and plant genotype.
Whether applied to the rhizosphere, peg zone, or bulk soil, this pipeline enables structured, reproducible, and scalable
analysis of amplicon sequencing data; from raw .fastq.gz reads all the way to curated OTU tables and taxonomic annotations,
ready for ecological inference and statistical modeling.

This particular analysis centers on fungal ITS1 sequences derived from single-end forward Illumina reads. For consistency
and comparability, the analysis pipeline exclusively targets the forward read (R1) using the ITS1F primer CTTGGTCATTTAGAGGAAGTAA.


# Flourishing_fungi ASAX HPC Bioinfomatic Pipeline

This pipeline was adapted from the Flourishing Fungi ITS workflow originally developed by Dr. Zachary Noel for fungal amplicon
sequencing on ASAX. This was published on February 9, 2024 and was designed to streamline high-throughput fungal community 
analysis. Our modifications extend and specialize this pipeline for mid season drought studies in peanut agro-ecosystems.


[**Flourishing Fungi** pipeline by Zachary Noel (2024); (https://zenodo.org/doi/10.5281/zenodo.10655178)**](https://zenodo.org/doi/10.5281/zenodo.10655178)

Bellow is the summary for this pipeline:

1. Sample ID Extraction
   The first step generates a 'samples.txt' file that lists unique sample identifiers.
   This is achieved by stripping Illumina-specific suffixes (e.g., _R1_001.fastq.gz) from
   file names (ITS reads). This standardized list of base sample names is critical for ensuring
   consistency across all processing loops.

2. Primer Removal
   Each sample is then processed with cutadapt to remove the forward primer
   (CTTGGTCATTTAGAGGAAGTAA) specific to the ITS-1 region. Only forward reads (R1) are used
   in this pipeline. Sequences that do not contain the primer are discarded. Trimmed reads
   are then stored in the 'trimmed/' directory.

3. Read Quality Statistics
   After primer removal, all trimmed FASTQ files are concatenated into a single file, and
   vsearch is used to generate quality statistics. These include total read count, average
   read length, and quality score distributions. The resulting file ('stats/stats_results.txt')
   informs the filtering thresholds applied in the next step.

4. Filtering and Truncation
   The pipeline uses vsearch to apply the following filtering criteria:
   - Maximum expected error: 1
   - No ambiguous bases (Ns): 0
   - Truncate reads to 264 bp
   - Optionally trim the first 44 bp from each read
   Filtered reads are then saved as both FASTA and FASTQ files in the 'filtered/' directory.

5. Dereplication and OTU Clustering
   Filtered reads are dereplicated using vsearch. Then usearch clusters sequences into
   Operational Taxonomic Units (OTUs) at 97% similarity, ignoring singletons (minsize = 2).
   All OTUs are labeled with a 'FOTU_' prefix and saved in the 'clustered/' directory.

6. Mapping Reads to OTUs
   Each sample's filtered FASTQ reads are converted to FASTA format using fastq_to_fasta.
   A Python script ('replace_fastaheaders2.py') modifies sample IDs to each header. Then all
   FASTA files are concatenated and mapped to clustered OTUs using vsearch. This generates
   the OTU abundance table:
     otu_table/otu_table_ITS_Fungi.txt

7. Taxonomy Assignment using CONSTAX2. 
   This step uses CONSTAX2 to assign taxonomy to each OTU based on consensus predictions from
   BLAST, SINTAX, and RDP classifiers, using the UNITE v10 dynamic ITS database with a confidence
   threshold of 0.8. The output is saved in the taxonomy_assignments/ directory and includes
   tabular CSV files suitable for downstream use in R. 
   
   > **Note:** This step is not yet publicly available in the original Flourishing Fungi pipeline
   by Dr. Zachary Noel, but is included here with permission. It is based on internal, unpublished
   code developed by Dr. Noel for CONSTAX2 integration on ASAX. This functionality may be included
   in future public releases of the original pipeline.
  
# R Analysis: Peanut Fungi

This R pipeline processes fungal ITS1 amplicon sequencing data from peanut samples
collected under drought and irrigated conditions across three years. It includes
steps for creating phyloseq objects, decontamination, diversity analysis, normalization,
PERMANOVA, and DESeq2-based differential abundance testing. This R analysis begins by
creating a phyloseq object. Before beginning, ensure that you have the following files
downloaded and in an appropriate directory so that R can utilize them:

---[OTU_TABLE_CSV_FUNGI.csv](Phyloseq_input/OTU_TABLE_CSV_FUNGI.csv)
---[In_use_constax_taxonomy.csv](Phyloseq_input/In_use_constax_taxonomy.csv)
---[FungalMetadata_2024.csv](Phyloseq_input/FungalMetadata_2024.csv)
---[otus.fasta](Phyloseq_input/otus.fasta)


  ## Load Dependencies
  Load required R libraries necessary for microbiome and statistical analysis. This pipeline includes :
  
  | Package        | Version  | Description                                                  |
  |----------------|----------|--------------------------------------------------------------|
  | phyloseq       | 1.48.0   | Microbiome data import, transformation, and analysis         |
  | DESeq2         | 1.44.0   | Differential abundance testing using negative binomial models|
  | vegan          | 2.6-8    | Community ecology analysis and PERMANOVA                     |
  | metagenomeSeq  | 1.46.0   | CSS normalization and abundance modeling for sparse count data|
  | decontam       | 1.24.0   | Statistical identification of contaminants                   |
  | Biostrings     | 2.72.1   | Efficient manipulation of biological sequences               |
  | ggplot2        | 3.5.2    | Core plotting package for data visualization                 |
  | ggpubr         | 0.6.0    | Publication-ready visualizations built on ggplot2            |
  | ggrepel        | 0.9.6    | Smart label placement in plots                               |
  | microbiome     | 1.26.0   | High-level tools for microbiome data exploration             |
  | tidyverse      | 2.0.0    | Collection of packages for tidy data science                 |
  | dplyr          | 1.1.4    | Data manipulation grammar (part of tidyverse)                |
  
  

  ## Create Phyloseq Object
   Load OTU table, taxonomy file, metadata, and fasta sequences from the pipeline output.
   Integrate all data into a single phyloseq object. Sample names are harmonized across
   OTU and metadata to ensure consistency.

  ## Decontamination
   Use 'decontam' to identify and remove potential contaminants based on prevalence in
   negative controls. Contaminants are visualized and removed before further analysis.

  ## Filtering and Cleaning
   In this step, only the samples labeled as 'True Sample' were retained, which excludes
   controls and blank samples to focus on the actual biological data. We specifically filtered
   for taxa classified under the Kingdom 'Fungi' to narrow the analysis to relevant fungal taxa.
   Any taxa with zero counts across all samples were removed to ensure that only those with 
   measurable presence were considered in subsequent analyses. Additionally, samples with fewer
   than 5,000 total reads were excluded, as they were deemed too low in depth for reliable analysis.
   These filtering steps ensure that the data used in the analysis is both relevant and of
   adequate quality.

  ## Summary Statistics
   After filtering the dataset, we calculated key statistics to assess sequencing depth across
   all retained samples. This included the total number of reads remaining, along with the mean
   and median read depths per sample. To visualize the distribution of sequencing effort,
   a histogram was generated showing the number of reads per sample, with a dashed vertical line
   indicating the median read depth. This step provides a quick overview of data quality before
   diversity and abundance analyses

  ## Rarefaction
   Rarefaction curves were generated to evaluate sequencing depth sufficiency and OTU richness
   across samples. Median read depth was overlaid to visualize typical sampling effort, and
   final plots were grouped by treatment and tissue to highlight differences in richness
   accumulation under each condition.

  ## Alpha Diversity
   Diversity indices including Shannon, Inverse Simpson, Observed Richness, and Evenness were
   calculated to assess within-sample fungal diversity. A 3-way ANOVA was performed separately
   by year to evaluate the effects of treatment, tissue type, and plant variety. Results were 
   visualized using faceted boxplots and saved for interpretation.

  ## CSS Normalization (metagenomeSeq)
   Read counts were normalized using cumulative sum scaling (CSS) via the metagenomeSeq package
   to correct for library size differences across samples. The resulting normalized phyloseq 
   object was saved for downstream beta diversity and differential abundance analyses.-
   
  ## Beta Diversity & PERMANOVA
   Bray-Curtis distance matrices were used to visualize community dissimilarities through PCoA.
   Global PERMANOVA tested the effects of Treatment, Tissue, and plant variety interactions, 
   while tissue and year specific PERMANOVA, beta-dispersion, and ANOSIM analyses provided
   fined tuned resolution. All outputs, including ordination plots and statistical results,
   were saved by year and tissue.

  ## Differential Abundance (DESeq2)
   For each year (2022–2024), DESeq2 was applied to Peg tissue samples under Drought vs Irrigated
   conditions. Significant OTUs were filtered, annotated with taxonomy and visualized using volcano
   plots colored by fungal class. Final figures were saved for each year to highlight key shifts 
   in community structure. 

# Output Summary

All plots are saved in the [`Plots/`](Plots/) directory. Statistical outputs (e.g., ANOVA tables, PERMANOVA results) are stored in the [`Tables/`](Tables/) directory, organized by analysis type.
Following is th quick view table for the reults.

| **Analysis Type**         | **Description**                                                 | **Results Table Link**                                      | **Plot Link**                                      |
|--------------------------|------------------------------------------------------------------|--------------------------------------------------------------|----------------------------------------------------|
| **Alpha Diversity**       | ANOVA tables by year for Shannon, Simpson, Richness, Evenness   | [`Tables/Alpha_diversity/`](Tables/Alpha_diversity/)         | [`Plots/Alpha_diversity/`](Plots/Alpha_diversity/) |
| **Global PERMANOVA**      | Bray-Curtis PERMANOVA on Treatment × Tissue × Variety           | [`Tables/Global_permanova/`](Tables/Global_permanova/)       | [`Plots/Beta_diversity_PCoA/`](Plots/Beta_diversity_PCoA/) |
| **PERMANOVA (by Tissue)** | Tissue-level PERMANOVA tests by year                           | [`Tables/Permanova/`](Tables/Permanova/)                     | [`Plots/Beta_diversity_PCoA/`](Plots/Beta_diversity_PCoA/) |
| **Beta Dispersion**       | Betadispersion test for homogeneity of dispersions              | [`Tables/Beta_dispersion/`](Tables/Beta_dispersion/)         | [`Plots/Beta_diversity_PCoA/`](Plots/Beta_diversity_PCoA/) |
| **ANOSIM**                | Analysis of similarity for treatment separation                 | [`Tables/Anosim/`](Tables/Anosim/)                           | [`Plots/Beta_diversity_PCoA/`](Plots/Beta_diversity_PCoA/) |
| **Differential Abundance**| DESeq2 differential analysis for Peg tissue by year             | *Not available (plot only)*                                 | [`Plots/Differential_abundance/`](Plots/Differential_abundance/) |

> **Note:** PCoA plots shown in the "Plot Link" column are available only for Peg tissue subsets by year, as Peg was identified
as the most significantly responsive tissue across years explaning the variation in fungal community compistion between our treatment.


# Citation and Contact

If this pipeline contributes to your research, please consider citing this repository in your publications or presentations.
For questions, bug reports, or collaboration inquiries, feel free to reach out:

**Bibek Dabargainya**  
Graduate Researcher, Auburn University  
bzd0052@auburn.edu

# File Tree

```
├──📁 HPC_script
│   ├──📄 Flourishing_fingi_pipeline.sh         # Bash pipeline script for ITS amplicon reads
│   ├──📄 replace_fastaheaders2.py              # Python script for renaming FASTA headers 
│   └──📄 samples.txt                           # Sample ID list used in loops
├──📄 Peanut_fungal_microbiome.Rproj            # RStudio project file
├──🌐 Peanut_Fungi_Final.html                   # Final HTML output
├──📝 Peanut_Fungi_Final.md                     # GitHub Markdown version of analysis
├──📘 Peanut_Fungi_Final.Rmd                    # Main RMarkdown file
├──📁 Phyloseq_input
│   ├──📊 FungalMetadata_2024.csv               # Metadata table (CSV format)
│   ├──💾 Fungi-phyloseq-clean-CSS.rds          # RDS object after CSS normalization
│   ├──💾 Fungi-phyloseq-clean.rds              # Filtered phyloseq object
│   ├──💾 Fungi-phyloseq.rds                    # Initial raw phyloseq object
│   ├──📊 In_use_constax_taxonomy.csv           # OTU taxonomy assignments
│   ├──🧬 otus.fasta                            # FASTA file of representative sequences
│   └──📊 OTU_TABLE_CSV_FUNGI.csv               # OTU abundance table
├──📁 Plots
│   ├──📁 Alpha_diversity                         # Boxplots and combined figures for alpha diversity
│   │   ├──🖼️ Eveness_final.png
│   │   ├──🖼 Final_combined_alpha_diversity_plot1.png
│   │   ├──🖼 Invsimpson_final.png
│   │   ├──🖼 Richness_final.png
│   │   └──🖼 shannon_final.png
│   ├──📁 Beta_diversity_PCoA                     # Ordination plots by year and tissue
│   │   ├──🖼 PCoA_Peg_AllYears.png
│   │   ├──🖼 PCoA_Root_AllYears.png
│   │   └──🖼 PCoA_Soil_AllYears.png
│   ├──🖼 Decontaminated_plot.png                 # Contaminant filtering summary plot
│   ├──📁 Differential_abundance                  # Volcano plots for peg tissue by year
│   │   ├──🖼 Diff_abundance_2022_peg2.png
│   │   ├──🖼 Diff_abundance_2023_peg.png
│   │   └──🖼️ Diff_abundance_2024_peg.png
│   ├──🖼 fungi.rareplot.png                      # Rarefaction curve
│   └──🖼 read.depths.plot.png                    # Histogram of read depths
├──📝 README.md                                # Project description notes
└──📁 Tables                                      # Statistical test outputs
    ├──📁 Alpha_diversity                         # ANOVA result tables for alpha diversity metrics
    │   ├──📄 even_2022.doc
    │   ├──📄 even_2023.doc
    │   ├──📄 even_2024.doc
    │   ├──📄 invsimpson_2022.doc
    │   ├──📄 invsimpson_2023.doc
    │   ├──📄 invsimpson_2024.doc
    │   ├──📄 richness_2022.doc
    │   ├──📄 richness_2023.doc
    │   ├──📄 richness_2024.doc
    │   ├──📄 shannon_2022.doc
    │   ├──📄 shannon_2023.doc
    │   └──📄 shannon_2024.doc
    ├──📁 Anosim                              # ANOSIM result tables
    │   ├──📄 ANOSIM_2022_Peg.doc
    │   ├──📄 ANOSIM_2022_Root.doc
    │   ├──📄 ANOSIM_2022_Soil.doc
    │   ├──📄 ANOSIM_2023_Peg.doc
    │   ├──📄 ANOSIM_2023_Root.doc
    │   ├──📄 ANOSIM_2023_Soil.doc
    │   ├──📄 ANOSIM_2024_Peg.doc
    │   ├──📄 ANOSIM_2024_Root.doc
    │   └──📄 ANOSIM_2024_Soil.doc
    ├──📁 Beta_dispersion                       # Betadispersion results
    │   ├──📄 BetaDispersion_2022_Peg.doc
    │   ├──📄 BetaDispersion_2022_Root.doc
    │   ├──📄 BetaDispersion_2022_Soil.doc
    │   ├──📄 BetaDispersion_2023_Peg.doc
    │   ├──📄 BetaDispersion_2023_Root.doc
    │   ├──📄 BetaDispersion_2023_Soil.doc
    │   ├──📄 BetaDispersion_2024_Peg.doc
    │   ├──📄 BetaDispersion_2024_Root.doc
    │   └──📄 BetaDispersion_2024_Soil.doc
    ├──📁 Global_permanova                      # Global PERMANOVA by year
    │   ├──📄 Permanova_2022.doc
    │   ├──📄 Permanova_2023.doc
    │   └──📄 Permanova_2024.doc
    └──📁 Permanova                             # Tissue-level PERMANOVA tests
        ├──📄 Permanova_2022_Peg.doc
        ├──📄 Permanova_2022_Root.doc
        ├──📄 Permanova_2022_Soil.doc
        ├──📄 Permanova_2023_Peg.doc
        ├──📄 Permanova_2023_Root.doc
        ├──📄 Permanova_2023_Soil.doc
        ├──📄 Permanova_2024_Peg.doc
        ├──📄 Permanova_2024_Root.doc
        └──📄 Permanova_2024_Soil.doc
        
```
