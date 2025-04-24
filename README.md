cat("
# Peanut Fungal ITS Amplicon Sequencing

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


[**Flourishing Fungi** pipeline by Zachary Noel (2024)](https://zenodo.org/doi/10.5281/zenodo.10655178)

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
   
   -- Note: This step is not yet publicly available in the original Flourishing Fungi pipeline
   by Dr. Zachary Noel, but is included here with permission. It is based on internal, unpublished
   code developed by Dr. Noel for CONSTAX2 integration on ASAX. This functionality may be included
   in future public releases of the original pipeline.
")
