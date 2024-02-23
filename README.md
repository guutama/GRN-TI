# GRN-TI: Gene Expression Prediction Using Gene Regulatory Networks

## Introduction

The GRN-TI (Gene Regulatory Network - Transcriptome Imputations) pipeline is designed for predicting gene expression levels in the context of transcriptome-wide association studies (TWAS). This pipeline incorporates a novel approach that leverages both local and distal genetic variants through gene regulatory networks (GRNs) in line with the omnigenic model of complex trait inheritance.

## Workflow Overview
Our pipeline consists of several distinct and interconnected stages, outlined in Figure 1 (see fig:gene_snp_matrices).

 
1. ### <span style="color:red;">Data Preprocessing</span>

    > **SNP Genotype Data Conversion**

    Converting raw genotype data in vcf format into a matrix of SNP genotypes. The resulting matrix contains values of 0, 1, 2, or 'NA', corresponding to the major homozygous genotype, the heterozygous genotype, the minor homozygous genotype, and missing data, respectively.

    ##### Conversion Tool and commnad: PLINK2
    We utilize the PLINK2 software for this conversion of snp genotype in the geuvadis study using the following command and flags.

    `plink2 --vcf "path/to/vcf.gz" --export Av --extract "path/file/snps.txt" \
        --update-name "path/to/idconvert.txt" 1 2 --geno 0.05 --out "/path/filename"
    `
    Key Flags Explained:

    * --export Av: This flag is used for converting the genotype data to biallelic SNP genotypes.

    * --extract: Utilize this flag to select a specific list of SNPs. For our purposes, we only include the SNPs found in the Geuvadis eQTL analysis results. The file snps.txt should list one SNP per line, and each SNP must be unique.

    * --update-name: This flag facilitates the conversion of genotype IDs to rs_id. The conversion file should contain two columns per line, the first being the rsid and the second the ID present in the file.

    * --geno: Sets the threshold for the minor allele frequency.

    * Output File
    The output of the command will be a file named `filename.traw` located in the specified path (/path/filename.traw). 

    > **Standardizing Gene and SNP ID Column Names Across Datasets**

        - The `standardize_genomic_column_names.py` script ensures consistency in the naming conventions of gene and SNP ID columns across different datasets (gene expression, SNP genotype, and mapping data containing eQTL statistics)

        Input Requirements
        To utilize this script, the following inputs are required:

        * CSV/TSV Files: The three datasets must be in either CSV (Comma-Separated Values) or TSV (Tab-Separated Values) formats. These files should be located within the directory specified by RAW_CSV_TXT_PATH in the path_config.py file and have name `expression.csv/tsv`,`genotype.csv/tsv` and `mapping.csv/tsv`.

        Configuration File: The column names to be standardized across the datasets need to be defined in  `config.yaml` file. This configuration includes:

        expression_gene_id_col_name: Specifies the current gene ID column name in the expression dataset.
        mapping_gene_id_col_name: Specifies the gene ID column name in mapping dataset and to be used across both the expression and mapping datasets.
        genotype_snp_id_col_name: Specifies the current SNP ID column name in the genotype dataset.
        mapping_snp_id_col_name: Specifies the  SNP ID column name in mapping dataset and to be used across both the genotype and mapping datasets.

        Output Description
        * Upon successful execution, the script outputs the processed datasets to the directory specified by RAW_CSV_GZ_PATH in the `path_config.py` file. Each dataset is compressed using the GZip format to minimize storage space while maintaining data integrity. The output files include:

        Compressed expression dataset with the same gene ID column name as the in mapping.
        Compressed genotype dataset with the same SNP ID column name as in mapping dateset.
        Compressed mapping dataset

	 > **Aligning genotype and expression data**

        Gene Expression Data: We begin with raw gene expression data for a specific set of genes (denoted A to F).
        Identification of cis-eQTLs: Top cis-eQTLs (E1 to E6) corresponding to these genes are identified, with matching colors between eQTLs and genes indicating their correspondence.

2. ### <span style="color:red;">Bayesian Inference and Network Reconstruction </span>

    Bayesian Posterior Probabilities: Utilizing the Findr-tool, we obtain Bayesian posterior probabilities, estimating the likelihood of interactions between genes.
    Network Reconstruction: These probabilities are used to reconstruct the GRN, employing a directed acyclic graph (DAG) structure.

3. ### <span style="color:red;">Model Development and Deployment</span>
    Prediction Model Training: Utilizing the GRN structure, we train a prediction model for each gene, incorporating both the gene's cis-eQTLs and the data from its parent genes in the network.
    Iterative Process: The prediction model undergoes an iterative refinement process, using gene expression and predicted gene expression inputs, represented by dashed lines in the workflow.
    Model Deployment: In the deployment phase, personal genotypes are used as inputs to predict individualized gene expression levels.
    Figure Description (Figure 1: fig:gene_snp_matrices)
    The figure provides a visual representation of our method, showcasing the raw data inputs, the Findr-tool process, and the subsequent stages leading to the reconstruction of the GRN and the development of the prediction model. Key elements include:

    Color Coding: Indicates the correspondence between genes and eQTLs, and the types of inputs used in different stages of the model.
    Arrows and Lines: Denote the flow of information and the relationships between different components of the pipeline.
    Getting Started
    Installation
    Instructions on how to install and set up the GRN-TI Pipeline.

##  Usage
    A detailed guide on how to use the pipeline, including command examples and configurations.

## Additional Sections
    Requirements: List of dependencies and environmental requirements.
    Contributing: Guidelines for contributing to the project.
    License: Information about the licensing of the GRN-TI Pipeline.
    Contact: Details for reaching out to the maintainers or contributors.
