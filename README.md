# GRN-TI Pipeline: Gene Expression Prediction Using Gene Regulatory Networks

## Introduction

The GRN-TI (Gene Regulatory Network - Transcriptome Imputations) Pipeline is designed for predicting gene expression levels in the context of transcriptome-wide association studies (TWAS). This pipeline incorporates a novel approach that leverages both local and distal genetic variants through gene regulatory networks (GRNs) in line with the omnigenic model of complex trait inheritance.

## Workflow Overview
Our pipeline consists of several distinct and interconnected stages, outlined in Figure 1 (see fig:gene_snp_matrices).

 1. ### <span style="color:red;">Data Preprocessing</span>

    > SNP Genotype Data Conversion

        Converting raw genotype data in vcf format into a matrix of SNP genotypes. The resulting matrix contains values of 0, 1, 2, or 'NA', corresponding to the major homozygous genotype, the heterozygous genotype, the minor homozygous genotype, and missing data, respectively.

        ##### Conversion Tool: PLINK2
        We utilize the PLINK2 software for this conversion process.

        ##### Example Command:
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


    > Aligning genotype and expression data


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
