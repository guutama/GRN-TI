# GRN-TI: Gene Expression Prediction Using Gene Regulatory Networks

## Introduction

The GRN-TI (Gene Regulatory Network - Transcriptome Imputations) pipeline is designed for predicting gene expression levels in the context of transcriptome-wide association studies (TWAS). This pipeline incorporates a novel approach that leverages both local and distal genetic variants through gene regulatory networks (GRNs) in line with the omnigenic model of complex trait inheritance.

## Workflow Overview

Our pipeline consists of several distinct and interconnected stages, outlined in Figure 1 (see fig:gene_snp_matrices).

### **Data Preprocessing**

> **SNP Genotype Data Conversion**

Converting raw genotype data in vcf format into a matrix of SNP genotypes. The resulting matrix contains values of 0, 1, 2, or 'NA', corresponding to the major homozygous genotype, the heterozygous genotype, the minor homozygous genotype, and missing data, respectively.

**Conversion Tool and Command: PLINK2**
We utilize the PLINK2 software for this conversion of SNP genotype in the Geuvadis study using the following command and flags.

plink2 --vcf "path/to/vcf.gz" --export Av --extract "path/file/snps.txt"
--update-name "path/to/idconvert.txt" 1 2 --geno 0.05 --out "/path/filename"


**Key Flags Explained:**

- `--export Av`: This flag is used for converting the genotype data to biallelic SNP genotypes.
- `--extract`: Utilize this flag to select a specific list of SNPs. For our purposes, we only include the SNPs found in the Geuvadis eQTL analysis results. The file snps.txt should list one SNP per line, and each SNP must be unique.
- `--update-name`: This flag facilitates the conversion of genotype IDs to rs_id. The conversion file should contain two columns per line, the first being the rsid and the second the ID present in the file.
- `--geno`: Sets the threshold for the minor allele frequency.

**Output File**
The output of the command will be a file named `filename.traw` located in the specified path (/path/filename.traw).

> **Standardizing Gene and SNP ID Column Names Across Datasets**

- The `standardize_column_names.py` script ensures consistency in the naming conventions of gene and SNP ID columns across different datasets (gene expression, SNP genotype, and mapping data containing eQTL statistics)

**Input Requirements**
To utilize this script, the following inputs are required:

- **CSV/TSV Files**: The three datasets must be in either CSV (Comma-Separated Values) or TSV (Tab-Separated Values) formats. These files should be located within the directory specified by RAW_CSV_TXT_PATH in the path_config.py file and have name `expression.csv/tsv`,`genotype.csv/tsv` and `mapping.csv/tsv`.
- **Configuration File**: The column names to be standardized across the datasets need to be defined in  `config.yaml` file. This configuration includes:
    - `expression_gene_id_col_name`: Specifies the current gene ID column name in the expression dataset.
    - `mapping_gene_id_col_name`: Specifies the gene ID column name in mapping dataset and to be used across both the expression and mapping datasets.
    - `genotype_snp_id_col_name`: Specifies the current SNP ID column name in the genotype dataset.
    - `mapping_snp_id_col_name`: Specifies the SNP ID column name in mapping dataset and to be used across both the genotype and mapping datasets.

**Output Description**
- Upon successful execution, the script outputs the processed datasets to the directory specified by RAW_CSV_GZ_PATH in the `path_config.py` file. Each dataset is compressed using the GZip format to minimize storage space while maintaining data integrity. The output files include:
    - Compressed expression dataset with the same gene ID column name as in mapping.
    - Compressed genotype dataset with the same SNP ID column name as in mapping dataset.
    - Compressed mapping dataset

### **Aligning Genotype and Expression Data**

In order to reconstruct the gene regulatory network using FINDR, aligning expression and genotype data is a crucial step. The `align_expression_genotype.py` script is specifically designed for this purpose, ensuring that genotype and expression datasets are precisely aligned to facilitate further analyses.

#### Input Format

The script expects input data in compressed CSV (.csv.gz) format, tailored to accommodate the following data types:

- **Expression Data (`expression.csv.gz`)**: Contains gene expression data, featuring at least two essential columns — one for gene identifiers (`GENE_ID`) and others for expression levels across various samples or conditions.

- **Genotype Data (`genotype.csv.gz`)**: Comprises SNP genotype information, including columns for SNP identifiers (`SNP_ID`) and genotype data across samples. The genotype data typically represent allele counts (0, 1, 2) for each SNP per sample.

- **Mapping Data (`mapping.csv.gz`)** (Optional): This file bridges genes to SNPs, offering insights into eQTLs. It includes columns for gene identifiers (`GENE_ID`), SNP identifiers (`SNP_ID`), and potentially statistical measures like the r-value, showcasing the gene-SNP association strength.

#### Output Format

After processing the input data, the script outputs several files in compressed CSV format, detailed as follows:

- **Top SNP Data (`top_snp.csv.gz`)**: Features genotype data for SNPs most strongly associated with genes, according to the mapping data. This file filters the genotype data to include only these top SNPs.

- **All SNP Data (`all_snp.csv.gz`)**: A comprehensive dataset that encompasses genotype information for all SNPs in the input genotype file, irrespective of their gene associations.

- **Top Mapping Data (`top_mapping.csv.gz`)**: Contains filtered mapping data for top-associated SNPs and genes, selecting only entries for SNPs with the highest r-values per gene.

- **All Mapping Data (`all_mapping.csv.gz`)**: Includes the entire set of mapping data from the input file without filtration, catering to analyses requiring a complete overview of gene-SNP associations.

- **Top Genes Data (`top_genes.csv.gz`)**: Consists of expression data for genes linked with SNPs in the top mapping data, focusing the analysis on genes with significant eQTLs.

- **All Genes Data (`all_genes.csv.gz`)**: Contains expression data for all genes in the input expression file, providing a dataset for more extensive analyses.

This alignment is fundamental for the GRN-TI pipeline's success, setting the stage for accurate gene regulatory network reconstruction and subsequent gene expression prediction.


### **Bayesian Inference and Network Reconstruction**

- Bayesian Posterior Probabilities: Utilizing the Findr-tool, we obtain Bayesian posterior probabilities, estimating the likelihood of interactions between genes.
- Network Reconstruction: These probabilities are used to reconstruct the GRN, employing a directed acyclic graph (DAG) structure.

### **Model Development and Deployment**

- Prediction Model Training: Utilizing the GRN


