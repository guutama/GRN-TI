# Predicting the genetic component of gene expression using gene regulatory networks

## Introduction

This README presents a detailed overview of our  pipeline for building prediction method and model for genetic component of gene expression. This document provides an overview of our data transformation  stages of preprocessing, structure learning, and modeling. 

---

## Pipeline Stages 

### Preprocessing Stage

The preprocessing stage prepares the raw data  by aligning expressions with genotypes and splitting the dataset into training and test sets.

### Structure Learning Stage

In the structure learning stage, the pipeline infers gene regulatory network and make input-output features for each gene in the network

### Modeling Stage

The final stage involves training the model on the featurized data and evaluating its performance.




## Path Configuration

Our pipeline utilizes a path configuration file `/src/path_config.py` where all necessary paths for data  are defined. 

Each Python file within our pipeline refers to these path configurations to read input data and store outputs. Here is an example demonstrating how scripts use the path configuration:
```python
#Example content of path_config.py file
from pathlib import Path
MY_PATH = Path("/cluster/projects/nn1015k/GRN-TI")  # TODO: Needs to be updated
DATA_PATH = MY_PATH / "data"
RAW_CSV_GZ_PATH = DATA_PATH / "raw_csv_gz" # Path to the raw CSV GZ files.
ALIGNED_PATH = DATA_PATH / "aligned" # Path where aligned data is stored.
SPLIT_PATH = DATA_PATH / "split" # Path for the split data.
PAIRWISE_PROBABILITY_PATH =  DATA_PATH / "pairwise_probability" # Path for the output of pairwise inference.
```

```python
#How the align_expression_genotype.py file called
if __name__ == '__main__':
    in_path = RAW_CSV_GZ_PATH
    out_path = ALIGNED_PATH
    if not out_path.exists():
        out_path.mkdir(parents=True)
    align_data(in_path=in_path,out_path=out_path)   
```
### Parameters Configuration

Additionally, the pipeline uses a `param_config.yaml` file to define all variables and parameters needed at each stage of the pipeline. This configuration file is structured to include settings for preprocessing, structure learning, and modeling stages, allowing for easy adjustments to the pipeline's behavior without modifying the code directly.

```yaml
preprocessing:
  align:
  split:
    test_size: 0.2
structure_learning:
  pairwise_inference:
  network_inference:
    fdr_prior: 0.6
  featurize:

modeling:
  train: 
    learning_rate: 0.01
```
```python
#How it is used in split_train_test.py
import yaml

params = yaml.safe_load(open("src/param_config.yaml"))["preprocessing"]
params_split = params['split']
test_size = float(params_split['test_size'])
```

## Pipeline Structure and Workflow 

The flow of the pipeline is shown below. Each stage depends on the next stage solely through the data file, meaning if the data is provided in the correct format, the code for each stage will be independent of the others.

<pre>
<code>
<span style="color: black;">PROJECT ROOT</span>
├── <span style="color: green;">src</span>
│   ├── <span style="color: lightgreen;">preprocessing</span>
│   │   ├── <span style="color: grey;">align_expression_genotype()</span> <!-- <br> -->
│   │   │   └── Uses data from: <span style="color: blue;">/data/raw_csv_gz</span> <!-- <br> -->
│   │   │   └── Outputs to: <span style="color: orange;">/data/aligned</span> <!-- <br> -->
│   │   │   <!-- <br> -->
│   │   └── <span style="color: grey;">split_train_test()</span> <!-- <br> -->
│   │       └── Uses data from: <span style="color: orange;">/data/aligned</span> <!-- <br> -->
│   │       └── Outputs to: <span style="color: purple;">/data/split</span> <!-- <br> -->
│   │       <!-- <br> -->
│   ├── <span style="color: skyblue;">structure_learning</span> <!-- <br> -->
│   │   ├── <span style="color: grey;">pairwise_inference()</span> <!-- <br> -->
│   │   │   └── Uses data from: <span style="color: purple;">/data/split</span> <!-- <br> -->
│   │   │   └── Outputs to: <span style="color: pink;">/data/pairwise_probability</span> <!-- <br> -->
│   │   │   <!-- <br> -->
│   │   ├── <span style="color: grey;">network_inference()</span> <!-- <br> -->
│   │   │   └── Uses data from: <span style="color: pink;">/data/pairwise_probability</span> <!-- <br> -->
│   │   │   └── Outputs to: <span style="color: red;">/data/networks</span> <!-- <br> -->
│   │   │   <!-- <br> -->
│   │   └── <span style="color: grey;">featurize()</span> <!-- <br> -->
│   │       └── Uses data from: <span style="color: red;">/data/networks</span> <!-- <br> -->
│   │       └── Outputs to: <span style="color: darkred;">/data/featurized</span> <!-- <br> -->
│   │       <!-- <br> -->
│   └── <span style="color: darkgreen;">modeling</span> <!-- <br> -->
│       └── <span style="color: grey;">train_evaluate()</span> <!-- <br> -->
│           └── Uses data from: <span style="color: darkred;">/data/featurized</span> <!-- <br> -->
│           └── Outputs to: <span style="color: gold;">/data/metrics</span> <!-- <br> -->
│           <!-- <br> -->
└── <span style="color: blue;">data</span> <!-- <br> -->
    ├── <span style="color: yellow;">raw_csv_gz</span> <!-- <br> -->
    │   ├── expression.csv.gz <!-- <br> -->
    │   ├── genotype.csv.gz <!-- <br> -->
    │   └── mapping.csv.gz <!-- <br> -->
    │   <!-- <br> -->
    ├── <span style="color: orange;">aligned</span> <!-- <br> -->
    │   └── [Generated by align_expression_genotype()] <!-- <br> -->
    │   <!-- <br> -->
    ├── <span style="color: purple;">split</span> <!-- <br> -->
    │   └── [Generated by split_train_test()] <!-- <br> -->
    │   <!-- <br> -->
    ├── <span style="color: pink;">pairwise_probability</span> <!-- <br> -->
    │   └── [Generated by pairwise_inference()] <!-- <br> -->
    │   <!-- <br> -->
    ├── <span style="color: red;">networks</span> <!-- <br> -->
    │   └── [Generated by network_inference()] <!-- <br> -->
    │   <!-- <br> -->
    ├── <span style="color: darkred;">featurized</span> <!-- <br> -->
    │   └── [Generated by featurize()] <!-- <br> -->
    │   <!-- <br> -->
    └── <span style="color: gold;">metrics</span> <!-- <br> -->
        └── [Generated by train_evaluate()] <!-- <br> -->
</code>
</pre> 






<!-- # Predicting the genetic component of gene expression using gene regulatory networks

<!-- ## Introduction

The GRN-TI (Gene Regulatory Network - Transcriptome Imputations) pipeline is designed for predicting gene expression levels in the context of transcriptome-wide association studies (TWAS). This pipeline incorporates a novel approach that leverages both local and distal genetic variants through gene regulatory networks (GRNs) in line with the omnigenic model of complex trait inheritance. -->
<!-- 
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
    - Compressed expression dataset with the GRN-TI: Gene Expression Prediction Using Gene Regulatory Networks> **Aligning Genotype and Expression Data**

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



> **Data Splitting for Training and Testing**

 The `split_expression_genotype.py` splits data into training and test


#### Input Format

The script expects the following input files in compressed CSV format:

- **Top SNP Data (`top_snp.csv.gz`)**: Genotype data for SNPs strongly associated with genes.
- **All SNP Data (`all_snp.csv.gz`)**: Genotype data for all SNPs.
- **Top Mapping Data (`top_mapping.csv.gz`)**: Mapping data for top-associated SNPs and genes.
- **Top Genes Data (`top_genes.csv.gz`)**: Expression data for genes with significant eQTLs.
- **All Genes Data (`all_genes.csv.gz`)**: Expression data for all genes.
- **Sample Names (`sample_names.csv.gz`)**: A list of sample names containing sample both in SNP and Gene data.

#### Output Format

The script generates the following output files in compressed CSV format:

- **Training and Testing Expression Data**: Separate files for top genes (`top_exp_train.csv.gz`, `top_exp_test.csv.gz`) and all genes (`all_exp_train.csv.gz`, `all_exp_test.csv.gz`).
- **Training and Testing Genotype Data**: Separate files for top SNPs (`top_eqtl_train.csv.gz`, `top_eqtl_test.csv.gz`) and all SNPs (`all_eqtl_train.csv.gz`, `all_eqtl_test.csv.gz`).
- **Sample Lists**: Lists of sample names used in the training (`train_sample.csv.gz`) and testing (`test_sample.csv.gz`) datasets.











### **Bayesian Inference and Network Reconstruction**

- Bayesian Posterior Probabilities: Utilizing the Findr-tool, we obtain Bayesian posterior probabilities, estimating the likelihood of interactions between genes.
- Network Reconstruction: These probabilities are used to reconstruct the GRN, employing a directed acyclic graph (DAG) structure.

### **Model Development and Deployment**

- Prediction Model Training: Utilizing the GRN

  -->
