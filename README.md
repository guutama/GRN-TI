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
│   │   ├── <span style="color: grey;">geuvadis</span>
│   │   │   ├── script.py
│   │   ├── <span style="color: grey;">dream</span>
│   │   │   ├── script.py
│   │   └── <span style="color: grey;">yeast</span>
│   │       ├── script.py
│   ├── <span style="color: skyblue;">structure_learning</span>
│   │   ├── <span style="color: grey;">geuvadis</span>
│   │   │   ├── script.py
│   │   ├── <span style="color: grey;">dream</span>
│   │   │   ├── script.py
│   │   └── <span style="color: grey;">yeast</span>
│   │       ├── script.py
│   └── <span style="color: darkgreen;">modeling</span>
│       ├── <span style="color: grey;">geuvadis</span>
│       │   └── script.py
│       ├── <span style="color: grey;">dream</span>
│       │   └── script.py
│       └── <span style="color: grey;">yeast</span>
│           └── script.py
│   ├── param_config.yaml
│   ├── path_config.py
│   └── utils.py
├── <span style="color: blue;">data</span>
│   ├── <span style="color: yellow;">raw</span>
│   │   ├── <span style="color: orange;">geuvadis</span>
│   │   ├── <span style="color: orange;">dream</span>
│   │   └── <span style="color: orange;">yeast</span>
│   └── <span style="color: orange;">processed</span>
│       ├── <span style="color: orange;">aligned</span>
|            └── [output directory  for preprocessing]
│       ├── <span style="color: orange;">split</span>
|            └── [output directory  for preprocessing]
│       ├── <span style="color: orange;">pairwise_inference</span>
|            └── [output directory  for structure learning]
│       ├── <span style="color: orange;">network_inference</span>
|            └── [output directory  for structure learning]
├── <span style="color: gold;">metrics</span>
|    └── [output directory  for modelling]
└── <span style="color: brown;">models</span>
|     └── [output directory  for modelling]
</code>
</pre>



































