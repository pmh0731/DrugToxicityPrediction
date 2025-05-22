# DrugToxicityPrediction

`drug_toxicity_predictor.py` is a Python script designed for evaluating drug toxicity based on biological target features (GPD) and chemical structure. It supports both performance evaluation on a dataset used in the paper (Park et al., 2025) and prediction for new drug candidates.

## Setup Instructions

Before running the code, create a Conda virtual environment using the provided `DrugToxicityPrediction.yaml` file. This ensures that all required dependencies (including R integration and machine learning libraries) are properly configured.

```bash
git clone https://github.com/pmh0731/DrugToxicityPrediction.git
cd ./DrugToxicityPrediction

conda env create -f DrugToxicityPrediction.yaml
conda activate DrugToxicityPrediction
```

## Script: `drug_toxicity_predictor.py`

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--mode` | Yes | Operation mode: `e` for evaluation using the dataset from the publication, `p` for predicting toxicity of new drugs |
| `--input` | Conditional | Path to the input file (required only for prediction mode `--mode p`) |
| `--hg` | No | Ensembl genome version of human target genes. Default is `hg38` |
| `--output` | Yes | Output directory path where results will be saved |

### Usage

#### 1. Reproduce Evaluation Results from the Study

This mode evaluates model performance using the dataset provided in `Supplementary_Table_3.xlsx` (from the original publication).

```bash
python ./code/drug_toxicity_predictor.py --mode e --output ./results_evaluation
```

This will generate a performance metrics file (`evaluation_performance.xlsx`) in the specified output directory.

#### 2. Predict Toxicity for New Drug Candidates

To predict the toxicity status of new compounds, prepare an input file containing target gene information and chemical structures. Use the format provided in `test_data.tsv` as a reference.

```bash
python ./code/drug_toxicity_predictor.py --mode p --input ./test_data.tsv --output ./results_prediction
```

Optional: specify a different Ensembl human genome version (default is `hg38`):

```bash
python ./code/drug_toxicity_predictor.py --mode p --input ./test_data.tsv --hg hg19 --output ./results_prediction
```

The prediction results, including toxicity probabilities and status (`Approved drug` or `Risky drug`), will be saved to the output directory.

##### Note: Predictions are not generated for compounds if their human or mouse target gene features have less than 50% coverage, or if SMILES featurization cannot be performed. For any drug that is excluded from prediction, the output file will indicate the specific reason.

### Notes

- The `input` argument is **only required for prediction mode** (`--mode p`).
- The Ensembl genome version specified by `--hg` affects mapping of human target genes.
- Ensure R is installed and accessible, as the script uses R functions via `rpy2`.

## Output Files

Depending on the mode, the following output files will be generated:

- **Evaluation mode (`e`)**: `evaluation_performance.xlsx`
- **Prediction mode (`p`)**: `[input_filename].xlsx` containing predicted toxicity outcomes
