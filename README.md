# B-RISA

Reproducibility materials for the study:

**Does Bayesian Posterior Uncertainty Improve Selective Clinical Prediction? Evidence from Geographic Validation of Appendicitis Models**

B-RISA (Bayesian Robust and Selective Artificial Intelligence) is an evaluation framework for testing whether posterior-derived uncertainty provides incremental decision value beyond confidence already encoded in the predicted probability. The study evaluates Bayesian logistic regression (BLR), Bayesian linear discriminant analysis (BLDA), and Bayesian additive regression trees (BART) using clinical appendicitis data, stationary simulation, nested selective-prediction analysis, locked geographic validation, and a paired source-to-target shift experiment.

## Repository structure

```text
B-RISA/
├── data/
│   ├── derived/
│   └── raw/
├── R/
├── results/
│   ├── adult/
│   ├── dusseldorf/
│   ├── nested_selector/
│   ├── paired_shift/
│   ├── rpa/
│   └── stationary_simulation/
├── .gitattributes
├── .gitignore
└── README.md
```

The `data/derived/` directory contains analysis-ready inputs and empirical parameter tables used by the released scripts. Raw clinical datasets are not redistributed in this repository.

The `results/` directory contains selected numerical outputs used to verify the manuscript results. Large stationary-simulation result files are stored using Git Large File Storage (Git LFS).

Scripts write newly generated computational outputs to separate output directories so that the committed manuscript-verification results are not overwritten during ordinary reproduction.

## Data sources

### Adult Appendicitis

The Adult Appendicitis dataset is publicly available from Kaggle:

https://www.kaggle.com/datasets/timrie/appendicitis

The analysis-ready predictors used in this study are included in `data/derived/adult_selected_cleaned_data.csv`.

### Regensburg Pediatric Appendicitis (RPA)

The RPA dataset is publicly available from Zenodo:

https://zenodo.org/records/7669442

The analysis-ready RPA dataset and the empirical parameters used for simulation calibration are included in `data/derived/`.

### Düsseldorf external validation

The independent Düsseldorf pediatric appendicitis dataset is available from:

https://github.com/i6092467/pediatric-appendicitis-ml-ext

The raw external dataset is not redistributed here. Download `app_data_ext.csv` from the source repository and place it at:

```text
data/raw/app_data_ext.csv
```

before running the external-validation or nested-selector scripts.

## Analysis scripts

The `R/` directory contains the seven analysis scripts used for the study.

```text
01_generate_stationary_simulation.R
02_analyze_stationary_simulation.R
03_adult_clinical_analysis.R
04_rpa_clinical_analysis.R
05_dusseldorf_external_validation.R
06_nested_incremental_selector.R
07_paired_source_target_shift.R
```

### 01_generate_stationary_simulation.R

Generates the 576-scenario clinically calibrated factorial simulation with 100 repetitions per scenario. The design varies generated sample size, predictor dimension, prevalence, boundary geometry, covariance structure, label noise, and predictor missingness.

### 02_analyze_stationary_simulation.R

Fits BLR, BLDA, and BART to the generated stationary simulation datasets and computes predictive performance, calibration, uncertainty, selective-prediction, and statistical-comparison outputs.

### 03_adult_clinical_analysis.R

Runs the supporting Adult Appendicitis analysis using 50 repeated stratified 10-fold cross-validation.

### 04_rpa_clinical_analysis.R

Runs the primary RPA base-model analysis using 50 repeated stratified 10-fold cross-validation with training-fold preprocessing.

### 05_dusseldorf_external_validation.R

Fits the final Bayesian models on RPA and applies the frozen RPA-trained models and preprocessing to the independent Düsseldorf cohort. The documented WBC harmonization rule converts values greater than 1000 by dividing by 1000.

### 06_nested_incremental_selector.R

Runs the fully nested RPA incremental-selector analysis. Ten outer folds define untouched test patients. Within each outer-training set, 50 repeated 10-fold inner out-of-fold predictions are used to train the margin-only and combined error-risk selectors. The script also performs the dependency audit, locked Düsseldorf selector application, WBC sensitivity analysis, false-negative capture analysis, and illustrative asymmetric-cost analyses.

### 07_paired_source_target_shift.R

Runs the 500-replicate paired source-to-target experiment. Each model is fitted once per RPA-like source sample and evaluated under no shift, prevalence shift, neutrophil-missingness shift, and combined shift.

## Main analysis settings

The principal reproducibility settings are:

- stationary simulation generation seed: `20260616`
- primary Bayesian fitting seed: `20260618`
- clinical statistical-analysis seed: `20260619`
- locked external-fit seed: `20260730`
- nested RPA seed: `20261115`
- nested patient-bootstrap seed: `20261116`
- RPA outer folds: `10`
- inner selector construction: `50` repeated `10`-fold cross-validation
- posterior draws per model: `300`
- BART trees: `50`
- BART burn-in: `150`
- primary classification threshold: `0.50`
- primary retained coverage: `0.80`
- selector ridge penalty: `1e-4`
- patient-bootstrap replicates: `2000`
- illustrative FN:FP cost ratios: `2:1` and `5:1`
- paired source-to-target replicates: `500`
- paired source sample size: `780`
- paired target sample size: `301`

## Recommended execution order

Run the scripts from the repository root.

```text
R/01_generate_stationary_simulation.R
R/02_analyze_stationary_simulation.R
R/03_adult_clinical_analysis.R
R/04_rpa_clinical_analysis.R
R/05_dusseldorf_external_validation.R
R/06_nested_incremental_selector.R
R/07_paired_source_target_shift.R
```

The stationary simulation is computationally intensive and produces large intermediate datasets. The committed `results/stationary_simulation/` files provide manuscript-verification outputs, so reproducing the full stationary simulation is not required simply to inspect the reported results.

## Software environment

The analyses reported in the manuscript were run in R 4.5.3. Recorded package versions include:

```text
BART       2.9.10
pROC       1.19.0.1
MASS       7.3-65
Matrix     1.7-4
data.table 1.18.2.1
dplyr      1.2.1
tidyr      1.3.2
ggplot2    4.0.2
```

The analysis scripts contain deterministic seeds for the principal experiments and save audit information needed to check model fitting and validation dependencies.

## Key reproducibility outputs

The committed results include:

- Adult and RPA cross-validation predictions and performance summaries;
- locked Düsseldorf model-performance and preprocessing audits;
- nested RPA patient predictions and dependency-audit output;
- nested Risk@80%, risk-coverage, error-prediction, and bootstrap contrasts;
- selector coefficients and out-of-fold manifests;
- WBC harmonization and sensitivity results;
- false-negative capture and illustrative asymmetric-cost results;
- paired source-to-target replication summaries and Monte Carlo confidence intervals;
- stationary simulation metrics, summaries, and inferential audit outputs.

## Interpretation

The repository supports the analyses reported in the manuscript. B-RISA is an evaluation framework rather than a new Bayesian classifier. Its primary question is whether posterior probability standard deviation improves selective prediction beyond confidence margin at matched coverage.

For BLR, BLDA, and BART in the evaluated appendicitis datasets, posterior probability SD did not demonstrate additional benefit over confidence margin for minimizing overall selective 0--1 error at a threshold of 0.50. The analysis does not establish that Bayesian uncertainty is generally ineffective, nor does it establish that one deferral strategy is clinically safer under every loss function or clinical workflow.

## Citation

If you use this code or the accompanying results, please cite the associated manuscript. Full publication information will be added after publication.

## Contact

Kazi Nur Uddin  
Department of Mathematical Sciences  
Kent State University  
Email: kuddin@kent.edu
