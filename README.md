# B-RISA

Reproducibility materials for:

**Does Bayesian Posterior Uncertainty Improve Selective Clinical Prediction? Evidence from Geographic Validation of Appendicitis Models**

B-RISA (Bayesian Robust and Selective Artificial Intelligence) is an evaluation framework for testing whether posterior-derived uncertainty provides incremental decision value beyond confidence already encoded in a predicted probability. The study evaluates Bayesian logistic regression (BLR), Bayesian linear discriminant analysis (BLDA), and Bayesian additive regression trees (BART) using clinical appendicitis datasets, stationary simulation, nested selective-prediction analysis, locked geographic validation, and a paired source-to-target shift experiment.

## Repository structure

```text
B-RISA/
├── data/
│   ├── derived/
│   └── raw/
├── R/
│   ├── 01_generate_stationary_simulation.R
│   ├── 02_analyze_stationary_simulation.R
│   ├── 03_adult_clinical_analysis.R
│   ├── 04_rpa_clinical_analysis.R
│   ├── 05_dusseldorf_external_validation.R
│   ├── 06_nested_incremental_selector.R
│   ├── 07_paired_source_target_shift.R
│   └── 08_locked_external_selector_bootstrap.R
├── results/
│   ├── adult/
│   ├── dusseldorf/
│   ├── locked_external_selector/
│   ├── nested_selector/
│   ├── paired_shift/
│   ├── rpa/
│   └── stationary_simulation/
├── .gitattributes
├── .gitignore
└── README.md
```

The `data/derived/` directory contains analysis-ready inputs and empirical parameter tables used by the released scripts. Raw clinical datasets are not redistributed in this repository.

The `results/` directory contains selected numerical outputs used to verify the manuscript results. Large stationary-simulation result files are stored with Git Large File Storage (Git LFS).

Newly generated outputs are written under `analysis_outputs/` so that committed verification files in `results/` are not overwritten during ordinary reproduction.

## Data sources

### Adult Appendicitis

The Adult Appendicitis dataset is publicly available from Kaggle:

https://www.kaggle.com/datasets/timrie/appendicitis

The analysis-ready dataset used by the released scripts is included as:

```text
data/derived/adult_selected_cleaned_data.csv
```

### Regensburg Pediatric Appendicitis (RPA)

The RPA dataset is publicly available from Zenodo:

https://zenodo.org/records/7669442

The analysis-ready RPA dataset, standardization parameters, and empirical mean/covariance files used by the released analyses are included under:

```text
data/derived/
```

### Düsseldorf external validation

The independent Düsseldorf pediatric appendicitis dataset is publicly available from:

https://github.com/i6092467/pediatric-appendicitis-ml-ext

The raw external dataset is not redistributed here. Download:

```text
app_data_ext.csv
```

and place it at:

```text
data/raw/app_data_ext.csv
```

before running the locked external-validation analysis.

## Analysis scripts

### `01_generate_stationary_simulation.R`

Generates the 576-scenario clinically calibrated factorial simulation with 100 repetitions per scenario. The design varies generated sample size, predictor dimension, prevalence, boundary geometry, covariance structure, label noise, and predictor missingness.

### `02_analyze_stationary_simulation.R`

Fits BLR, BLDA, and BART to the generated stationary simulation datasets and computes predictive performance, calibration, uncertainty, selective-prediction, and statistical-comparison outputs.

### `03_adult_clinical_analysis.R`

Runs the supporting Adult Appendicitis analysis using 50 repeated stratified 10-fold cross-validation.

### `04_rpa_clinical_analysis.R`

Runs the RPA base-model analysis using 50 repeated stratified 10-fold cross-validation with training-fold preprocessing. It produces the repeated held-out predictions used in subsequent selector analyses.

### `05_dusseldorf_external_validation.R`

Fits the prespecified Bayesian models on RPA and applies frozen RPA-derived preprocessing and fitted models to the independent Düsseldorf cohort.

The documented WBC harmonization rule converts observed WBC values greater than 1000 by dividing by 1000. No Düsseldorf outcomes are used to fit, tune, recalibrate, or select the prediction models.

### `06_nested_incremental_selector.R`

Runs the fully nested RPA incremental-selector analysis used for the primary internal incremental-value inference.

Ten outer folds define untouched test patients. Within each outer-training set, 50 repeated 10-fold inner out-of-fold base-model predictions are used to train the margin-only and combined error-risk selectors. The script produces the nested dependency audit, Risk@80% contrasts, error-prediction metrics, false-negative capture, WBC sensitivity analysis, and illustrative asymmetric-cost analyses.

**This is the authoritative script for the primary RPA selector results.**

### `07_paired_source_target_shift.R`

Runs the 500-replicate paired source-to-target experiment. Each model is fitted once per RPA-like source sample and evaluated under four matched target conditions:

- no shift;
- prevalence shift;
- neutrophil-missingness shift;
- combined prevalence and missingness shift.

Source-derived preprocessing is applied unchanged to all targets.

### `08_locked_external_selector_bootstrap.R`

Reproduces the exploratory, post-hoc Düsseldorf selective-prediction comparisons reported in the manuscript.

The script:

1. reads the saved 50-repeat RPA held-out predictions;
2. constructs the transport-compatible RPA posterior SD as the square root of the mean within-fit posterior variance across repetitions;
3. fits the final ridge-logistic error-risk selector using RPA only;
4. applies the frozen selector to the locked Düsseldorf predictions;
5. computes posterior-SD, confidence-margin, and combined-selector selective risk;
6. computes risk-coverage curves and normalized AURC;
7. reproduces the paired 2,000-replicate Düsseldorf patient-bootstrap contrasts and Holm-adjusted p values;
8. checks the reproduced Düsseldorf values against the manuscript.

This script does **not** replace the fully nested RPA analysis. The primary RPA incremental-value results remain those from `06_nested_incremental_selector.R`.

## Recommended execution order

Run the scripts from the repository root or from the `R/` directory.

```text
R/01_generate_stationary_simulation.R
R/02_analyze_stationary_simulation.R
R/03_adult_clinical_analysis.R
R/04_rpa_clinical_analysis.R
R/05_dusseldorf_external_validation.R
R/06_nested_incremental_selector.R
R/07_paired_source_target_shift.R
R/08_locked_external_selector_bootstrap.R
```

Scripts 06 and 08 have different roles:

- **Script 06:** fully nested RPA evaluation of incremental selector value.
- **Script 08:** frozen RPA-only selector applied to Düsseldorf, including the exploratory external paired bootstrap.

The stationary simulation is computationally intensive and generates large intermediate datasets. The committed stationary-simulation outputs allow inspection of the reported results without rerunning all 57,600 simulated datasets.

## Main analysis settings

Key settings used in the manuscript include:

- stationary simulation data-generation seed: `20260616`
- primary Bayesian fitting seed: `20260618`
- clinical statistical-analysis seed: `20260619`
- locked external model-fit seed: `20260730`
- fully nested RPA seed: `20261115`
- nested RPA patient-bootstrap seed: `20261116`
- Düsseldorf exploratory selector-bootstrap seed: `20261031`
- paired source-to-target base seed: `20261101`
- RPA outer folds: `10`
- inner nested selector construction: `50` repeated `10`-fold cross-validation
- posterior draws per Bayesian fit: `300`
- BART trees: `50`
- BART burn-in: `150`
- classification threshold: `0.50`
- primary retained coverage: `0.80`
- selector ridge penalty: `1e-4`
- patient-bootstrap replicates: `2000`
- illustrative false-negative:false-positive cost ratios: `2:1` and `5:1`
- paired source-to-target replicates: `500`
- paired source sample size: `780`
- paired target sample size: `301`

The final full-RPA selector used for locked Düsseldorf application is deterministic conditional on the saved RPA held-out predictions. Seed `20261031` controls the exploratory Düsseldorf patient-bootstrap comparisons.

## Software environment

The analyses reported in the manuscript were run in R 4.5.3 on Apple silicon under macOS Tahoe 26.5.1.

Recorded package versions include:

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

The released scripts use deterministic seeds where stochastic operations affect reported results and retain audit information needed to trace the validation design.

## Verification outputs

### Adult analysis

Selected Adult Appendicitis verification outputs are stored in:

```text
results/adult/
```

These include aggregated held-out predictions, overall metrics, bootstrap comparisons, repeated-CV summaries, and statistical-analysis manifests.

### RPA analysis

Selected RPA verification outputs are stored in:

```text
results/rpa/
```

These include repeated held-out predictions, aggregated predictions, model-status audits, overall metrics, bootstrap comparisons, repeated-CV summaries, and statistical-analysis manifests.

### Düsseldorf locked external validation

Selected locked external-validation outputs are stored in:

```text
results/dusseldorf/
```

These include cohort summaries, model-status checks, preprocessing parameters, WBC harmonization audits, predictor-shift summaries, and external performance metrics.

### Fully nested RPA selector

Selected fully nested incremental-selector outputs are stored in:

```text
results/nested_selector/
```

These include patient-level nested predictions, dependency audits, selector coefficients, Risk@80% results, risk-coverage results, error-prediction metrics, bootstrap contrasts, WBC sensitivity outputs, asymmetric-cost analyses, and false-negative capture.

### Locked external selector

The reproducibility outputs for the frozen RPA selector applied to Düsseldorf are stored in:

```text
results/locked_external_selector/
```

The folder contains:

```text
01_final_rpa_selector_for_external_application.csv
02_dusseldorf_frozen_selector_predictions.csv
03_dusseldorf_selective_risk_primary_80pct.csv
04_dusseldorf_selective_risk_coverage.csv
05_dusseldorf_selective_risk_naurc.csv
06_dusseldorf_paired_bootstrap_rule_comparisons_80pct.csv
07_method_manifest.csv
08_manuscript_value_audit.csv
sessionInfo.txt
```

`08_manuscript_value_audit.csv` verifies that the reproduced Düsseldorf Risk@80 values and paired-bootstrap Table 4 results agree with the values reported in the manuscript at the stated reporting precision.

### Paired source-to-target experiment

Selected paired-shift outputs are stored in:

```text
results/paired_shift/
```

These include replication-level metrics, target manifests, model-fit status, within-replication contrasts, Monte Carlo confidence intervals, and publication-ready selective-risk summaries.

### Stationary simulation

Selected stationary-simulation outputs are stored in:

```text
results/stationary_simulation/
```

These contain simulation-unit model metrics, summary tables, and inferential/audit outputs. Large files in this directory are tracked using Git LFS.

## Interpretation

B-RISA is an evaluation framework rather than a new Bayesian classifier. Its main question is whether posterior probability standard deviation adds decision value beyond confidence already encoded in the predicted probability.

For BLR, BLDA, and BART in the evaluated appendicitis datasets, posterior probability SD did not demonstrate additional benefit over confidence margin for minimizing overall selective 0–1 error at a classification threshold of 0.50.

This result does not establish that Bayesian uncertainty is generally ineffective, nor does it establish that one deferral rule is clinically safer under every loss function, target population, or workflow.

The Düsseldorf rule comparisons are exploratory and post hoc. The asymmetric false-negative cost analyses use illustrative sensitivity parameters rather than clinically established utilities.

## Citation

If you use this code or the accompanying results, please cite the associated manuscript:

**Does Bayesian Posterior Uncertainty Improve Selective Clinical Prediction? Evidence from Geographic Validation of Appendicitis Models**

Full publication information will be added after publication.

## Contact

Kazi Nur Uddin  
Department of Mathematical Sciences  
Kent State University  
Email: kuddin@kent.edu