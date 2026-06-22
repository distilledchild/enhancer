# Exploratory scikit-learn Workflow for ATAC-supported Enhancer Anchors

## Purpose

This folder contains an exploratory machine-learning workflow using regulatory genomics data from the Hi-C promoter-enhancer interaction project.

The goal is to practice and document a reproducible scikit-learn workflow that tests whether loop-derived regulatory features can model ATAC-seq support for Hi-C-derived enhancer anchors.

This is an applied machine-learning and workflow-development exercise. It should not be described as evidence for a strong predictive biological signal or as a manuscript-level biological discovery.

## Biological Framing

Question:

Can simple regulatory genomics features from Hi-C-derived promoter-enhancer interactions help classify whether an enhancer anchor overlaps ATAC-seq accessible chromatin?

Label:

- `1`: enhancer anchor overlaps an ATAC-seq peak
- `0`: enhancer anchor does not overlap an ATAC-seq peak

Main caveat:

The label is highly imbalanced. In the current feature table, 13,428 of 15,085 enhancer anchors are ATAC-supported. Therefore, raw accuracy is not very informative. Balanced accuracy, ROC-AUC, average precision, and comparison against a dummy baseline are more useful.

## Workflow Structure

```text
h4_sckit_learn/
  data/
    atac_enhancer_anchor_features.csv
    atac_enhancer_anchor_features_refined.csv
  logs/
    run_log.md
    refined_run_log.md
  results/
    model_metrics.csv
    feature_importance_random_forest.csv
    logistic_regression_coefficients.csv
    model_result_summary.md
    model_metric_comparison.png
    top_random_forest_features.png
    top_logistic_coefficients.png
    refined_model_metrics.csv
    refined_feature_importance_random_forest.csv
    refined_logistic_regression_coefficients.csv
    refined_vs_initial_model_metrics.csv
  scripts/
    01_export_atac_anchor_features.R
    02_train_atac_support_models.py
    03_summarize_atac_support_results.py
    04_train_refined_atac_support_models.py
  tests/
    test_export_atac_anchor_features.R
    test_train_atac_support_models.py
    test_summarize_atac_support_results.py
    test_train_refined_atac_support_models.py
```

## Step 01: Export Feature Table

Script:

```text
scripts/01_export_atac_anchor_features.R
```

Purpose:

Read the ATAC overlap table and the final Hi-C promoter-enhancer interaction RDS, then export one row per enhancer anchor.

Inputs:

```text
/Users/pete/Desktop/playground/enhancer/r_files/atac_validation/atac_loop_anchor_overlap_detail.csv
/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/r_files/rds/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds
```

Output:

```text
data/atac_enhancer_anchor_features.csv
```

Key operations:

- Filter ATAC overlap table to enhancer anchors.
- Convert ATAC overlap status to `label_atac_overlap`.
- Collapse loop annotations to one row per loop.
- Join loop-level regulatory features with the ATAC label.

Run:

```bash
Rscript scripts/01_export_atac_anchor_features.R
```

## Step 02: Train Baseline Models

Script:

```text
scripts/02_train_atac_support_models.py
```

Purpose:

Train baseline and interpretable scikit-learn classifiers.

Models:

- `DummyClassifier(strategy="most_frequent")`
- `LogisticRegression(class_weight="balanced")`
- `RandomForestClassifier(class_weight="balanced")`

Evaluation splits:

- Stratified 80/20 train/test split
- Chromosome group holdout split by `chr1`

Outputs:

```text
results/model_metrics.csv
results/feature_importance_random_forest.csv
results/logistic_regression_coefficients.csv
logs/run_log.md
```

Run:

```bash
python3 scripts/02_train_atac_support_models.py
```

## Step 03: Summarize Results

Script:

```text
scripts/03_summarize_atac_support_results.py
```

Purpose:

Read model outputs from Step 02 and create a compact interpretation summary and figures.

Outputs:

```text
results/model_result_summary.md
results/model_metric_comparison.png
results/top_random_forest_features.png
results/top_logistic_coefficients.png
```

Run:

```bash
python3 scripts/03_summarize_atac_support_results.py
```

## Step 04: Refined Feature Modeling

Script:

```text
scripts/04_train_refined_atac_support_models.py
```

Purpose:

Add leakage-free derived features from loop geometry and existing annotations, then compare refined model performance with the initial run.

Added features:

- `log_distance`
- `log_loop_span`
- `distance_per_loop_span`
- `anchor_width_sum`
- `anchor_width_diff`
- `anchor_width_ratio`
- `where_classification`

Outputs:

```text
data/atac_enhancer_anchor_features_refined.csv
results/refined_model_metrics.csv
results/refined_feature_importance_random_forest.csv
results/refined_logistic_regression_coefficients.csv
results/refined_vs_initial_model_metrics.csv
logs/refined_run_log.md
```

Run:

```bash
python3 scripts/04_train_refined_atac_support_models.py
```

## Current Results

Initial model run:

| Model | Split | Balanced Accuracy | Average Precision |
|---|---|---:|---:|
| DummyMostFrequent | Stratified 80/20 | 0.5000 | 0.8903 |
| LogisticRegression | Stratified 80/20 | 0.5559 | 0.9099 |
| RandomForest | Stratified 80/20 | 0.4984 | 0.9006 |
| DummyMostFrequent | Chromosome holdout | 0.5000 | 0.8975 |
| LogisticRegression | Chromosome holdout | 0.5501 | 0.9196 |
| RandomForest | Chromosome holdout | 0.5049 | 0.9087 |

Refined feature run:

| Model | Split | Balanced Accuracy | Average Precision |
|---|---|---:|---:|
| DummyMostFrequent | Stratified 80/20 | 0.5000 | 0.8903 |
| LogisticRegression | Stratified 80/20 | 0.5531 | 0.9175 |
| RandomForest | Stratified 80/20 | 0.5020 | 0.9003 |
| DummyMostFrequent | Chromosome holdout | 0.5000 | 0.8975 |
| LogisticRegression | Chromosome holdout | 0.5583 | 0.9201 |
| RandomForest | Chromosome holdout | 0.5036 | 0.9081 |

Conservative interpretation:

Logistic regression showed a modest improvement over the dummy baseline. Random forest performed close to baseline in this first pass. The refined feature set slightly improved chromosome holdout balanced accuracy for logistic regression, but the signal remains modest.

## How To Validate

Run all validation checks:

```bash
Rscript tests/test_export_atac_anchor_features.R
python3 tests/test_train_atac_support_models.py
python3 tests/test_summarize_atac_support_results.py
python3 tests/test_train_refined_atac_support_models.py
```

Expected:

All four checks should pass.

## Safe Resume Or Interview Wording

Safe wording:

```text
Built an exploratory scikit-learn workflow using regulatory genomics data to model ATAC-seq support for Hi-C-derived enhancer anchors, including feature engineering, dummy baseline comparison, logistic regression, random forest, chromosome holdout evaluation, balanced metrics, and feature-importance interpretation.
```

Avoid saying:

```text
Developed deep-learning models for regulatory genomics.
Discovered predictive regulatory mechanisms.
Built a strong production-grade biological predictor.
```

The correct framing is:

```text
exploratory scikit-learn workflow
applied machine-learning practice
interpretable baseline modeling
regulatory genomics feature engineering
balanced evaluation under label imbalance
```
