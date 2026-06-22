# scikit-learn Workflow Summary: First Exploratory Model

## Goal

Build one small, reproducible, interpretable scikit-learn workflow using the Hi-C promoter-enhancer interaction project.

The first model will test whether simple regulatory-genomics features from Hi-C-derived enhancer-promoter interactions can classify whether an enhancer anchor overlaps ATAC-seq accessible chromatin.

This is an exploratory machine-learning practice workflow, not a manuscript claim.

## Project To Implement First

**Project 1: Enhancer Anchor Activity Prediction**

Biological question:

Can we predict whether a Hi-C-derived enhancer anchor is supported by open chromatin using genomic and regulatory features?

Machine-learning framing:

- Binary classification
- Label:
  - `1`: enhancer anchor overlaps ATAC-seq peak
  - `0`: enhancer anchor does not overlap ATAC-seq peak

Important caution:

The current ATAC overlap rate is high. In the existing ATAC validation output, enhancer anchors show about 89% ATAC overlap. This means a naive model may look accurate even if it is not biologically useful. Therefore, use balanced metrics and compare against a dummy baseline.

## Working Directory

Main workflow directory:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn
```

Existing planning document:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/workflow.md
```

Existing Korean translation:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/workflow-kr.md
```

## Input Files

ATAC overlap detail table:

```text
/Users/pete/Desktop/playground/enhancer/r_files/atac_validation/atac_loop_anchor_overlap_detail.csv
```

Current structure:

```text
loop_id, anchor_type, category, atac_overlap
```

Main Hi-C promoter-enhancer interaction RDS:

```text
/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/r_files/rds/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds
```

Useful columns confirmed in the RDS:

```text
distance
loop.id
chr1
x1
x2
chr2
y1
y2
resolution
gene_id
gene_name
component
WHERE
classification
```

## Output Files To Create

Create these folders:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/logs
```

Create these files:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/01_export_atac_anchor_features.R
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/02_train_atac_support_models.py
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/model_metrics.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/feature_importance_random_forest.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/logistic_regression_coefficients.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/logs/run_log.md
```

## Step 1: Create Folder Structure

Run:

```bash
mkdir -p \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts" \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data" \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results" \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/logs"
```

Expected:

The four folders exist under:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn
```

## Step 2: Export A Clean Feature Table

Create:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/01_export_atac_anchor_features.R
```

Purpose:

Read the existing ATAC overlap table and the Hi-C promoter-enhancer interaction RDS, then export one row per enhancer anchor.

Input:

```text
/Users/pete/Desktop/playground/enhancer/r_files/atac_validation/atac_loop_anchor_overlap_detail.csv
/Users/pete/Library/CloudStorage/Dropbox-UTHSCGGI/K P/Gateway_to_Hao/enhancer/r_files/rds/df_final_up_down_directional_point_decision_COMBINED_OK_filtered_lt_Q3_final_200kb.rds
```

Output:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv
```

Feature table design:

- Filter ATAC table to `anchor_type == "enhancer"`.
- Convert `atac_overlap` to integer label:
  - `TRUE` -> `1`
  - `FALSE` -> `0`
- Join loop-level features from the RDS by `loop_id == loop.id`.
- Collapse RDS rows to one row per loop before joining, because one loop may be assigned to more than one gene.
- Keep only non-ATAC features as predictors.

Candidate columns:

```text
loop_id
label_atac_overlap
chr1
chr2
resolution
distance
anchor1_width
anchor2_width
loop_span
n_assigned_genes
n_components
has_tss_component
has_promoter_component
WHERE
classification
```

Do not include as predictors:

```text
atac_overlap
anchor_type
```

Leakage rule:

Do not use any feature directly derived from the ATAC overlap label.

Run:

```bash
Rscript "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/01_export_atac_anchor_features.R"
```

Expected:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv
```

The exported CSV should have about 15,085 rows before removing incomplete rows.

## Step 3: Train Baseline And Interpretable Models

Create:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/02_train_atac_support_models.py
```

Purpose:

Train simple baseline and interpretable scikit-learn classifiers to predict `label_atac_overlap`.

Input:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv
```

Models:

```text
DummyClassifier(strategy="most_frequent")
LogisticRegression(class_weight="balanced", max_iter=5000)
RandomForestClassifier(n_estimators=300, class_weight="balanced", random_state=42)
```

Primary split:

```text
Stratified train/test split, 80/20, random_state=42
```

Secondary split:

```text
Group-aware split by chromosome using chr1 as group
```

Recommended features for the first pass:

Numeric:

```text
resolution
distance
anchor1_width
anchor2_width
loop_span
n_assigned_genes
n_components
has_tss_component
has_promoter_component
```

Categorical:

```text
WHERE
classification
```

Avoid using chromosome as a predictor in the first pass. Chromosome can be used for group-aware splitting, but using it as a predictor may encourage genomic-location shortcuts.

Metrics:

```text
accuracy
balanced_accuracy
precision
recall
f1
roc_auc
average_precision
confusion_matrix
```

Outputs:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/model_metrics.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/feature_importance_random_forest.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/logistic_regression_coefficients.csv
```

Run:

```bash
python3 "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/02_train_atac_support_models.py"
```

Expected:

The script prints model metrics and writes result files under:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results
```

## Step 4: Interpret Results Conservatively

Interpretation checklist:

- Compare LogisticRegression and RandomForest against DummyClassifier.
- Focus on `balanced_accuracy`, `f1`, `roc_auc`, and `average_precision`, not raw accuracy alone.
- Check whether the model actually improves over the dummy baseline.
- Inspect feature importance and coefficients.
- If performance is weak, describe it as a useful hands-on workflow rather than a biological discovery.
- If performance is strong, still treat it cautiously because labels are imbalanced and nearby genomic regions may be correlated.

## Step 5: Record The Run

Create or update:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/logs/run_log.md
```

Each run should record:

```text
date
input feature table
label definition
feature columns
model list
split strategy
metrics summary
main caveat
next action
```

Example entry:

```markdown
## Run 1

Date: 2026-06-21

Goal: First exploratory classifier for ATAC-supported enhancer anchors.

Input:
`/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv`

Label:
`label_atac_overlap`

Models:
DummyClassifier, LogisticRegression, RandomForestClassifier

Main caveat:
The ATAC overlap label is imbalanced, so raw accuracy is not sufficient.
```

## Step 6: Update Main Workflow Document

After the first successful run, update:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/workflow.md
```

Add a short implementation note under Project 1:

```text
First exploratory implementation completed using ATAC overlap labels, loop-derived features, LogisticRegression, RandomForestClassifier, and baseline comparison.
```

Do not add a resume bullet until the scripts and result files actually exist.

## Acceptance Criteria For The First Model

The first scikit-learn workflow is complete when all of the following exist:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/01_export_atac_anchor_features.R
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/02_train_atac_support_models.py
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/model_metrics.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/feature_importance_random_forest.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results/logistic_regression_coefficients.csv
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/logs/run_log.md
```

Minimum successful result:

- Feature table created without missing required columns.
- Dummy baseline trained.
- Logistic regression trained.
- Random forest trained.
- Metrics saved.
- At least one feature-importance or coefficient table saved.

## Important Language For Applications Or Interviews

Safe phrasing after implementation:

```text
I built an exploratory scikit-learn workflow using regulatory genomics data to model ATAC-seq support for Hi-C-derived enhancer anchors with interpretable baseline classifiers, including logistic regression and random forest.
```

Avoid saying:

```text
I developed deep-learning models for genomics.
```

Avoid saying:

```text
I discovered predictive regulatory mechanisms.
```

The correct framing is:

```text
hands-on applied machine-learning workflow
exploratory regulatory genomics modeling
interpretable baseline models
reproducible feature table and metrics
```

## Immediate Next Command

Start implementation by creating folders:

```bash
mkdir -p \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts" \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data" \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/results" \
  "/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/logs"
```

Then create:

```text
/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/scripts/01_export_atac_anchor_features.R
```
