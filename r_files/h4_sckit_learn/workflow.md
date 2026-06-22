# scikit-learn Workflow: Hi-C Regulatory Genomics Projects

## Purpose

This folder records exploratory scikit-learn analyses based on the Hi-C promoter-enhancer interaction project. The goal is to gain hands-on machine-learning experience using existing regulatory genomics datasets while staying biologically interpretable and closely connected to the PE-interaction research.

The analyses below are intended as training and exploratory research, not as primary manuscript claims unless they are later validated carefully.

## Available Data Sources

Potential inputs from the current project include:

- Hi-C loop calls across HRDP rat strains
- Loop anchor coordinates and anchor midpoint information
- Promoter/TSS assignment results
- Enhancer-anchor and promoter-anchor classifications
- Final CTCF-supported regulatory interaction set
- Loop length and anchor-to-feature distance
- CTCF motif presence, orientation, count, or motif score
- ATAC-seq overlap at candidate enhancer anchors
- AlphaGenome enhancer prediction scores, if available
- Gene-level summaries from Hi-C-based variant-to-gene or H-MAGMA follow-up analysis

## Project 1: Enhancer Anchor Activity Prediction

### Biological Question

Can we predict whether a Hi-C-derived enhancer anchor is supported by open chromatin using genomic and regulatory features?

### Machine-Learning Framing

Binary classification.

### Possible Label

- `1`: enhancer anchor overlaps ATAC-seq accessible chromatin
- `0`: matched background anchor or non-ATAC-supported enhancer anchor

A random-shift or matched negative set should be used carefully to avoid making the classification task artificially easy.

### Candidate Features

- Loop length
- Anchor size
- Distance from enhancer anchor to assigned promoter/TSS anchor
- Distance from anchor midpoint to nearest TSS/promoter
- Number of loops connected to the anchor
- CTCF motif count or score near anchor
- CTCF motif orientation category
- AlphaGenome enhancer score, if available
- Chromosome or genomic context features, if useful

### Candidate Models

- LogisticRegression
- RandomForestClassifier
- GradientBoostingClassifier
- HistGradientBoostingClassifier

### Evaluation

- Train/test split
- Cross-validation
- Precision, recall, F1, ROC-AUC, PR-AUC
- Feature importance or permutation importance
- Chromosome-based holdout split if possible

### Main Caution

Random genomic splitting can overestimate performance because nearby genomic regions are not independent. A chromosome-based split is more conservative and more convincing.

### Implementation Note

First exploratory implementation completed using ATAC overlap labels, loop-derived features, DummyClassifier, LogisticRegression, RandomForestClassifier, stratified train/test split, and chromosome group holdout split. Results are recorded in `results/model_metrics.csv`, `results/feature_importance_random_forest.csv`, `results/logistic_regression_coefficients.csv`, and `logs/run_log.md`.

Result summarization was added in `scripts/03_summarize_atac_support_results.py`. Interpretation outputs are recorded in `results/model_result_summary.md`, `results/model_metric_comparison.png`, `results/top_random_forest_features.png`, and `results/top_logistic_coefficients.png`.

Refined feature modeling was added in `scripts/04_train_refined_atac_support_models.py`. This run adds leakage-free derived features from loop geometry and existing annotations, including `log_distance`, `log_loop_span`, `distance_per_loop_span`, `anchor_width_sum`, `anchor_width_diff`, `anchor_width_ratio`, and `where_classification`. Outputs are recorded in `data/atac_enhancer_anchor_features_refined.csv`, `results/refined_model_metrics.csv`, `results/refined_feature_importance_random_forest.csv`, `results/refined_logistic_regression_coefficients.csv`, `results/refined_vs_initial_model_metrics.csv`, and `logs/refined_run_log.md`.

## Project 2: Regulatory Gene Clustering

### Biological Question

Do genes with similar Hi-C regulatory-loop profiles cluster into biologically meaningful groups?

### Machine-Learning Framing

Unsupervised learning.

### Unit Of Analysis

Gene-level summary table.

### Candidate Features

- Number of assigned enhancer-promoter loops per gene
- Number of CTCF-supported loops per gene
- Number of ATAC-supported enhancer anchors per gene
- Median loop length per gene
- Minimum/median anchor-to-TSS distance
- Maximum or mean AlphaGenome enhancer score, if available
- Whether the gene appears in GWAS/H-MAGMA candidate lists

### Candidate Methods

- StandardScaler
- PCA
- KMeans
- AgglomerativeClustering
- GaussianMixture

### Interpretation

- Compare clusters by regulatory complexity
- Perform GO or pathway enrichment per cluster
- Check whether development-related genes or highly connected genes concentrate in specific clusters

### Main Caution

Clustering does not produce a true label. Results should be described as exploratory grouping, not discovery of definitive gene classes.

## Project 3: CTCF-Supported Loop Classification

### Biological Question

Which non-CTCF-derived loop features distinguish CTCF-supported regulatory loops from non-CTCF-supported loops?

### Machine-Learning Framing

Binary classification.

### Possible Label

- `1`: final CTCF-supported loop or interaction
- `0`: loop not supported by the CTCF filtering criteria

### Candidate Features

- Loop length
- Anchor-to-TSS/promoter distance
- Anchor type combination
- Number of assigned genes
- ATAC support at enhancer anchor
- Direction of assigned TSS/promoter relative to anchor midpoint
- Chromosome-level or gene-level summaries

### Important Leakage Warning

Do not use CTCF motif count, CTCF motif score, or CTCF motif orientation as predictors if the label itself is defined by CTCF support. That would create label leakage and make the model trivially predictive.

CTCF-related features may be used only in a separate demonstration model where the leakage issue is explicitly stated.

## Project 4: GWAS Candidate Gene Prioritization

### Biological Question

Can regulatory features from Hi-C loops help prioritize candidate genes from GWAS-linked promoter-enhancer interactions?

### Machine-Learning Framing

Exploratory ranking or weakly supervised classification.

### Possible Unit Of Analysis

Gene-level table from H-MAGMA or related gene-based association results.

### Candidate Features

- H-MAGMA gene p-value or rank
- Number of enhancer loops assigned to gene
- Number of CTCF-supported interactions assigned to gene
- ATAC-supported enhancer anchor count
- Mean or max AlphaGenome enhancer score, if available
- Median loop length
- Minimum anchor-to-TSS/promoter distance
- Functional enrichment category or pathway label, if used only for interpretation

### Candidate Methods

- LogisticRegression or RandomForestClassifier for top-candidate vs non-top-candidate comparison
- RandomForestRegressor or GradientBoostingRegressor for ranking-score prediction
- PCA or clustering for exploratory gene-prioritization patterns

### Main Caution

There is no gold-standard causal gene label in the current dataset. This project should be described as candidate prioritization or exploratory modeling, not causal-gene prediction.

## Recommended Starting Order

1. Build a clean loop-anchor or gene-level feature table.
2. Start with Project 1: enhancer anchor activity prediction.
3. Add Project 2: regulatory gene clustering.
4. Try Project 3 only after defining leakage-free predictors.
5. Try Project 4 after H-MAGMA/GWAS gene-level tables are organized.

## Practical Implementation Plan

1. Create one R script or Python notebook to export a clean CSV feature table from existing R objects.
2. Use Python/scikit-learn for modeling.
3. Keep all model inputs in a `data/` subfolder.
4. Save figures and metrics in a `results/` subfolder.
5. Record each model run, feature set, label definition, and evaluation split in this workflow file or a separate log.
