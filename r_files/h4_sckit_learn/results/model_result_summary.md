# ATAC-supported enhancer anchors: model result summary

Date: 2026-06-22

## Input

Feature table: `/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv`

Rows: 15,085
Positive labels: 13,428
Negative labels: 1,657
Positive label rate: 0.8902

The label is imbalanced, so raw accuracy is not the main metric. Balanced accuracy, ROC-AUC, and average precision are more useful for judging whether the model improves over the dummy baseline.

## Model Metrics

| model              | split                    |   n_test |   positive_rate_test |   accuracy |   balanced_accuracy |   precision |   recall |     f1 |   roc_auc |   average_precision |   tn |   fp |   fn |   tp |
|:-------------------|:-------------------------|---------:|---------------------:|-----------:|--------------------:|------------:|---------:|-------:|----------:|--------------------:|-----:|-----:|-----:|-----:|
| DummyMostFrequent  | stratified_80_20         |     3017 |               0.8903 |     0.8903 |              0.5    |      0.8903 |   1      | 0.942  |    0.5    |              0.8903 |    0 |  331 |    0 | 2686 |
| LogisticRegression | stratified_80_20         |     3017 |               0.8903 |     0.5277 |              0.5559 |      0.9118 |   0.5197 | 0.6621 |    0.5776 |              0.9099 |  196 |  135 | 1290 | 1396 |
| RandomForest       | stratified_80_20         |     3017 |               0.8903 |     0.8545 |              0.4984 |      0.89   |   0.9546 | 0.9211 |    0.5329 |              0.9006 |   14 |  317 |  122 | 2564 |
| DummyMostFrequent  | chromosome_group_holdout |     5553 |               0.8975 |     0.8975 |              0.5    |      0.8975 |   1      | 0.946  |    0.5    |              0.8975 |    0 |  569 |    0 | 4984 |
| LogisticRegression | chromosome_group_holdout |     5553 |               0.8975 |     0.5278 |              0.5501 |      0.9156 |   0.5221 | 0.665  |    0.5719 |              0.9196 |  329 |  240 | 2382 | 2602 |
| RandomForest       | chromosome_group_holdout |     5553 |               0.8975 |     0.8658 |              0.5049 |      0.8985 |   0.9589 | 0.9277 |    0.5344 |              0.9087 |   29 |  540 |  205 | 4779 |

## Best Model By Split

| split                    | model              |   balanced_accuracy |   average_precision |
|:-------------------------|:-------------------|--------------------:|--------------------:|
| chromosome_group_holdout | LogisticRegression |              0.5501 |              0.9196 |
| stratified_80_20         | LogisticRegression |              0.5559 |              0.9099 |

## Feature Interpretation

Top random forest features:

| feature                |   importance |
|:-----------------------|-------------:|
| distance               |   0.622795   |
| loop_span              |   0.328083   |
| resolution             |   0.0064552  |
| anchor1_width          |   0.00642109 |
| anchor2_width          |   0.00591098 |
| has_tss_component      |   0.00582979 |
| has_promoter_component |   0.00540426 |
| WHERE_UP               |   0.00522969 |

Top logistic regression coefficients:

| feature                |   coefficient |   abs_coefficient |
|:-----------------------|--------------:|------------------:|
| WHERE_UP               |        0.1007 |            0.1007 |
| resolution             |        0.0867 |            0.0867 |
| anchor1_width          |        0.0867 |            0.0867 |
| anchor2_width          |        0.0867 |            0.0867 |
| loop_span              |       -0.0843 |            0.0843 |
| WHERE_DOWN             |       -0.0758 |            0.0758 |
| classification_Both_OK |        0.062  |            0.062  |
| distance               |        0.0612 |            0.0612 |

## Conservative Interpretation

LogisticRegression showed a modest improvement over DummyMostFrequent in balanced accuracy and average precision. RandomForest performed close to the baseline in this first pass.

This should be described as an exploratory scikit-learn workflow rather than evidence for a strong predictive biological signal.

Safe wording:

> Built an exploratory scikit-learn workflow using regulatory genomics features to model ATAC-seq support for Hi-C-derived enhancer anchors, including dummy baseline comparison, logistic regression, random forest, balanced metrics, chromosome holdout evaluation, and feature-importance interpretation.
