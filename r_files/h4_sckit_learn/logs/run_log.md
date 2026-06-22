## Run 1

Date: 2026-06-21

Goal: First exploratory classifier for ATAC-supported enhancer anchors.

Input: `/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv`

Label: `label_atac_overlap`

Features: `resolution, distance, anchor1_width, anchor2_width, loop_span, n_assigned_genes, n_components, has_tss_component, has_promoter_component, WHERE, classification`

Models: DummyMostFrequent, LogisticRegression, RandomForest

Splits: stratified 80/20 and chromosome group holdout by chr1

Main caveat: The ATAC overlap label is imbalanced, so raw accuracy is not sufficient.

Best model by split:

| split                    | model              |   balanced_accuracy |   average_precision |
|:-------------------------|:-------------------|--------------------:|--------------------:|
| chromosome_group_holdout | LogisticRegression |            0.550139 |            0.919642 |
| stratified_80_20         | LogisticRegression |            0.555938 |            0.90995  |
