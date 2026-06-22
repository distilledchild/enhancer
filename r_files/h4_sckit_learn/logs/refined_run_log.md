## Refined Feature Run

Date: 2026-06-22

Input: `/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features.csv`
Refined feature table: `/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn/data/atac_enhancer_anchor_features_refined.csv`

Added features:

- `log_distance`
- `log_loop_span`
- `distance_per_loop_span`
- `anchor_width_sum`
- `anchor_width_diff`
- `anchor_width_ratio`
- `where_classification`

Best refined model by split:

| split                    | model              |   balanced_accuracy |   average_precision |
|:-------------------------|:-------------------|--------------------:|--------------------:|
| chromosome_group_holdout | LogisticRegression |            0.55831  |            0.920113 |
| stratified_80_20         | LogisticRegression |            0.553076 |            0.917478 |

Change from initial run:

| model              | split                    |   delta_balanced_accuracy |   delta_average_precision |
|:-------------------|:-------------------------|--------------------------:|--------------------------:|
| LogisticRegression | chromosome_group_holdout |                0.00817061 |               0.000471058 |
| RandomForest       | stratified_80_20         |                0.00351549 |              -0.000308868 |
| DummyMostFrequent  | stratified_80_20         |                0          |               0           |
| DummyMostFrequent  | chromosome_group_holdout |                0          |               0           |
| RandomForest       | chromosome_group_holdout |               -0.00133573 |              -0.000628968 |
| LogisticRegression | stratified_80_20         |               -0.00286255 |               0.00752839  |

Main caveat: These features are derived from existing loop geometry and annotations only. This remains exploratory and should not be framed as biological discovery.
