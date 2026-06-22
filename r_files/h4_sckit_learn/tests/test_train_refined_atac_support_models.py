#!/usr/bin/env python3
"""Validate refined feature and model outputs for ATAC support workflow."""

from pathlib import Path

import pandas as pd


PROJECT_DIR = Path("/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn")
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
LOGS_DIR = PROJECT_DIR / "logs"

REFINED_FEATURE_PATH = DATA_DIR / "atac_enhancer_anchor_features_refined.csv"
REFINED_METRICS_PATH = RESULTS_DIR / "refined_model_metrics.csv"
REFINED_RF_IMPORTANCE_PATH = RESULTS_DIR / "refined_feature_importance_random_forest.csv"
REFINED_LOGISTIC_COEF_PATH = RESULTS_DIR / "refined_logistic_regression_coefficients.csv"
REFINED_COMPARISON_PATH = RESULTS_DIR / "refined_vs_initial_model_metrics.csv"
REFINED_LOG_PATH = LOGS_DIR / "refined_run_log.md"


def require_nonempty_file(path: Path) -> None:
    if not path.exists():
        raise AssertionError(f"Missing expected output: {path}")
    if path.stat().st_size == 0:
        raise AssertionError(f"Output file is empty: {path}")


def main() -> None:
    for path in [
        REFINED_FEATURE_PATH,
        REFINED_METRICS_PATH,
        REFINED_RF_IMPORTANCE_PATH,
        REFINED_LOGISTIC_COEF_PATH,
        REFINED_COMPARISON_PATH,
        REFINED_LOG_PATH,
    ]:
        require_nonempty_file(path)

    refined_features = pd.read_csv(REFINED_FEATURE_PATH)
    required_refined_columns = {
        "log_distance",
        "log_loop_span",
        "distance_per_loop_span",
        "anchor_width_sum",
        "anchor_width_diff",
        "anchor_width_ratio",
        "where_classification",
    }
    missing_refined_columns = required_refined_columns.difference(refined_features.columns)
    if missing_refined_columns:
        raise AssertionError(f"Missing refined feature columns: {sorted(missing_refined_columns)}")
    if len(refined_features) != 15085:
        raise AssertionError(f"Unexpected refined feature row count: {len(refined_features)}")

    refined_metrics = pd.read_csv(REFINED_METRICS_PATH)
    expected_models = {"DummyMostFrequent", "LogisticRegression", "RandomForest"}
    observed_models = set(refined_metrics["model"])
    if not expected_models.issubset(observed_models):
        raise AssertionError(f"Missing refined models: {sorted(expected_models.difference(observed_models))}")

    comparison = pd.read_csv(REFINED_COMPARISON_PATH)
    required_comparison_columns = {
        "model",
        "split",
        "initial_balanced_accuracy",
        "refined_balanced_accuracy",
        "delta_balanced_accuracy",
        "initial_average_precision",
        "refined_average_precision",
        "delta_average_precision",
    }
    missing_comparison_columns = required_comparison_columns.difference(comparison.columns)
    if missing_comparison_columns:
        raise AssertionError(f"Missing comparison columns: {sorted(missing_comparison_columns)}")

    print(f"PASS refined feature rows: {len(refined_features)}")
    print(f"PASS refined metrics rows: {len(refined_metrics)}")
    print(f"PASS comparison rows: {len(comparison)}")


if __name__ == "__main__":
    main()
