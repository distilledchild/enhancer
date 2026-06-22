#!/usr/bin/env python3
"""Validate outputs from the ATAC support scikit-learn workflow."""

from pathlib import Path

import pandas as pd


PROJECT_DIR = Path("/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn")
RESULTS_DIR = PROJECT_DIR / "results"

METRICS_PATH = RESULTS_DIR / "model_metrics.csv"
RF_IMPORTANCE_PATH = RESULTS_DIR / "feature_importance_random_forest.csv"
LOGISTIC_COEF_PATH = RESULTS_DIR / "logistic_regression_coefficients.csv"


def require_file(path: Path) -> None:
    if not path.exists():
        raise AssertionError(f"Missing expected output: {path}")


def main() -> None:
    for path in [METRICS_PATH, RF_IMPORTANCE_PATH, LOGISTIC_COEF_PATH]:
        require_file(path)

    metrics = pd.read_csv(METRICS_PATH)
    required_metric_columns = {
        "model",
        "split",
        "accuracy",
        "balanced_accuracy",
        "precision",
        "recall",
        "f1",
        "roc_auc",
        "average_precision",
        "tn",
        "fp",
        "fn",
        "tp",
    }
    missing_metric_columns = required_metric_columns.difference(metrics.columns)
    if missing_metric_columns:
        raise AssertionError(f"Missing metric columns: {sorted(missing_metric_columns)}")

    expected_models = {"DummyMostFrequent", "LogisticRegression", "RandomForest"}
    observed_models = set(metrics["model"])
    if not expected_models.issubset(observed_models):
        raise AssertionError(f"Missing models: {sorted(expected_models.difference(observed_models))}")

    rf_importance = pd.read_csv(RF_IMPORTANCE_PATH)
    if not {"feature", "importance"}.issubset(rf_importance.columns):
        raise AssertionError("Random forest importance file must include feature and importance columns")
    if rf_importance.empty:
        raise AssertionError("Random forest importance table is empty")

    logistic_coef = pd.read_csv(LOGISTIC_COEF_PATH)
    if not {"feature", "coefficient"}.issubset(logistic_coef.columns):
        raise AssertionError("Logistic coefficient file must include feature and coefficient columns")
    if logistic_coef.empty:
        raise AssertionError("Logistic coefficient table is empty")

    print(f"PASS metrics rows: {len(metrics)}")
    print(f"PASS random forest features: {len(rf_importance)}")
    print(f"PASS logistic features: {len(logistic_coef)}")


if __name__ == "__main__":
    main()
