#!/usr/bin/env python3
"""Train baseline scikit-learn models for ATAC-supported enhancer anchors.

This script uses the feature table exported by 01_export_atac_anchor_features.R.
It trains a dummy baseline, logistic regression, and random forest classifier.
Outputs are intended for exploratory machine-learning practice, not biological
claims.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.dummy import DummyClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GroupShuffleSplit, train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler


PROJECT_DIR = Path("/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn")
DATA_PATH = PROJECT_DIR / "data" / "atac_enhancer_anchor_features.csv"
RESULTS_DIR = PROJECT_DIR / "results"
LOGS_DIR = PROJECT_DIR / "logs"

METRICS_PATH = RESULTS_DIR / "model_metrics.csv"
RF_IMPORTANCE_PATH = RESULTS_DIR / "feature_importance_random_forest.csv"
LOGISTIC_COEF_PATH = RESULTS_DIR / "logistic_regression_coefficients.csv"
RUN_LOG_PATH = LOGS_DIR / "run_log.md"

LABEL_COLUMN = "label_atac_overlap"
GROUP_COLUMN = "chr1"

NUMERIC_FEATURES = [
    "resolution",
    "distance",
    "anchor1_width",
    "anchor2_width",
    "loop_span",
    "n_assigned_genes",
    "n_components",
    "has_tss_component",
    "has_promoter_component",
]

CATEGORICAL_FEATURES = [
    "WHERE",
    "classification",
]

FEATURE_COLUMNS = NUMERIC_FEATURES + CATEGORICAL_FEATURES


def build_preprocessor() -> ColumnTransformer:
    """Create preprocessing steps for numeric and categorical features."""
    return ColumnTransformer(
        transformers=[
            ("numeric", StandardScaler(), NUMERIC_FEATURES),
            (
                "categorical",
                OneHotEncoder(handle_unknown="ignore", sparse_output=False),
                CATEGORICAL_FEATURES,
            ),
        ]
    )


def build_models() -> dict[str, Pipeline]:
    """Create baseline and interpretable scikit-learn pipelines."""
    return {
        "DummyMostFrequent": Pipeline(
            steps=[
                ("preprocessor", build_preprocessor()),
                ("classifier", DummyClassifier(strategy="most_frequent")),
            ]
        ),
        "LogisticRegression": Pipeline(
            steps=[
                ("preprocessor", build_preprocessor()),
                (
                    "classifier",
                    LogisticRegression(class_weight="balanced", max_iter=5000, random_state=42),
                ),
            ]
        ),
        "RandomForest": Pipeline(
            steps=[
                ("preprocessor", build_preprocessor()),
                (
                    "classifier",
                    RandomForestClassifier(
                        n_estimators=300,
                        class_weight="balanced",
                        random_state=42,
                        n_jobs=-1,
                    ),
                ),
            ]
        ),
    }


def safe_score(metric_func: Any, y_true: pd.Series, y_score: np.ndarray) -> float:
    """Return a metric value or NA when the metric is undefined."""
    try:
        return float(metric_func(y_true, y_score))
    except ValueError:
        return float("nan")


def evaluate_model(model_name: str, split_name: str, model: Pipeline, x_test: pd.DataFrame, y_test: pd.Series) -> dict[str, Any]:
    """Evaluate a fitted model using balanced classification metrics."""
    y_pred = model.predict(x_test)
    y_score = model.predict_proba(x_test)[:, 1]
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()

    return {
        "model": model_name,
        "split": split_name,
        "n_test": int(len(y_test)),
        "positive_rate_test": float(np.mean(y_test)),
        "accuracy": float(accuracy_score(y_test, y_pred)),
        "balanced_accuracy": float(balanced_accuracy_score(y_test, y_pred)),
        "precision": float(precision_score(y_test, y_pred, zero_division=0)),
        "recall": float(recall_score(y_test, y_pred, zero_division=0)),
        "f1": float(f1_score(y_test, y_pred, zero_division=0)),
        "roc_auc": safe_score(roc_auc_score, y_test, y_score),
        "average_precision": safe_score(average_precision_score, y_test, y_score),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def make_splits(df_features: pd.DataFrame) -> dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]]:
    """Create stratified and chromosome-grouped train/test splits."""
    x_all = df_features[FEATURE_COLUMNS]
    y_all = df_features[LABEL_COLUMN]

    x_train, x_test, y_train, y_test = train_test_split(
        x_all,
        y_all,
        test_size=0.2,
        random_state=42,
        stratify=y_all,
    )

    group_splitter = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=42)
    group_train_idx, group_test_idx = next(group_splitter.split(x_all, y_all, groups=df_features[GROUP_COLUMN]))

    return {
        "stratified_80_20": (x_train, x_test, y_train, y_test),
        "chromosome_group_holdout": (
            x_all.iloc[group_train_idx],
            x_all.iloc[group_test_idx],
            y_all.iloc[group_train_idx],
            y_all.iloc[group_test_idx],
        ),
    }


def extract_feature_names(model: Pipeline) -> list[str]:
    """Return transformed feature names from the fitted preprocessor."""
    preprocessor = model.named_steps["preprocessor"]
    return list(preprocessor.get_feature_names_out())


def write_run_log(metrics: pd.DataFrame) -> None:
    """Write a compact markdown log for the model run."""
    best_rows = metrics.sort_values(["split", "balanced_accuracy"], ascending=[True, False])
    best_summary = best_rows.groupby("split").head(1)[["split", "model", "balanced_accuracy", "average_precision"]]

    lines = [
        f"## Run 1",
        "",
        f"Date: {date.today().isoformat()}",
        "",
        "Goal: First exploratory classifier for ATAC-supported enhancer anchors.",
        "",
        f"Input: `{DATA_PATH}`",
        "",
        "Label: `label_atac_overlap`",
        "",
        f"Features: `{', '.join(FEATURE_COLUMNS)}`",
        "",
        "Models: DummyMostFrequent, LogisticRegression, RandomForest",
        "",
        "Splits: stratified 80/20 and chromosome group holdout by chr1",
        "",
        "Main caveat: The ATAC overlap label is imbalanced, so raw accuracy is not sufficient.",
        "",
        "Best model by split:",
        "",
        best_summary.to_markdown(index=False),
        "",
    ]

    RUN_LOG_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run the ATAC support classification workflow."""
    if not DATA_PATH.exists():
        raise FileNotFoundError(f"Missing input feature table: {DATA_PATH}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    df_features = pd.read_csv(DATA_PATH)
    required_columns = set(FEATURE_COLUMNS + [LABEL_COLUMN, GROUP_COLUMN])
    missing_columns = required_columns.difference(df_features.columns)
    if missing_columns:
        raise ValueError(f"Missing required columns: {sorted(missing_columns)}")

    df_features = df_features.dropna(subset=FEATURE_COLUMNS + [LABEL_COLUMN, GROUP_COLUMN]).copy()
    df_features[LABEL_COLUMN] = df_features[LABEL_COLUMN].astype(int)

    metrics_rows: list[dict[str, Any]] = []
    fitted_stratified_models: dict[str, Pipeline] = {}

    for split_name, (x_train, x_test, y_train, y_test) in make_splits(df_features).items():
        for model_name, model in build_models().items():
            model.fit(x_train, y_train)
            metrics_rows.append(evaluate_model(model_name, split_name, model, x_test, y_test))

            if split_name == "stratified_80_20":
                fitted_stratified_models[model_name] = model

    df_metrics = pd.DataFrame(metrics_rows)
    df_metrics.to_csv(METRICS_PATH, index=False)

    rf_model = fitted_stratified_models["RandomForest"]
    rf_features = extract_feature_names(rf_model)
    rf_importance = rf_model.named_steps["classifier"].feature_importances_
    df_rf_importance = pd.DataFrame({"feature": rf_features, "importance": rf_importance})
    df_rf_importance = df_rf_importance.sort_values("importance", ascending=False)
    df_rf_importance.to_csv(RF_IMPORTANCE_PATH, index=False)

    logistic_model = fitted_stratified_models["LogisticRegression"]
    logistic_features = extract_feature_names(logistic_model)
    logistic_coef = logistic_model.named_steps["classifier"].coef_[0]
    df_logistic_coef = pd.DataFrame({"feature": logistic_features, "coefficient": logistic_coef})
    df_logistic_coef = df_logistic_coef.assign(abs_coefficient=lambda df: df["coefficient"].abs())
    df_logistic_coef = df_logistic_coef.sort_values("abs_coefficient", ascending=False)
    df_logistic_coef.to_csv(LOGISTIC_COEF_PATH, index=False)

    write_run_log(df_metrics)

    print(f"Wrote metrics: {METRICS_PATH}")
    print(df_metrics[["model", "split", "balanced_accuracy", "average_precision"]].to_string(index=False))
    print(f"Wrote random forest importance: {RF_IMPORTANCE_PATH}")
    print(f"Wrote logistic coefficients: {LOGISTIC_COEF_PATH}")
    print(f"Wrote run log: {RUN_LOG_PATH}")


if __name__ == "__main__":
    main()
