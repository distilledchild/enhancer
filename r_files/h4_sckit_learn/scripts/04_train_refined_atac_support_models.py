#!/usr/bin/env python3
"""Train refined scikit-learn models for ATAC-supported enhancer anchors.

This script keeps the first model workflow conservative and reproducible. It
adds leakage-free derived features from existing loop geometry and categorical
annotations, then compares refined model metrics against the initial run.
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
DATA_DIR = PROJECT_DIR / "data"
RESULTS_DIR = PROJECT_DIR / "results"
LOGS_DIR = PROJECT_DIR / "logs"

INPUT_FEATURE_PATH = DATA_DIR / "atac_enhancer_anchor_features.csv"
REFINED_FEATURE_PATH = DATA_DIR / "atac_enhancer_anchor_features_refined.csv"
INITIAL_METRICS_PATH = RESULTS_DIR / "model_metrics.csv"

REFINED_METRICS_PATH = RESULTS_DIR / "refined_model_metrics.csv"
REFINED_RF_IMPORTANCE_PATH = RESULTS_DIR / "refined_feature_importance_random_forest.csv"
REFINED_LOGISTIC_COEF_PATH = RESULTS_DIR / "refined_logistic_regression_coefficients.csv"
REFINED_COMPARISON_PATH = RESULTS_DIR / "refined_vs_initial_model_metrics.csv"
REFINED_LOG_PATH = LOGS_DIR / "refined_run_log.md"

LABEL_COLUMN = "label_atac_overlap"
GROUP_COLUMN = "chr1"

BASE_NUMERIC_FEATURES = [
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

REFINED_NUMERIC_FEATURES = [
    "log_distance",
    "log_loop_span",
    "distance_per_loop_span",
    "anchor_width_sum",
    "anchor_width_diff",
    "anchor_width_ratio",
]

NUMERIC_FEATURES = BASE_NUMERIC_FEATURES + REFINED_NUMERIC_FEATURES

CATEGORICAL_FEATURES = [
    "WHERE",
    "classification",
    "where_classification",
]

FEATURE_COLUMNS = NUMERIC_FEATURES + CATEGORICAL_FEATURES


def add_refined_features(df_features: pd.DataFrame) -> pd.DataFrame:
    """Add leakage-free geometry and annotation-derived features."""
    refined = df_features.copy()
    refined["log_distance"] = np.log1p(refined["distance"])
    refined["log_loop_span"] = np.log1p(refined["loop_span"])
    refined["distance_per_loop_span"] = refined["distance"] / refined["loop_span"].replace(0, np.nan)
    refined["anchor_width_sum"] = refined["anchor1_width"] + refined["anchor2_width"]
    refined["anchor_width_diff"] = (refined["anchor1_width"] - refined["anchor2_width"]).abs()
    refined["anchor_width_ratio"] = (
        refined[["anchor1_width", "anchor2_width"]].max(axis=1)
        / refined[["anchor1_width", "anchor2_width"]].min(axis=1).replace(0, np.nan)
    )
    refined["where_classification"] = refined["WHERE"].astype(str) + "_" + refined["classification"].astype(str)

    refined = refined.replace([np.inf, -np.inf], np.nan)
    refined = refined.dropna(subset=FEATURE_COLUMNS + [LABEL_COLUMN, GROUP_COLUMN]).copy()
    return refined


def build_preprocessor() -> ColumnTransformer:
    """Create preprocessing steps for refined numeric and categorical features."""
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


def compare_with_initial_metrics(refined_metrics: pd.DataFrame) -> pd.DataFrame:
    """Compare refined model metrics against the initial model run."""
    if not INITIAL_METRICS_PATH.exists():
        raise FileNotFoundError(f"Missing initial metrics: {INITIAL_METRICS_PATH}")

    initial_metrics = pd.read_csv(INITIAL_METRICS_PATH)
    initial_small = initial_metrics[
        ["model", "split", "balanced_accuracy", "average_precision"]
    ].rename(
        columns={
            "balanced_accuracy": "initial_balanced_accuracy",
            "average_precision": "initial_average_precision",
        }
    )
    refined_small = refined_metrics[
        ["model", "split", "balanced_accuracy", "average_precision"]
    ].rename(
        columns={
            "balanced_accuracy": "refined_balanced_accuracy",
            "average_precision": "refined_average_precision",
        }
    )

    comparison = initial_small.merge(refined_small, on=["model", "split"], how="inner")
    comparison["delta_balanced_accuracy"] = (
        comparison["refined_balanced_accuracy"] - comparison["initial_balanced_accuracy"]
    )
    comparison["delta_average_precision"] = (
        comparison["refined_average_precision"] - comparison["initial_average_precision"]
    )
    return comparison


def write_run_log(refined_metrics: pd.DataFrame, comparison: pd.DataFrame) -> None:
    """Write markdown log for refined feature run."""
    best_rows = refined_metrics.sort_values(["split", "balanced_accuracy"], ascending=[True, False])
    best_summary = best_rows.groupby("split").head(1)[["split", "model", "balanced_accuracy", "average_precision"]]
    delta_summary = comparison.sort_values("delta_balanced_accuracy", ascending=False)

    lines = [
        "## Refined Feature Run",
        "",
        f"Date: {date.today().isoformat()}",
        "",
        f"Input: `{INPUT_FEATURE_PATH}`",
        f"Refined feature table: `{REFINED_FEATURE_PATH}`",
        "",
        "Added features:",
        "",
        "- `log_distance`",
        "- `log_loop_span`",
        "- `distance_per_loop_span`",
        "- `anchor_width_sum`",
        "- `anchor_width_diff`",
        "- `anchor_width_ratio`",
        "- `where_classification`",
        "",
        "Best refined model by split:",
        "",
        best_summary.to_markdown(index=False),
        "",
        "Change from initial run:",
        "",
        delta_summary[
            [
                "model",
                "split",
                "delta_balanced_accuracy",
                "delta_average_precision",
            ]
        ].to_markdown(index=False),
        "",
        "Main caveat: These features are derived from existing loop geometry and annotations only. This remains exploratory and should not be framed as biological discovery.",
        "",
    ]

    REFINED_LOG_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Run refined feature training and compare against the initial run."""
    if not INPUT_FEATURE_PATH.exists():
        raise FileNotFoundError(f"Missing input feature table: {INPUT_FEATURE_PATH}")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    df_features = pd.read_csv(INPUT_FEATURE_PATH)
    refined_features = add_refined_features(df_features)
    refined_features[LABEL_COLUMN] = refined_features[LABEL_COLUMN].astype(int)
    refined_features.to_csv(REFINED_FEATURE_PATH, index=False)

    metrics_rows: list[dict[str, Any]] = []
    fitted_stratified_models: dict[str, Pipeline] = {}

    for split_name, (x_train, x_test, y_train, y_test) in make_splits(refined_features).items():
        for model_name, model in build_models().items():
            model.fit(x_train, y_train)
            metrics_rows.append(evaluate_model(model_name, split_name, model, x_test, y_test))

            if split_name == "stratified_80_20":
                fitted_stratified_models[model_name] = model

    refined_metrics = pd.DataFrame(metrics_rows)
    refined_metrics.to_csv(REFINED_METRICS_PATH, index=False)

    rf_model = fitted_stratified_models["RandomForest"]
    rf_features = extract_feature_names(rf_model)
    rf_importance = rf_model.named_steps["classifier"].feature_importances_
    df_rf_importance = pd.DataFrame({"feature": rf_features, "importance": rf_importance})
    df_rf_importance = df_rf_importance.sort_values("importance", ascending=False)
    df_rf_importance.to_csv(REFINED_RF_IMPORTANCE_PATH, index=False)

    logistic_model = fitted_stratified_models["LogisticRegression"]
    logistic_features = extract_feature_names(logistic_model)
    logistic_coef = logistic_model.named_steps["classifier"].coef_[0]
    df_logistic_coef = pd.DataFrame({"feature": logistic_features, "coefficient": logistic_coef})
    df_logistic_coef = df_logistic_coef.assign(abs_coefficient=lambda df: df["coefficient"].abs())
    df_logistic_coef = df_logistic_coef.sort_values("abs_coefficient", ascending=False)
    df_logistic_coef.to_csv(REFINED_LOGISTIC_COEF_PATH, index=False)

    comparison = compare_with_initial_metrics(refined_metrics)
    comparison.to_csv(REFINED_COMPARISON_PATH, index=False)
    write_run_log(refined_metrics, comparison)

    print(f"Wrote refined features: {REFINED_FEATURE_PATH}")
    print(f"Wrote refined metrics: {REFINED_METRICS_PATH}")
    print(refined_metrics[["model", "split", "balanced_accuracy", "average_precision"]].to_string(index=False))
    print(f"Wrote comparison: {REFINED_COMPARISON_PATH}")
    print(comparison[["model", "split", "delta_balanced_accuracy", "delta_average_precision"]].to_string(index=False))
    print(f"Wrote refined run log: {REFINED_LOG_PATH}")


if __name__ == "__main__":
    main()
