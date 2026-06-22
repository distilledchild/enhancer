#!/usr/bin/env python3
"""Summarize exploratory ATAC support model results.

This script reads the model metrics and feature-importance outputs from
02_train_atac_support_models.py, then writes a compact markdown summary and
three figures for interpretation.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_DIR = Path("/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn")
DATA_PATH = PROJECT_DIR / "data" / "atac_enhancer_anchor_features.csv"
RESULTS_DIR = PROJECT_DIR / "results"

METRICS_PATH = RESULTS_DIR / "model_metrics.csv"
RF_IMPORTANCE_PATH = RESULTS_DIR / "feature_importance_random_forest.csv"
LOGISTIC_COEF_PATH = RESULTS_DIR / "logistic_regression_coefficients.csv"

SUMMARY_PATH = RESULTS_DIR / "model_result_summary.md"
METRIC_PLOT_PATH = RESULTS_DIR / "model_metric_comparison.png"
RF_PLOT_PATH = RESULTS_DIR / "top_random_forest_features.png"
LOGISTIC_PLOT_PATH = RESULTS_DIR / "top_logistic_coefficients.png"


def clean_feature_name(feature: str) -> str:
    """Return a reader-friendly feature name."""
    return feature.replace("numeric__", "").replace("categorical__", "")


def require_inputs() -> None:
    """Validate that all required result files exist before summarizing."""
    for path in [DATA_PATH, METRICS_PATH, RF_IMPORTANCE_PATH, LOGISTIC_COEF_PATH]:
        if not path.exists():
            raise FileNotFoundError(f"Missing required input: {path}")


def plot_metric_comparison(metrics: pd.DataFrame) -> None:
    """Plot balanced metrics for each model and split."""
    metric_columns = ["balanced_accuracy", "roc_auc", "average_precision"]
    model_order = ["DummyMostFrequent", "LogisticRegression", "RandomForest"]
    split_order = ["stratified_80_20", "chromosome_group_holdout"]
    split_labels = {
        "stratified_80_20": "Stratified",
        "chromosome_group_holdout": "Chr holdout",
    }

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharey=False)

    for axis, metric in zip(axes, metric_columns):
        x_positions = np.arange(len(model_order))
        bar_width = 0.36

        for index, split_name in enumerate(split_order):
            values = []
            for model_name in model_order:
                row = metrics[(metrics["model"] == model_name) & (metrics["split"] == split_name)]
                values.append(float(row.iloc[0][metric]))

            offset = (index - 0.5) * bar_width
            axis.bar(x_positions + offset, values, width=bar_width, label=split_labels[split_name])

        axis.set_title(metric.replace("_", " ").title())
        axis.set_xticks(x_positions)
        axis.set_xticklabels(["Dummy", "Logistic", "RF"], rotation=0)
        axis.set_ylim(0, 1)
        axis.grid(axis="y", alpha=0.25)

    axes[0].set_ylabel("Score")
    axes[-1].legend(frameon=False, loc="lower right")
    fig.suptitle("ATAC Support Model Metrics", y=1.03)
    fig.tight_layout()
    fig.savefig(METRIC_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_random_forest_importance(rf_importance: pd.DataFrame) -> None:
    """Plot top random forest feature importances."""
    plot_df = rf_importance.copy()
    plot_df["feature"] = plot_df["feature"].map(clean_feature_name)
    plot_df = plot_df.sort_values("importance", ascending=False).head(10)
    plot_df = plot_df.sort_values("importance", ascending=True)

    fig, axis = plt.subplots(figsize=(7, 4.5))
    axis.barh(plot_df["feature"], plot_df["importance"], color="#4C78A8")
    axis.set_xlabel("Importance")
    axis.set_title("Top Random Forest Features")
    axis.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(RF_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_logistic_coefficients(logistic_coef: pd.DataFrame) -> None:
    """Plot top logistic regression coefficients by absolute magnitude."""
    plot_df = logistic_coef.copy()
    plot_df["feature"] = plot_df["feature"].map(clean_feature_name)
    plot_df = plot_df.sort_values("abs_coefficient", ascending=False).head(10)
    plot_df = plot_df.sort_values("coefficient", ascending=True)
    colors = np.where(plot_df["coefficient"] >= 0, "#59A14F", "#E15759")

    fig, axis = plt.subplots(figsize=(7, 4.5))
    axis.barh(plot_df["feature"], plot_df["coefficient"], color=colors)
    axis.axvline(0, color="black", linewidth=0.8)
    axis.set_xlabel("Coefficient")
    axis.set_title("Top Logistic Regression Coefficients")
    axis.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(LOGISTIC_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_summary(feature_table: pd.DataFrame, metrics: pd.DataFrame, rf_importance: pd.DataFrame, logistic_coef: pd.DataFrame) -> None:
    """Write markdown interpretation of the first model run."""
    total_rows = len(feature_table)
    positive_rows = int(feature_table["label_atac_overlap"].sum())
    negative_rows = total_rows - positive_rows
    positive_rate = positive_rows / total_rows

    metric_table = metrics[
        [
            "model",
            "split",
            "n_test",
            "positive_rate_test",
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
        ]
    ].copy()

    rounded_metric_table = metric_table.copy()
    for column in [
        "positive_rate_test",
        "accuracy",
        "balanced_accuracy",
        "precision",
        "recall",
        "f1",
        "roc_auc",
        "average_precision",
    ]:
        rounded_metric_table[column] = rounded_metric_table[column].round(4)

    top_rf = rf_importance.copy()
    top_rf["feature"] = top_rf["feature"].map(clean_feature_name)
    top_rf = top_rf.head(8)

    top_logistic = logistic_coef.copy()
    top_logistic["feature"] = top_logistic["feature"].map(clean_feature_name)
    top_logistic = top_logistic.head(8)
    top_logistic = top_logistic[["feature", "coefficient", "abs_coefficient"]].copy()
    top_logistic["coefficient"] = top_logistic["coefficient"].round(4)
    top_logistic["abs_coefficient"] = top_logistic["abs_coefficient"].round(4)

    best_by_split = metrics.sort_values(["split", "balanced_accuracy"], ascending=[True, False])
    best_by_split = best_by_split.groupby("split").head(1)[["split", "model", "balanced_accuracy", "average_precision"]].copy()
    best_by_split["balanced_accuracy"] = best_by_split["balanced_accuracy"].round(4)
    best_by_split["average_precision"] = best_by_split["average_precision"].round(4)

    lines = [
        "# ATAC-supported enhancer anchors: model result summary",
        "",
        f"Date: {date.today().isoformat()}",
        "",
        "## Input",
        "",
        f"Feature table: `{DATA_PATH}`",
        "",
        f"Rows: {total_rows:,}",
        f"Positive labels: {positive_rows:,}",
        f"Negative labels: {negative_rows:,}",
        f"Positive label rate: {positive_rate:.4f}",
        "",
        "The label is imbalanced, so raw accuracy is not the main metric. Balanced accuracy, ROC-AUC, and average precision are more useful for judging whether the model improves over the dummy baseline.",
        "",
        "## Model Metrics",
        "",
        rounded_metric_table.to_markdown(index=False),
        "",
        "## Best Model By Split",
        "",
        best_by_split.to_markdown(index=False),
        "",
        "## Feature Interpretation",
        "",
        "Top random forest features:",
        "",
        top_rf.to_markdown(index=False),
        "",
        "Top logistic regression coefficients:",
        "",
        top_logistic.to_markdown(index=False),
        "",
        "## Conservative Interpretation",
        "",
        "LogisticRegression showed a modest improvement over DummyMostFrequent in balanced accuracy and average precision. RandomForest performed close to the baseline in this first pass.",
        "",
        "This should be described as an exploratory scikit-learn workflow rather than evidence for a strong predictive biological signal.",
        "",
        "Safe wording:",
        "",
        "> Built an exploratory scikit-learn workflow using regulatory genomics features to model ATAC-seq support for Hi-C-derived enhancer anchors, including dummy baseline comparison, logistic regression, random forest, balanced metrics, chromosome holdout evaluation, and feature-importance interpretation.",
        "",
    ]

    SUMMARY_PATH.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    """Create summary files and figures for the first model run."""
    require_inputs()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    feature_table = pd.read_csv(DATA_PATH)
    metrics = pd.read_csv(METRICS_PATH)
    rf_importance = pd.read_csv(RF_IMPORTANCE_PATH)
    logistic_coef = pd.read_csv(LOGISTIC_COEF_PATH)

    plot_metric_comparison(metrics)
    plot_random_forest_importance(rf_importance)
    plot_logistic_coefficients(logistic_coef)
    write_summary(feature_table, metrics, rf_importance, logistic_coef)

    print(f"Wrote summary: {SUMMARY_PATH}")
    print(f"Wrote metric plot: {METRIC_PLOT_PATH}")
    print(f"Wrote random forest plot: {RF_PLOT_PATH}")
    print(f"Wrote logistic plot: {LOGISTIC_PLOT_PATH}")


if __name__ == "__main__":
    main()
