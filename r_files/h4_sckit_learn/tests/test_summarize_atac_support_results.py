#!/usr/bin/env python3
"""Validate summary outputs from the ATAC support model workflow."""

from pathlib import Path


PROJECT_DIR = Path("/Users/pete/Desktop/playground/enhancer/r_files/h4_sckit_learn")
RESULTS_DIR = PROJECT_DIR / "results"

SUMMARY_PATH = RESULTS_DIR / "model_result_summary.md"
METRIC_PLOT_PATH = RESULTS_DIR / "model_metric_comparison.png"
RF_PLOT_PATH = RESULTS_DIR / "top_random_forest_features.png"
LOGISTIC_PLOT_PATH = RESULTS_DIR / "top_logistic_coefficients.png"


def require_nonempty_file(path: Path) -> None:
    if not path.exists():
        raise AssertionError(f"Missing expected output: {path}")
    if path.stat().st_size == 0:
        raise AssertionError(f"Output file is empty: {path}")


def main() -> None:
    for path in [SUMMARY_PATH, METRIC_PLOT_PATH, RF_PLOT_PATH, LOGISTIC_PLOT_PATH]:
        require_nonempty_file(path)

    summary_text = SUMMARY_PATH.read_text(encoding="utf-8")
    required_phrases = [
        "ATAC-supported enhancer anchors",
        "DummyMostFrequent",
        "LogisticRegression",
        "RandomForest",
        "imbalanced",
        "exploratory scikit-learn workflow",
    ]
    missing_phrases = [phrase for phrase in required_phrases if phrase not in summary_text]
    if missing_phrases:
        raise AssertionError(f"Summary is missing phrases: {missing_phrases}")

    for path in [METRIC_PLOT_PATH, RF_PLOT_PATH, LOGISTIC_PLOT_PATH]:
        if path.stat().st_size < 1000:
            raise AssertionError(f"Plot file may be invalid or too small: {path}")

    print(f"PASS summary: {SUMMARY_PATH}")
    print(f"PASS metric plot: {METRIC_PLOT_PATH}")
    print(f"PASS random forest plot: {RF_PLOT_PATH}")
    print(f"PASS logistic plot: {LOGISTIC_PLOT_PATH}")


if __name__ == "__main__":
    main()
