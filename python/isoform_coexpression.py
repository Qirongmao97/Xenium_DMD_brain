"""Summarize long and short dystrophin probe co-detection by region and class.

The denominator is the number of cells with a signal in either probe group.
The final short-only panel preserves the notebook's two display thresholds.
"""

import argparse

from xenium_io import WT, add_xenium_arguments, feature_frame, merge_samples, save_heatmap

LONG = ("Dp427c", "Dp427m", "Dp427p")
SHORT = ("Dp71_Dp40_firstexon", "Dp40")


def summarize(expression, metadata):
    import pandas as pd

    a = expression[list(LONG)].gt(0).any(axis=1)
    b = expression[list(SHORT)].gt(0).any(axis=1)
    flags = pd.DataFrame({
        "Coexpression_Both": a & b, "Coexpression_Either": a | b,
        "Only_A": a & ~b, "Only_B": b & ~a, "expressionA": a, "expressionB": b,
    }).join(metadata[["Region", "Class"]])
    grouped = flags.groupby(["Region", "Class"], observed=True)
    stats = grouped.sum()
    stats["n_cells"] = grouped.size()
    denominator = stats["Coexpression_Either"].where(stats["Coexpression_Either"] > 0)
    for name, col in (("percent_A", "Only_A"), ("percent_B", "Only_B"),
                      ("percent_coexpressed", "Coexpression_Both")):
        stats[name] = 100 * stats[col] / denominator
    cls_total = stats.groupby(level="Class", observed=True)["n_cells"].transform("sum")
    stats["short_only_display"] = stats["percent_B"].where(
        (stats["n_cells"] > 0.005 * cls_total) & (stats["Only_B"] > 10)
    ).fillna(0)
    return stats


def main():
    p = argparse.ArgumentParser(description=__doc__)
    add_xenium_arguments(p, "coexpression", WT)
    args = p.parse_args()
    data = merge_samples(args.input_dir, args.samples, required=("Region", "Class"))
    stats = summarize(feature_frame(data, LONG + SHORT), data.obs)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    stats.to_csv(args.output_dir / "isoform_coexpression.csv")
    for col, name, palette, label, fmt in (
        ("Coexpression_Both", "Coexpressed_cells_boxplot", "Purples", "Cells positive for both", ".0f"),
        ("Coexpression_Either", "DMD_positive_boxplot", "Greens", "Cells positive for either", ".0f"),
        ("short_only_display", "short_isoforms", "Reds", "Short-only cells (% of either-positive)", ".1f"),
    ):
        matrix = stats[col].unstack("Class").fillna(0)
        save_heatmap(matrix, args.output_dir / f"{name}.svg", palette=palette, label=label,
                     annotate=True, fmt=fmt)


if __name__ == "__main__":
    main()
