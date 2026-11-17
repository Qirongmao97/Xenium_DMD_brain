"""Plot the PCR and Western blot values supplied in Skipping_WB.ipynb.

Bars show means with standard errors, matching the notebook's plotting calls.
Dots are individual recorded values, not paired measurements across assays.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def plot_assay(data, ylabel, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    regions = list(dict.fromkeys(data["Region"]))
    groups = ["ASO", "U7ex51"]
    colors = {"ASO": "#8DA0CB", "U7ex51": "#E78AC3"}
    rng = np.random.default_rng(42)
    fig, ax = plt.subplots(figsize=(4, 3))
    xs = np.arange(len(regions))
    for i, group in enumerate(groups):
        values = [data.loc[data["Region"].eq(r) & data["Group"].eq(group), "Value"].to_numpy()
                  for r in regions]
        if any(len(v) < 2 for v in values):
            raise ValueError("Each assay/region/group needs at least two values for standard errors.")
        offset = xs + (i - 0.5) * 0.34
        ax.bar(offset, [v.mean() for v in values], 0.33,
               yerr=[v.std(ddof=1) / np.sqrt(len(v)) for v in values],
               color=colors[group], edgecolor="black", capsize=2, label=group)
        for x, y in zip(offset, values):
            ax.scatter(x + rng.uniform(-0.06, 0.06, len(y)), y, s=14, c="black", zorder=3)
    ax.set_xticks(xs, regions)
    ax.set(xlabel="Region", ylabel=ylabel)
    ax.legend(title="Group", frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, default=Path("inputs/pcr_western_blot.csv"))
    p.add_argument("--output-dir", type=Path, default=Path("results/assay_validation"))
    args = p.parse_args()
    data = pd.read_csv(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for assay, ylabel, name in (
        ("PCR", "Exon skipping rate (%)", "PCR_exon_skipping"),
        ("Western blot", "Dystrophin protein expression (%)", "Western_blot_expression"),
    ):
        sub = data.loc[data["Assay"].eq(assay)]
        if sub.empty:
            p.error(f"Missing measurements for {assay}.")
        plot_assay(sub, ylabel, args.output_dir / f"{name}.svg")


if __name__ == "__main__":
    main()
