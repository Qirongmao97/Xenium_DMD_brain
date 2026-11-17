"""Export sample-level isoform expression and WT–mdx52 expression contrasts.

Combines Isoform_comparison and Expression_abundance. Gene expression is the
mean of log-normalized cell values, not a summed-count pseudobulk.
"""

import argparse

from xenium_io import (
    CORE_SAMPLES, EXON_FEATURES, ISOFORM_FEATURES, add_xenium_arguments,
    feature_frame, log_normalize, merge_samples, save_heatmap,
)


def summarize(expression, columns):
    keys = ["Sample", "Group", *columns]
    return expression.groupby(keys, observed=True, as_index=False).mean()


def sample_tests(sample_means, features):
    import numpy as np
    import pandas as pd
    from scipy.stats import ttest_ind

    rows = []
    for gene in features:
        wt = sample_means.loc[sample_means["Group"].eq("WT"), gene].dropna()
        mdx = sample_means.loc[sample_means["Group"].eq("mdx52"), gene].dropna()
        p = ttest_ind(wt, mdx, equal_var=False).pvalue if min(len(wt), len(mdx)) >= 2 else np.nan
        rows.append((gene, len(wt), len(mdx), mdx.mean() - wt.mean(), p))
    out = pd.DataFrame(rows, columns=["Feature", "n_WT", "n_mdx52", "mdx52_minus_WT", "p_Welch"])
    out["p_Bonferroni"] = (out["p_Welch"] * out["p_Welch"].notna().sum()).clip(upper=1)
    return out


def plot_sample_expression(means, features, title, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    if means.empty:
        return
    data = means.melt(id_vars=["Sample", "Group"], value_vars=features,
                      var_name="Feature", value_name="Expression")
    order = list(dict.fromkeys(data["Group"]))
    fig, ax = plt.subplots(figsize=(max(6, len(features) * 0.9), 4))
    sns.boxplot(data=data, x="Feature", y="Expression", hue="Group", hue_order=order,
                palette="Set2", showfliers=False, ax=ax)
    handles, labels = ax.get_legend_handles_labels()
    sns.stripplot(data=data, x="Feature", y="Expression", hue="Group", hue_order=order,
                  dodge=True, jitter=False, palette=dict.fromkeys(order, "black"),
                  size=3, alpha=0.65, ax=ax)
    ax.legend(handles, labels, title="Group", frameon=False)
    ax.set(title=title, ylabel="Mean log-normalized expression per sample", xlabel=None)
    ax.tick_params(axis="x", labelrotation=45)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main():
    import pandas as pd
    import re

    p = argparse.ArgumentParser(description=__doc__)
    add_xenium_arguments(p, "isoform_expression", CORE_SAMPLES)
    p.add_argument("--delta-feature", default="Dystrophin_Exon2_29")
    p.add_argument("--plot-region", action="append", default=[], help="Repeat for additional regions.")
    args = p.parse_args()
    data = merge_samples(args.input_dir, args.samples, required=("Region", "Class"))
    if "Group" not in data.obs or data.obs["Group"].isna().any():
        p.error("Each sample must have Group metadata or a known study sample ID.")
    features = [g for g in ISOFORM_FEATURES + EXON_FEATURES if g in data.var_names]
    if args.delta_feature not in features:
        p.error("--delta-feature must be an available isoform or exon-set feature.")
    # The original baseline comparison normalized WT and mdx52 together,
    # separately from the analyses that also included treated sections.
    baseline_mask = data.obs["Group"].isin(["WT", "mdx52"])
    baseline = None
    if baseline_mask.any():
        comparison_data = data[baseline_mask].copy()
        log_normalize(comparison_data)
        baseline = feature_frame(comparison_data, features).join(
            comparison_data.obs[["Sample", "Group", "Region", "Class"]]
        )
        del comparison_data
    log_normalize(data)
    expr = feature_frame(data, features).join(data.obs[["Sample", "Group", "Region", "Class"]])
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    sample_means = summarize(expr.drop(columns=["Region", "Class"]), [])
    region_means = summarize(expr.drop(columns="Class"), ["Region"])
    stratum_means = summarize(expr, ["Region", "Class"])
    sample_means.to_csv(out / "sample_expression.csv", index=False)
    region_means.to_csv(out / "sample_region_expression.csv", index=False)
    stratum_means.to_csv(out / "sample_region_class_expression.csv", index=False)
    plot_sample_expression(sample_means, features, "Whole brain", out / "sample_expression.svg")
    for region in args.plot_region:
        region_subset = region_means.loc[region_means["Region"].eq(region)]
        if region_subset.empty:
            p.error(f"No region labeled {region!r} in the selected samples.")
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", region)
        plot_sample_expression(region_subset, features, region, out / f"{safe}_sample_expression.svg")
    if "Dp427c" in expr:
        (expr.assign(detected=expr["Dp427c"].gt(0))
             .groupby(["Region", "Class"], observed=True)["detected"].sum()
             .unstack(fill_value=0).to_csv(out / "Dp427c_positive_cell_counts.csv"))

    if baseline is None:
        return
    sample_tests(summarize(baseline.drop(columns=["Region", "Class"]), []), features).to_csv(
        out / "WT_mdx52_sample_tests.csv", index=False
    )
    subset = summarize(baseline, ["Region", "Class"])
    group_means = subset.groupby(["Region", "Class", "Group"], observed=True)[args.delta_feature].mean()
    if group_means.empty:
        return
    wide = group_means.unstack("Group")
    if not {"WT", "mdx52"}.issubset(wide.columns):
        return
    delta = (wide["mdx52"] - wide["WT"]).unstack("Class")
    # The notebook masks strata at or below 0.5% of the total cells in their class.
    coverage = pd.crosstab(baseline["Region"], baseline["Class"])
    mask = coverage.gt(coverage.sum(axis=0) * 0.005, axis=1)
    delta = delta.where(mask.reindex_like(delta).fillna(False))
    delta.to_csv(out / f"{args.delta_feature}_mdx52_minus_WT.csv")
    if delta.notna().any().any():
        save_heatmap(delta, out / f"{args.delta_feature}_mdx52_minus_WT.svg",
                     palette="vlag", label="Mean log-normalized expression difference", center=0)


if __name__ == "__main__":
    main()
