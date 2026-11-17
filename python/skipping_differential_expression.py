"""Compare junction-positive and junction-negative cells by matched resampling.

Reuses the notebook's Class/Region match, class/region/global fallbacks,
re-use cap, and cell-level Wilcoxon tests. The output is a sensitivity
summary, not a replicate-level or permutation p-value.
"""

import argparse

import numpy as np
import pandas as pd

from xenium_io import (
    EXON_FEATURES, JUNCTION_FEATURES, TCDNA, U7,
    add_xenium_arguments, feature_frame, log_normalize, merge_samples,
)

JUNCTION = JUNCTION_FEATURES[0]


def negative_pools(obs, negatives, keys):
    strata = obs.iloc[negatives].groupby(keys, observed=True, sort=False).indices
    return {key: negatives[positions] for key, positions in strata.items()}


def matched_draw(strata, exact, by_class, by_region, negatives, usage, max_usage, rng):
    """Take one set of controls, without duplicate cells within the draw."""
    used = np.zeros(len(usage), dtype=bool)
    selected = []
    matched_counts = dict(exact=0, class_only=0, region_only=0, global_pool=0)
    empty = np.empty(0, dtype=int)
    for (cls, reg), requested in strata.items():
        needed = int(requested)
        for tier, pool in (
            ("exact", exact.get((cls, reg), empty)), ("class_only", by_class.get(cls, empty)),
            ("region_only", by_region.get(reg, empty)), ("global_pool", negatives),
        ):
            candidates = pool[(~used[pool]) & (usage[pool] < max_usage)]
            count = min(len(candidates), needed)
            if count:
                chosen = rng.choice(candidates, size=count, replace=False)
                used[chosen] = True
                selected.extend(chosen)
                matched_counts[tier] += count
                needed -= count
            if not needed:
                break
        if needed:
            return None, matched_counts
    return np.asarray(selected, dtype=int), matched_counts


def save_effect_plot(summary, treatment, iterations, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if treatment == "ASO":
        display = summary.reindex([*JUNCTION_FEATURES, "Dp427c", "Dp427m", "Dp427p",
                                   "Dp71_Dp40_firstexon", "Dp40",
                                   *[g for g in EXON_FEATURES if g != "Dystrophin_Exon45_49"]])
    else:
        eligible = summary[summary["Significance_Count"] > 0.9 * iterations]
        if eligible.empty:
            return
        display = eligible.sort_values("Avg_LogFC", key=lambda x: x.abs(), ascending=False).head(30)
        display = display.sort_values("Avg_LogFC", ascending=False)
    display = display.dropna(subset=["Avg_LogFC"])
    if display.empty:
        return
    fig, (a, b) = plt.subplots(1, 2, sharey=True, figsize=(8, max(3, len(display) * 0.25)))
    a.scatter(display["Significance_Count"], display.index, c="darkred", s=15)
    b.barh(display.index, display["Avg_LogFC"],
           color=np.where(display["Avg_LogFC"] < 0, "#f8766d", "#00bfc4"))
    a.invert_yaxis()
    a.set_xlabel("Draws with nominal p < 0.05")
    b.set_xlabel("Mean Scanpy log fold change")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    add_xenium_arguments(p, "skipping_de", None)
    p.add_argument("--treatment", choices=("ASO", "U7ex51"), required=True)
    p.add_argument("--iterations", type=int, default=1000)
    p.add_argument("--max-control-uses", type=int, default=5)
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()
    if min(args.iterations, args.max_control_uses) < 1:
        p.error("Iteration and control re-use limits must be positive.")
    import scanpy as sc

    ids = args.samples or (TCDNA if args.treatment == "ASO" else U7)
    data = merge_samples(args.input_dir, ids, required=("Class", "Region"))
    if data.obs[["Class", "Region"]].isna().any().any():
        p.error("Matching requires nonmissing Class and Region for every cell.")
    positive = feature_frame(data, [JUNCTION])[JUNCTION].gt(0).to_numpy()
    pos, neg = np.flatnonzero(positive), np.flatnonzero(~positive)
    if min(len(pos), len(neg)) < 2:
        p.error("At least two junction-positive and two junction-negative cells are required.")
    log_normalize(data, target_sum=1e4)
    metadata = data.obs
    strata = metadata.iloc[pos].groupby(["Class", "Region"], observed=True).size()
    exact = negative_pools(metadata, neg, ["Class", "Region"])
    by_class = negative_pools(metadata, neg, "Class")
    by_region = negative_pools(metadata, neg, "Region")
    usage = np.zeros(data.n_obs, dtype=int)
    rng = np.random.default_rng(args.seed)
    genes = pd.Index(data.var_names, name="Gene")
    effect_sum = np.zeros(len(genes))
    finite_counts = np.zeros(len(genes), dtype=int)
    sig_counts = np.zeros(len(genes), dtype=int)
    diagnostics = []
    convergence = []
    previous_significant, previous_top = set(), set()
    completed = 0
    for iteration in range(1, args.iterations + 1):
        selected, counts = matched_draw(strata, exact, by_class, by_region, neg,
                                        usage, args.max_control_uses, rng)
        diagnostics.append({"iteration": iteration, "completed": selected is not None, **counts})
        if selected is None:
            continue
        usage[selected] += 1
        draw = data[np.hstack((pos, selected))].copy()
        draw.obs["pop_vs_rest"] = pd.Categorical(np.repeat(["pop1", "rest"], [len(pos), len(selected)]))
        sc.tl.rank_genes_groups(draw, "pop_vs_rest", groups=["pop1"], reference="rest",
                               method="wilcoxon", n_genes=draw.n_vars, use_raw=False)
        rank = sc.get.rank_genes_groups_df(draw, group="pop1").set_index("names").reindex(genes)
        effects = rank["logfoldchanges"].to_numpy()
        finite = np.isfinite(effects)
        effect_sum[finite] += effects[finite]
        finite_counts += finite
        sig_counts += rank["pvals"].lt(0.05).to_numpy()
        completed += 1
        if iteration % 10 == 0:
            significant = set(genes[sig_counts >= 10])
            average = pd.Series(np.divide(effect_sum, finite_counts,
                                out=np.full(len(genes), np.nan), where=finite_counts > 0), index=genes)
            top = set(average.dropna().nlargest(30).index)

            def jaccard(a, b):
                return len(a & b) / len(a | b) if a or b else 1.0

            convergence.append({
                "iteration": iteration, "jaccard_significant": jaccard(previous_significant, significant),
                "jaccard_logfc": jaccard(previous_top, top), "n_sig_genes": len(significant),
            })
            previous_significant, previous_top = significant, top

    out = args.output_dir / args.treatment
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(diagnostics).to_csv(out / "matching_diagnostics.csv", index=False)
    pd.DataFrame(convergence, columns=["iteration", "jaccard_significant", "jaccard_logfc",
                                      "n_sig_genes"]).to_csv(out / "jaccard_convergence.csv", index=False)
    counts = pd.crosstab(metadata["Class"], pd.Series(
        np.where(positive, "Positive", "Negative"), index=metadata.index)
    )
    counts.to_csv(out / "class_counts_by_junction_status.csv")
    counts.div(counts.sum(axis=0), axis=1).mul(100).to_csv(out / "class_percent_by_junction_status.csv")
    if not completed:
        raise RuntimeError(f"No control draws completed. See {out / 'matching_diagnostics.csv'}.")
    mean = np.divide(effect_sum, finite_counts, out=np.full(len(genes), np.nan), where=finite_counts > 0)
    summary = pd.DataFrame({"Avg_LogFC": mean, "Significance_Count": sig_counts,
                            "Completed_Draws": completed}, index=genes)
    summary = summary.sort_values("Significance_Count", ascending=False, kind="stable")
    summary.to_csv(out / "skipping_de_summary.csv")
    save_effect_plot(summary, args.treatment, args.iterations, out / f"Skipping_DE_{args.treatment}.svg")
    print(f"{args.treatment}: {completed}/{args.iterations} matched draws completed.")


if __name__ == "__main__":
    main()
