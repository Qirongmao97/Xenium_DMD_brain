"""Generate class and region marker dot plots and cell-composition tables."""

import argparse

from xenium_io import CORE_SAMPLES, add_xenium_arguments, log_normalize, merge_samples, open_sample

CLASS_ORDER = (
    "IT-ET Glut", "NP-CT-L6b Glut", "DG-IMN Glut", "CTX-CGE GABA", "CTX-MGE GABA",
    "CNU-LGE GABA", "CNU-HYa GABA", "CNU-HYa Glut", "HY GABA", "HY Glut",
    "MH-LH Glut", "TH Glut", "Astro-Epen", "OPC-Oligo", "Vascular", "Immune",
)
REGION_ORDER = (
    "Cortex L1", "Cortex L2/3", "Cortex L4", "Cortex L5", "Cortex L6a", "Cortex L6b",
    "Retrosplenial cortex", "Olfactory areas", "CA1", "CA2/3", "Dendate gyrus",
    "HPF field", "Thalamus 1", "Thalamus 2", "Reticular nucleus", "Lateral habenula",
    "Medial habenula", "Hypothalamus", "Basolateral amygdalar nucleus",
    "Central amygdalar nucleus", "Medial amygdalar nucleus", "Caudoputamen",
    "Ventricular system", "Fiber tracts", "Unassigned",
)


def run_markers(data, key, preferred, min_lfc, min_fraction, out):
    import pandas as pd
    import scanpy as sc

    subset = data[data.obs[key].notna()].copy()
    seen = list(pd.unique(subset.obs[key]))
    order = [name for name in preferred if name in seen]
    order += [name for name in seen if name not in order]
    subset.obs[key] = pd.Categorical(subset.obs[key], categories=order, ordered=True)
    if len(order) < 2 or subset.obs[key].value_counts().min() < 2:
        raise ValueError(f"{key} marker testing requires at least two cells per category and two categories.")
    sc.tl.rank_genes_groups(subset, key, method="wilcoxon", pts=True, use_raw=False)
    markers = sc.get.rank_genes_groups_df(subset, group=None)
    markers.to_csv(out / f"ranked_markers_{key.lower()}.csv", index=False)
    selected = markers.loc[
        markers["logfoldchanges"].gt(min_lfc)
        & markers["pct_nz_group"].gt(min_fraction)
        & markers["pct_nz_reference"].lt(0.5)
    ].sort_values("logfoldchanges", ascending=False)
    selected = selected.groupby("group", observed=True, sort=False).head(3).copy()
    selected["group"] = pd.Categorical(selected["group"], categories=order, ordered=True)
    selected = selected.sort_values("group", kind="stable")
    selected.to_csv(out / f"selected_markers_{key.lower()}.csv", index=False)
    genes = selected["names"].drop_duplicates().tolist()
    if genes:
        plot = sc.pl.dotplot(
            subset, var_names=genes, groupby=key, categories_order=order,
            standard_scale="var", dot_max=1.0, dot_min=0.0, smallest_dot=0.1,
            figsize=(max(4, 0.3 * len(genes)), 6 if key == "Class" else 8),
            show=False, return_fig=True,
        )
        plot.savefig(out / f"Dotplot_marker_{key.lower()}.svg")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    add_xenium_arguments(p, "markers", CORE_SAMPLES)
    p.add_argument("--composition-samples", nargs="+", metavar="ID",
                   help="Optional sample list for composition tables only.")
    args = p.parse_args()
    import matplotlib

    matplotlib.use("Agg")
    data = merge_samples(args.input_dir, args.samples, required=("Region", "Class"))
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)

    obs = data.obs
    if args.composition_samples:
        import pandas as pd

        if len(set(args.composition_samples)) != len(args.composition_samples):
            p.error("Composition sample IDs must be unique.")
        obs = pd.concat([
            open_sample(args.input_dir, sid, required=("Region", "Class")).obs
            for sid in args.composition_samples
        ])
    composition = (obs.groupby(["Sample", "Region", "Class"], observed=True)
                   .size().rename("n_cells").reset_index())
    totals = composition.groupby(["Sample", "Region"], observed=True)["n_cells"].transform("sum")
    composition["percent_of_sample_region"] = 100 * composition["n_cells"] / totals
    composition.to_csv(out / "sample_region_class_composition.csv", index=False)
    immune = composition.loc[composition["Class"].eq("Immune")].copy()
    if len(immune):
        totals = immune.groupby("Sample", observed=True)["n_cells"].transform("sum")
        immune["percent_of_sample_immune_cells"] = 100 * immune["n_cells"] / totals
    immune.to_csv(out / "immune_cells_by_sample_region.csv", index=False)

    log_normalize(data)
    run_markers(data, "Class", CLASS_ORDER, 1.0, 0.7, out)
    run_markers(data, "Region", REGION_ORDER, 0.5, 0.6, out)


if __name__ == "__main__":
    main()
