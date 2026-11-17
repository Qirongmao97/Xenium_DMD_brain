"""Exploratory Scanpy ingest of Allen 10Xv2 regions into an HPF reference.

Preprocess the raw matrices and refit the reference using common variable
features before transferring its Leiden labels. This is not a Xenium
annotation workflow.
"""

import argparse
from pathlib import Path


def preprocess(path):
    import numpy as np
    import scanpy as sc

    data = sc.read_h5ad(path)
    if "gene_symbol" in data.var:
        data.var_names = data.var["gene_symbol"].astype(str).to_numpy()
    data.var_names_make_unique()
    sc.pp.filter_cells(data, min_genes=200)
    sc.pp.filter_genes(data, min_cells=3)
    if not data.n_obs or not data.n_vars:
        raise ValueError(f"No cells or genes remain after initial filtering in {path.name}.")
    data.var["mt"] = data.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(data, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    genes = data.obs["n_genes_by_counts"]
    lower, upper = np.quantile(genes, [0.02, 0.98])
    filtered = data[(genes > lower) & (genes < upper) & (data.obs["pct_counts_mt"] < 20)].copy()
    if filtered.n_obs < 3:
        raise ValueError(f"Too few cells after QC in {path.name}.")
    sc.pp.normalize_total(filtered, target_sum=1e4)
    sc.pp.log1p(filtered)
    sc.pp.highly_variable_genes(filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
    return filtered[:, filtered.var["highly_variable"]].copy()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True, help="Raw HPF H5AD.")
    p.add_argument("--queries", type=Path, nargs="+", required=True)
    p.add_argument("--output-dir", type=Path, default=Path("results/allen_reference/ingest"))
    args = p.parse_args()
    import anndata as ad
    import scanpy as sc

    paths = [args.reference, *args.queries]
    labels = [path.stem for path in paths]
    if len(set(labels)) != len(labels):
        p.error("Reference and query filenames must have distinct stems.")
    prepared = [preprocess(path) for path in paths]
    common = prepared[0].var_names
    for query in prepared[1:]:
        common = common.intersection(query.var_names, sort=False)
    if len(common) < 3:
        p.error("Fewer than three shared variable genes after preprocessing.")
    prepared = [data[:, common].copy() for data in prepared]
    for name, data in zip(labels, prepared):
        data.obs["region"] = name
        sc.pp.regress_out(data, ["total_counts", "pct_counts_mt"])
        sc.pp.scale(data, max_value=10)
    ref, *queries = prepared
    n_pcs = min(50, ref.n_obs - 1, ref.n_vars - 1)
    sc.tl.pca(ref, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(ref, n_neighbors=min(10, ref.n_obs - 1), n_pcs=min(20, n_pcs))
    sc.tl.leiden(ref, resolution=0.25, flavor="leidenalg")
    sc.tl.umap(ref)
    for query in queries:
        sc.tl.ingest(query, ref, obs="leiden")

    merged = ad.concat(prepared, join="inner", index_unique=":")
    composition = merged.obs.groupby(["region", "leiden"], observed=True).size().rename("cells").reset_index()
    composition["percent_of_region"] = 100 * composition["cells"] / composition.groupby(
        "region", observed=True)["cells"].transform("sum")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(args.output_dir / "allen_ingest.h5ad")
    composition.to_csv(args.output_dir / "allen_ingest_composition.csv", index=False)


if __name__ == "__main__":
    main()
