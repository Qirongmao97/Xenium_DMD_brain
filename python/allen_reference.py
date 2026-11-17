"""Prepare or summarize the Allen whole-mouse-brain 10Xv2 reference.

Subcommands assemble release metadata, sample 10% of cells within each class,
merge prepared reference shards, or summarize measured Dmd expression.
"""

import argparse
import json
from pathlib import Path

DEFAULT_ROOT = Path("data/allen")
DEFAULT_OUT = Path("results/allen_reference")
NT_GENES = ("Slc17a7", "Slc17a6", "Slc17a8", "Slc32a1", "Slc6a5",
            "Slc18a3", "Slc6a3", "Slc6a4", "Slc6a2")


def build_metadata(data_dir, manifest_path):
    import pandas as pd

    with manifest_path.open() as handle:
        release = json.load(handle)["file_listing"]
    meta = release["WMB-10X"]["metadata"]
    taxonomy = release["WMB-taxonomy"]["metadata"]

    def read(section, name, index, **kwargs):
        path = data_dir / section[name]["files"]["csv"]["relative_path"]
        return pd.read_csv(path, index_col=index, **kwargs)

    cells = read(meta, "cell_metadata", "cell_label", dtype={"cell_label": str})
    details = read(taxonomy, "cluster_to_cluster_annotation_membership_pivoted",
                   "cluster_alias", keep_default_na=False)
    colors = read(taxonomy, "cluster_to_cluster_annotation_membership_color", "cluster_alias")
    roi = read(meta, "region_of_interest_metadata", "acronym").rename(
        columns={"order": "region_of_interest_order", "color_hex_triplet": "region_of_interest_color"}
    )
    return (cells.join(details, on="cluster_alias", validate="m:1")
            .join(colors, on="cluster_alias", validate="m:1")
            .join(roi[["region_of_interest_order", "region_of_interest_color"]],
                  on="region_of_interest_acronym", validate="m:1"))


def class_sample(metadata, fraction, seed):
    import pandas as pd

    metadata = metadata.loc[metadata["library_method"].eq("10Xv2")]
    groups = [
        group.sample(frac=fraction, random_state=seed)
        for _, group in metadata.groupby("class", observed=True)
    ]
    if not groups:
        raise ValueError("No 10Xv2 cells with class annotations in the metadata.")
    return pd.concat(groups)


def merge_reference(paths, selected=None):
    import anndata as ad

    matrices = []
    for path in paths:
        backed = ad.read_h5ad(path, backed="r")
        try:
            if selected is None:
                matrix = backed.to_memory()
            else:
                keep = backed.obs_names.isin(selected.index)
                if not keep.any():
                    continue
                matrix = backed[keep, :].to_memory()
        finally:
            backed.file.close()
        if "gene_symbol" in matrix.var:
            matrix.var_names = matrix.var["gene_symbol"].astype(str).to_numpy()
        matrix.var_names_make_unique()
        if selected is not None:
            matrix.obs["cluster_ann"] = True
        matrix.obs["source_matrix"] = path.stem
        matrices.append(matrix)
    if not matrices:
        raise ValueError("No matching cells were found in the input H5AD files.")
    out = ad.concat(matrices, join="inner", index_unique=None)
    if not out.obs_names.is_unique or not out.n_vars:
        raise ValueError("Reference matrices must have disjoint cell IDs and shared genes.")
    return out


def summarize_expression(matrix_path, metadata_path, output_dir, matrix_label=None):
    import anndata as ad
    import pandas as pd

    from xenium_io import save_heatmap

    metadata = pd.read_csv(metadata_path, index_col=0)
    metadata.index = metadata.index.astype(str)
    if not metadata.index.is_unique:
        raise ValueError("Atlas metadata cell IDs must be unique.")
    if matrix_label:
        metadata = metadata.loc[metadata["feature_matrix_label"].eq(matrix_label)]
    data = ad.read_h5ad(matrix_path, backed="r")
    try:
        symbols = data.var["gene_symbol"] if "gene_symbol" in data.var else data.var_names
        gene_ids = pd.Series(data.var_names, index=pd.Index(symbols))
        wanted = [gene for gene in ("Dmd", *NT_GENES) if gene in gene_ids]
        if "Dmd" not in wanted or gene_ids.loc[wanted].index.has_duplicates:
            raise ValueError("The input must contain a unique Dmd and unique selected gene symbols.")
        keep = data.obs_names.isin(metadata.index)
        if not keep.any():
            raise ValueError("No expression cells match the selected metadata.")
        expr = data[keep, gene_ids.loc[wanted].tolist()].to_memory().to_df()
    finally:
        data.file.close()
    expr.columns = wanted
    joined = metadata.join(expr, how="inner", validate="1:1")
    output_dir.mkdir(parents=True, exist_ok=True)
    for field in ("neurotransmitter", "class", "subclass", "supertype"):
        if field not in joined:
            continue
        means = joined.groupby(field, observed=True)[wanted].mean().sort_values("Dmd", ascending=False)
        means.to_csv(output_dir / f"expression_by_{field}.csv")
        if len(means):
            n_show = 8 if field in ("class", "subclass") else 10
            save_heatmap(means[["Dmd"]].head(n_show), output_dir / f"Dmd_by_{field}.svg",
                         palette="magma_r", label="Mean input expression")
            markers = [g for g in wanted if g != "Dmd"]
            if field == "neurotransmitter" and markers:
                save_heatmap(means[markers], output_dir / "neurotransmitter_markers.svg",
                             palette="magma_r", label="Mean input expression")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    subs = p.add_subparsers(dest="command", required=True)
    m = subs.add_parser("metadata")
    m.add_argument("--data-dir", type=Path, default=DEFAULT_ROOT)
    m.add_argument("--manifest", type=Path)
    m.add_argument("--output", type=Path, default=DEFAULT_OUT / "cell_metadata_with_cluster_annotation.csv")

    s = subs.add_parser("subsample")
    s.add_argument("--input-dir", type=Path, required=True)
    s.add_argument("--metadata", type=Path, default=DEFAULT_OUT / "cell_metadata_with_cluster_annotation.csv")
    s.add_argument("--fraction", type=float, default=0.1)
    s.add_argument("--seed", type=int, default=42)
    s.add_argument("--output", type=Path, default=DEFAULT_OUT / "file_sampled.h5ad")

    c = subs.add_parser("merge")
    c.add_argument("--inputs", type=Path, nargs="+", required=True)
    c.add_argument("--output", type=Path, default=DEFAULT_OUT / "merged_sc.h5ad")

    v = subs.add_parser("summary")
    v.add_argument("--matrix", type=Path, required=True)
    v.add_argument("--metadata", type=Path, default=DEFAULT_OUT / "cell_metadata_with_cluster_annotation.csv")
    v.add_argument("--matrix-label", help="Optional Allen feature_matrix_label, e.g. WMB-10Xv2-TH.")
    v.add_argument("--output-dir", type=Path, default=DEFAULT_OUT / "expression_summary")
    args = p.parse_args()
    if args.command == "metadata":
        manifest = args.manifest or args.data_dir / "releases" / "20230830" / "manifest.json"
        result = build_metadata(args.data_dir, manifest)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(args.output)
    elif args.command == "summary":
        summarize_expression(args.matrix, args.metadata, args.output_dir, args.matrix_label)
    else:
        if args.command == "subsample":
            import pandas as pd

            if not 0 < args.fraction <= 1:
                p.error("--fraction must be in (0, 1].")
            metadata = pd.read_csv(args.metadata, index_col=0)
            metadata.index = metadata.index.astype(str)
            selected = class_sample(metadata, args.fraction, args.seed)
            paths = sorted(args.input_dir.glob("*raw.h5ad"))
            if not paths:
                p.error("No *raw.h5ad files in --input-dir.")
        else:
            paths, selected = args.inputs, None
        result = merge_reference(paths, selected)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        result.write_h5ad(args.output)


if __name__ == "__main__":
    main()
