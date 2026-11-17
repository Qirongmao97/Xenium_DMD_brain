"""Exploratory SpaGE prediction of Dmd in Allen MERFISH cells.

Provide normalized cells-by-genes input and a local upstream SpaGE checkout.
Predictions are not direct Xenium measurements.
"""

import argparse
from pathlib import Path
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", type=Path, required=True, help="Normalized Allen scRNA-seq H5AD.")
    p.add_argument("--spatial", type=Path, required=True, help="Normalized MERFISH cells-by-genes CSV.")
    p.add_argument("--spage-dir", type=Path, required=True, help="Upstream SpaGE/ directory containing main.py.")
    p.add_argument("--metadata", type=Path, help="Optional MERFISH metadata CSV indexed by cell_label.")
    p.add_argument("--section", help="Select a brain_section_label for the coordinate plots.")
    p.add_argument("--genes", nargs="+", default=["Dmd"])
    p.add_argument("--n-pv", type=int, default=460)
    p.add_argument("--output-dir", type=Path, default=Path("results/allen_reference/spage"))
    args = p.parse_args()
    if args.section and not args.metadata:
        p.error("--section requires --metadata.")
    if not (args.spage_dir / "main.py").is_file():
        p.error("--spage-dir must contain the upstream SpaGE main.py.")
    import anndata as ad
    import pandas as pd

    sys.path.insert(0, str(args.spage_dir.resolve()))
    from main import SpaGE

    spatial = pd.read_csv(args.spatial, index_col=0)
    spatial.index = spatial.index.astype(str)
    if not spatial.index.is_unique or not spatial.columns.is_unique:
        p.error("Spatial cell IDs and gene columns must be unique.")
    reference = ad.read_h5ad(args.reference, backed="r")
    try:
        shared = spatial.columns.intersection(reference.var_names)
        missing = set(args.genes) - set(reference.var_names)
        if missing:
            p.error(f"Requested genes missing from the reference: {', '.join(sorted(missing))}")
        if not 1 <= args.n_pv <= min(len(shared), len(spatial), reference.n_obs):
            p.error("--n-pv exceeds available shared genes or cells.")
        columns = list(dict.fromkeys([*shared, *args.genes]))
        ref = reference[:, columns].to_memory().to_df()
    finally:
        reference.file.close()

    prediction = SpaGE(spatial.loc[:, shared], ref, n_pv=args.n_pv, genes_to_predict=args.genes)
    if len(prediction) != len(spatial):
        raise ValueError("SpaGE did not return a prediction for each spatial cell.")
    prediction = pd.DataFrame(prediction)
    prediction.index = spatial.index
    prediction.index.name = "cell_label"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prediction.to_csv(args.output_dir / "spage_predictions.csv")
    if not args.metadata:
        return
    metadata = pd.read_csv(args.metadata, dtype={"cell_label": str}).set_index("cell_label")
    if args.section:
        metadata = metadata.loc[metadata["brain_section_label"].eq(args.section)]
    mapped = prediction.join(metadata[["x", "y"]], how="inner", validate="1:1").dropna(subset=["x", "y"])
    if mapped.empty:
        p.error("No predicted cells match the selected metadata coordinates.")
    mapped.to_csv(args.output_dir / "spage_predictions_with_coordinates.csv")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import re

    for gene in args.genes:
        fig, ax = plt.subplots(figsize=(8, 6))
        points = ax.scatter(-mapped["x"], -mapped["y"], c=mapped[gene], cmap="viridis", alpha=0.7, s=2)
        ax.set(aspect="equal", title=f"SpaGE prediction of {gene}")
        ax.set_axis_off()
        fig.colorbar(points, ax=ax, label="Predicted expression")
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", gene)
        fig.savefig(args.output_dir / f"{safe}_spatial_prediction.png", dpi=300, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
