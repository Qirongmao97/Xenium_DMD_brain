"""Recreate the exploratory Dp140 proxy from two Xenium exon-set signals.

This is a clipped difference in counts, not a direct isoform measurement.
The violin shows positive scaled values, as in the source notebook.
"""

import argparse
from pathlib import Path

from xenium_io import DEFAULT_INPUT, DEFAULT_OUTPUT, feature_frame, log_normalize, open_sample

FEATURES = ("Dystrophin_Exon56_62", "Dystrophin_Exon45_49")
PROXY = "Dp140_exon_difference_proxy"


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT / "dp140_proxy")
    p.add_argument("--sample", default="895")
    args = p.parse_args()

    import anndata as ad
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    from scipy.sparse import csr_matrix, hstack

    data = open_sample(args.input_dir, args.sample, required=("Region",))
    signal = feature_frame(data, FEATURES)
    proxy = (signal[FEATURES[0]] - signal[FEATURES[1]]).clip(lower=0)
    if PROXY in data.var_names:
        p.error(f"Input already contains the reserved feature {PROXY}.")
    matrix = hstack([csr_matrix(data.X), csr_matrix(proxy.to_numpy()[:, None])], format="csr")
    augmented = ad.AnnData(matrix, obs=data.obs.copy(),
                           var=pd.DataFrame(index=[*data.var_names, PROXY]))
    log_normalize(augmented)
    # Scaling a single column is equivalent to scaling the same column in the
    # full panel and avoids densifying the complete expression matrix.
    import scanpy as sc
    scored = augmented[:, [PROXY]].copy()
    sc.pp.scale(scored)
    output = data.obs[["Region"]].assign(
        proxy_counts=proxy, proxy_scaled=feature_frame(scored, [PROXY])[PROXY]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output_dir / f"Sample_{args.sample}_dp140_proxy.csv", index_label="cell_id")
    positive = output.loc[output["proxy_scaled"] > 0]
    if len(positive):
        fig, ax = plt.subplots(figsize=(10, 4))
        sns.violinplot(data=positive, x="Region", y="proxy_scaled", ax=ax, cut=0)
        ax.tick_params(axis="x", labelrotation=90)
        fig.savefig(args.output_dir / f"Sample_{args.sample}_dp140_proxy.svg", bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
