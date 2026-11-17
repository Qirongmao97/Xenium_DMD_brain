"""Shared input, study sample IDs, and plotting helpers for Xenium analyses."""

from pathlib import Path

WT = ("895", "894", "898")
MDX52 = ("896", "897", "18")
TCDNA = ("886", "890", "891")
U7 = ("4", "5", "3")
CORE_SAMPLES = WT + MDX52 + TCDNA + U7
GROUP_BY_SAMPLE = {
    **dict.fromkeys(WT, "WT"),
    **dict.fromkeys(MDX52, "mdx52"),
    **dict.fromkeys(TCDNA, "tcDNA"),
    **dict.fromkeys(U7, "U7ex51"),
}
ISOFORM_FEATURES = ("Dp427c", "Dp427m", "Dp427p", "Dp427p2", "Dp71_Dp40_firstexon", "Dp40")
EXON_FEATURES = tuple(
    "Dystrophin_Exon" + suffix for suffix in ("2_29", "30_44", "45_49", "56_62", "79")
)
JUNCTION_FEATURES = ("DystrophinExon50_53", "DystrophinExon50_51", "DystrophinExon51_53")
DEFAULT_INPUT = Path("data/xenium")
DEFAULT_OUTPUT = Path("results/xenium")


def add_xenium_arguments(parser, out_name, samples):
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT / out_name)
    parser.add_argument("--samples", nargs="+", default=samples, metavar="ID")


def open_sample(input_dir, sample_id, required=()):
    import anndata as ad

    path = input_dir / f"Sample_{sample_id}.h5ad"
    if not path.is_file():
        raise FileNotFoundError(f"Missing sample H5AD: {path}")
    sample = ad.read_h5ad(path)
    missing = set(required) - set(sample.obs.columns)
    if missing:
        raise ValueError(f"{path.name} is missing obs columns: {', '.join(sorted(missing))}")
    if not sample.n_obs or not sample.var_names.is_unique or not sample.obs_names.is_unique:
        raise ValueError(f"{path.name} must have cells and unique cell IDs and feature names.")
    sample.obs["Sample"] = sample_id
    if sample_id in GROUP_BY_SAMPLE:
        sample.obs["Group"] = GROUP_BY_SAMPLE[sample_id]
    return sample


def merge_samples(input_dir, sample_ids, required=()):
    import anndata as ad

    if len(set(sample_ids)) != len(sample_ids):
        raise ValueError("Supply each sample ID only once.")
    samples = {sid: open_sample(input_dir, sid, required) for sid in sample_ids}
    data = ad.concat(samples, join="inner", index_unique=":")
    if not data.n_vars:
        raise ValueError("No common Xenium features across the selected samples.")
    return data


def require_features(data, features):
    absent = [gene for gene in features if gene not in data.var_names]
    if absent:
        raise ValueError(f"Missing Xenium features: {', '.join(absent)}")


def feature_frame(data, features):
    import numpy as np
    import pandas as pd
    from scipy.sparse import issparse

    require_features(data, features)
    values = data[:, list(features)].X
    values = values.toarray() if issparse(values) else np.asarray(values)
    return pd.DataFrame(values, index=data.obs_names, columns=features)


def log_normalize(data, target_sum=None):
    import scanpy as sc

    sc.pp.normalize_total(data, target_sum=target_sum)
    sc.pp.log1p(data)


def save_heatmap(matrix, path, *, palette, label, center=None, annotate=False, fmt=".1f"):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(max(5, 0.35 * len(matrix.columns)), max(3, 0.25 * len(matrix))))
    sns.heatmap(matrix, ax=ax, cmap=palette, center=center, annot=annotate,
                fmt=fmt, cbar_kws={"label": label})
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
