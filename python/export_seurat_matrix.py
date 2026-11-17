"""Export an H5AD matrix in feature-by-cell Matrix Market format for Seurat."""

import argparse
import gzip
from pathlib import Path


def export_matrix(data, output_dir, layer=None):
    from scipy import io
    from scipy.sparse import csr_matrix

    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError("Cell barcodes and feature names must be unique.")
    matrix = data.X if layer is None else data.layers[layer]
    if matrix is None:
        raise ValueError("The selected expression matrix is missing.")
    output_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_dir / "matrix.mtx.gz", "wb") as handle:
        io.mmwrite(handle, csr_matrix(matrix.T))
    with gzip.open(output_dir / "barcodes.tsv.gz", "wt", encoding="utf-8") as handle:
        handle.writelines(f"{barcode}\n" for barcode in data.obs_names)
    with gzip.open(output_dir / "features.tsv.gz", "wt", encoding="utf-8") as handle:
        handle.writelines(f"{gene}\t{gene}\tGene Expression\n" for gene in data.var_names)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, default=Path("results/seurat_matrix"))
    p.add_argument("--layer", help="AnnData layer to export instead of X, e.g. counts.")
    args = p.parse_args()
    import anndata as ad

    export_matrix(ad.read_h5ad(args.input), args.output_dir, args.layer)


if __name__ == "__main__":
    main()
