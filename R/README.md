# R spatial clustering

[Banksy_spatial_clustering.R](Banksy_spatial_clustering.R) clusters the four supplied sections on Xenium slide 1. It normalizes each Xenium assay separately, matches analyzer centroids by cell ID, runs BANKSY, and computes cluster markers and sample-by-cluster counts.

## Inputs

Place the inputs under `data/banksy/slide1/` or pass their paths on the command line.

| File | Required content |
| --- | --- |
| `slide1_g896.rds`, `slide1_g895.rds`, `slide1_g886.rds`, `slide1_g706.rds` | Prepared Seurat objects with a `Xenium` assay. Original cell IDs must match the analyzer metadata, either as column names or in a `cell_id` metadata field. |
| `cells.csv.gz` | Xenium Analyzer table for that slide containing unique `cell_id` values and numeric `x_centroid` and `y_centroid` coordinates. |

The script checks coordinate coverage before merging. It uses separate sample labels for BANKSY. The clustering settings from the original script are lambda 0.8, `k_geom = 10`, 30 PCs, and resolution 0.9.

## Run

Install `Seurat`, `data.table`, `R.utils` (for reading gzipped inputs), `Banksy`, and `SeuratWrappers` with BANKSY support. Use versions of SeuratObject and SeuratWrappers that are compatible with one another, particularly for the BANKSY assay data argument.

~~~r
install.packages(c("Seurat", "data.table", "R.utils", "BiocManager", "remotes"))
BiocManager::install("Banksy")
remotes::install_github("satijalab/seurat-wrappers")
~~~

From the repository root:

~~~bash
Rscript R/Banksy_spatial_clustering.R \
  data/banksy/slide1 data/banksy/slide1/cells.csv.gz results/banksy/slide1
~~~

The positional arguments are the RDS directory, the analyzer metadata file, and the output directory. They can all be omitted to use the paths above.

## Outputs

- `slide1_banksy_clustered.rds` contains the clustered Seurat object.
- `slide1_banksy_markers.csv` contains positive markers from the Xenium assay.
- `slide1_cluster_counts.csv` contains the number of cells per cluster and sample.

For optional Matrix Market conversion, see [the Python workflow](../python/README.md#xenium-and-assay-outputs). The exported files can be read with Seurat's `Read10X()`, but do not by themselves supply the Xenium coordinates or annotations required by this BANKSY script.
