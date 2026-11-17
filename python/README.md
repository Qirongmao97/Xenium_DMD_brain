# Python analysis workflows

These command-line scripts replace the mixed Jupyter notebooks. Run them from the repository root. Each script accepts `--help`. The [repository README](../README.md) links the paper and data records.

## Environment

Create an isolated Python environment and install the listed dependencies:

~~~bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r python/requirements.txt
~~~

On Windows, activate the environment with `.venv\Scripts\activate`. SpaGE has an additional external source requirement. See the Atlas section below. For a reproducible rerun, record the exact versions you used rather than assuming the latest packages match the original notebook environment.

## Xenium inputs

Unless a script says otherwise, place each prepared H5AD in `data/xenium/` with the name `Sample_<ID>.h5ad`. `_input-dir_` examples below use `--input-dir` to choose another folder. The scripts that compute signal detection or normalize counts expect **raw, unlogged Xenium counts in `adata.X`**, unique cell and feature names, and the named probe features used by the analysis. `Region` and `Class` are required `obs` columns for most Xenium workflows. The loader derives `Sample` from the filename and uses the cohort mapping below for the 12 known IDs.

| Study group | Sample IDs |
| --- | --- |
| Wild-type | 895, 894, 898 |
| Untreated mdx52 | 896, 897, 18 |
| tcDNA / ASO | 886, 890, 891 |
| U7ex51 | 4, 5, 3 |

Options such as `--samples` can select a subset of sections. For sample-level expression analysis, extra sections must have `Group` metadata. The marker script can take a separate `--composition-samples` list to include additional sections without changing the marker analysis cohort.

## Xenium and assay outputs

| Script | Source notebook(s) | What it produces |
| --- | --- | --- |
| [isoform_expression.py](isoform_expression.py) | `Isoform_comparison.ipynb`, `Expression_abundance.ipynb` | Sample, sample-by-region, and sample-by-region-by-class mean expression tables. Also produces plots, a WT–mdx52 difference heatmap, and sample-level Welch tests with Bonferroni adjustment. |
| [isoform_coexpression.py](isoform_coexpression.py) | `Coexpression_Isoforms.ipynb` | WT long/short probe co-detection counts and percentages by region and class, plus heatmaps. The denominator for the percentages is cells with either signal. |
| [xenium_markers.py](xenium_markers.py) | `Xenium_marker_dotplot.ipynb` | Scanpy Wilcoxon marker rankings and class/region dot plots, plus cell composition and immune-cell distributions. |
| [skipping_differential_expression.py](skipping_differential_expression.py) | `Skipping_DE.ipynb` | ASO or U7ex51 matched-control comparisons, nominal-significance frequency, mean log fold changes, and matching diagnostics. |
| [pcr_western_blot.py](pcr_western_blot.py) | `Skipping_WB.ipynb` | Two SVGs with mean ± standard error and individual values from the bundled [assay table](../inputs/pcr_western_blot.csv). |
| [dp140_proxy.py](dp140_proxy.py) | `Dp140.ipynb` | Exploratory clipped exon-set difference and scaled proxy by region for a selected section. This is not a direct Dp140 measurement. |
| [export_seurat_matrix.py](export_seurat_matrix.py) | `Conversion_h5ad_to_seurat.ipynb` | Gzipped feature × cell Matrix Market file, barcodes, and features. It does not export AnnData annotations or create a Seurat object. |

~~~bash
python python/isoform_expression.py --input-dir data/xenium --plot-region CA1
python python/isoform_coexpression.py
python python/xenium_markers.py
python python/skipping_differential_expression.py --treatment ASO
python python/skipping_differential_expression.py --treatment U7ex51
python python/pcr_western_blot.py
~~~

Default outputs go under `results/xenium/` by analysis, except the assay figures, which go to `results/assay_validation/`. The matrix export defaults to `results/seurat_matrix/`. Use `--output-dir` on each script to redirect files.

### Interpretation and reproducibility

- The expression summaries average log-normalized cell values within each sample or stratum. They are not summed-count pseudobulk profiles. The WT–mdx52 branch normalizes the baseline samples together and the whole-brain tests use **one value per sample**, not repeated regions from a sample as independent replicates.
- Junction-positive status uses a nonzero `DystrophinExon50_53` signal. The matched comparison uses the original class/region-first logic, with class, region, and global fallbacks. The defaults request 1,000 draws and permit each negative cell in at most five completed draws. Always check `matching_diagnostics.csv` for completed versus attempted draws and fallback use. The frequency of nominal `p < 0.05` cell-level tests is descriptive, not a permutation p-value or an animal-level test.
- The Dp140 script preserves the clipped difference of the two supplied exon-set signals but stores it under a distinct proxy name. The per-cell CSV includes both the clipped counts and the scaled plotting values.
- The two validation plots use standard errors, as specified in the notebook plotting calls. The exported table does not assign a paired identity to points across assays.

## Allen Atlas reference and exploratory work

| Script | Source notebook(s) | Purpose |
| --- | --- | --- |
| [allen_reference.py](allen_reference.py) | `whole-brain.ipynb`, `Untitled.ipynb`, `merging_allen_brain.ipynb`, reference-merge portion of `merge.ipynb` | Join Atlas cell metadata with taxonomy and ROI tables, subsample 10Xv2 cells within class, merge prepared inputs, and summarize measured Dmd and neurotransmitter marker expression. |
| [allen_ingest.py](allen_ingest.py) | `simple_scanpy_integration.ipynb` | Exploratory QC, alignment on shared variable genes, and transfer of HPF reference Leiden labels into other Allen regions. Not used to assign Xenium cell classes here. |
| [allen_spage_imputation.py](allen_spage_imputation.py) | SpaGE portion of `merge.ipynb` | Exploratory expression prediction in MERFISH cells from a normalized Allen scRNA-seq reference. The imputed output is distinct from measured Dmd expression and from the Xenium workflow. |

For the metadata step, download the Allen `20230830` manifest and the CSVs it references into a local directory, preserving the release-relative file paths. For subsampling, provide a directory of raw WMB-10Xv2 H5ADs. For expression summaries, provide an expression H5AD and the metadata table. The `summary` command uses the matrix's expression values as supplied and does not renormalize them.

~~~bash
python python/allen_reference.py metadata --data-dir data/allen
python python/allen_reference.py subsample --input-dir data/allen/raw
python python/export_seurat_matrix.py --input results/allen_reference/file_sampled.h5ad
~~~

For the optional SpaGE command, provide normalized expression, a spatial CSV with cells as rows and genes as columns, and the directory containing `main.py` from the [upstream SpaGE repository](https://github.com/tabdelaal/SpaGE). Install the upstream dependencies in the environment. Select `--n-pv` to fit the number of shared genes and cells. The notebook used 460. If you supply MERFISH metadata, it must contain `cell_label`, `x`, and `y`. To filter coordinate plots to one section, also provide `brain_section_label` in the metadata and set `--section`.

~~~bash
python python/allen_spage_imputation.py \
  --reference data/allen/normalized_reference.h5ad \
  --spatial data/allen/merfish_expression.csv \
  --spage-dir external/SpaGE/SpaGE
~~~

## Smoke tests

~~~bash
python -m unittest discover -s tests
~~~

The smoke tests use synthetic inputs to check co-detection denominators, sample-level averaging, matched-draw rules, reference subsampling, and Matrix Market orientation. They do not replace a rerun on the study data.
