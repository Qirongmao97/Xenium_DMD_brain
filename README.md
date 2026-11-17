# Evaluating exon skipping in the central nervous system in Duchenne muscular dystrophy using spatial transcriptomics

Analysis code for the mouse brain Xenium study in **iScience**. The repository includes regional and cell-class analyses of dystrophin probe signals, an exon-junction detection comparison, supporting assay plots, and spatial clustering with BANKSY.

## Start here

| Location | Contents |
| --- | --- |
| [Python workflows](python/README.md) | Xenium expression and co-detection, marker dot plots, junction-positive versus junction-negative comparisons, PCR/Western blot plots. |
| [R workflow](R/README.md) | BANKSY spatial clustering for one Xenium slide, using four prepared Seurat objects and Xenium Analyzer metadata. |
| [Allen Atlas workflows](python/README.md#allen-atlas-reference-and-exploratory-work) | Reference metadata and downsampling, plus separate exploratory Scanpy ingest and SpaGE scripts. These are not required for analyses that start from annotated Xenium H5ADs. |

Run commands from the repository root. Large data inputs and generated results are not included. The scripts use project-relative paths by default and accept overrides. The small PCR/Western blot input table is bundled so those two plots can run without Xenium data.

## Data

- **Xenium study data and Processed study H5ADs**. Start with [GEO GSE318507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE318507), the accession listed in the published paper.
- **Adult mouse scRNA-seq reference**. The paper links the [Allen Brain Cell Atlas WMB-10Xv2 matrices](https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#expression_matrices/WMB-10Xv2/20230630/). The Atlas scripts also use the metadata from release `20230830`.

The Python study analyses start from per-section H5ADs with expression and annotations already prepared. Confirm that their probe names and metadata match the [Python input schema](python/README.md#xenium-inputs). The BANKSY script starts from prepared RDSs. The code does not implement raw Xenium segmentation or a complete raw-to-annotated-Xenium pipeline.

## Scope of the cleanup

The original notebooks contained repeat analyses and alternative plots. The Python scripts retain the main workflows and export source tables, but do not attempt to reproduce every interactive cell or the exact styling of the published figures. Where the old drafts used different statistical units, the script documentation states which unit is used. The Atlas ingest and SpaGE examples remain exploratory.

## Citation and license

Mao Q, Ahmadi A, de Vries S, Heezen LGM, Vacca O, Doisy M, Aartsma-Rus A, van Putten M, Goyenvalle A, Mahfouz A, Spitali P. *Evaluating exon skipping in the central nervous system in Duchenne muscular dystrophy using spatial transcriptomics*. **iScience** 29(8), 116906 (2026). [doi:10.1016/j.isci.2026.116906](https://doi.org/10.1016/j.isci.2026.116906).

The repository software is distributed under the included [GNU GPL v3 license](LICENSE). Refer to the source repositories for terms governing external datasets and tools.
