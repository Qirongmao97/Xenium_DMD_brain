# Cluster sections from one Xenium slide with BANKSY.
#
# Inputs: four Seurat RDS files with a Xenium assay and the slide's cells.csv.gz
# from Xenium Analyzer. Cell IDs in the Seurat objects or their cell_id
# metadata must match the Analyzer cell_id field.
#
# Usage from the repository root:
#   Rscript R/Banksy_spatial_clustering.R [rds_dir] [cells_file] [output_dir]
# By default, inputs are under data/banksy/slide1 and outputs under
# results/banksy/slide1. Saves the clustered Seurat object, cluster markers,
# and sample-by-cluster cell counts.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 3L || any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript R/Banksy_spatial_clustering.R [rds_dir] [cells_file] [output_dir]\n")
  quit(status = if (length(args) > 3L) 1L else 0L)
}
rds_dir <- if (length(args) >= 1L) args[[1L]] else file.path("data", "banksy", "slide1")
cells_file <- if (length(args) >= 2L) args[[2L]] else file.path(rds_dir, "cells.csv.gz")
output_dir <- if (length(args) >= 3L) args[[3L]] else file.path("results", "banksy", "slide1")

sample_ids <- c("g896", "g895", "g886", "g706")
paths <- setNames(file.path(rds_dir, paste0("slide1_", sample_ids, ".rds")), sample_ids)
if (any(!file.exists(paths)) || !file.exists(cells_file)) {
  missing <- c(paths[!file.exists(paths)], cells_file[!file.exists(cells_file)])
  stop("Missing inputs: ", paste(missing, collapse = ", "), call. = FALSE)
}

cells <- data.table::fread(cells_file, colClasses = c(cell_id = "character"))
needed <- c("cell_id", "x_centroid", "y_centroid")
if (!all(needed %in% names(cells)) || anyNA(cells$cell_id) || anyDuplicated(cells$cell_id)) {
  stop("Analyzer metadata needs unique cell_id values and both centroid columns.", call. = FALSE)
}
if (!is.numeric(cells$x_centroid) || !is.numeric(cells$y_centroid)) {
  stop("Analyzer centroid columns must be numeric.", call. = FALSE)
}

samples <- lapply(sample_ids, function(sid) {
  obj <- readRDS(paths[[sid]])
  if (!inherits(obj, "Seurat") || !"Xenium" %in% Assays(obj)) {
    stop(sid, " must be a Seurat object with a Xenium assay.", call. = FALSE)
  }
  ids <- if ("cell_id" %in% names(obj[[]])) as.character(obj$cell_id) else colnames(obj)
  index <- match(ids, cells$cell_id)
  if (anyNA(index) || anyDuplicated(ids) ||
      any(!is.finite(cells$x_centroid[index])) || any(!is.finite(cells$y_centroid[index]))) {
    stop("Cannot match unique, finite Analyzer coordinates for every cell in ", sid, ".", call. = FALSE)
  }
  md <- data.frame(
    sample_id = sid,
    analyzer_cell_id = ids,
    x_centroid = cells$x_centroid[index],
    y_centroid = cells$y_centroid[index],
    row.names = colnames(obj)
  )
  obj <- AddMetaData(obj, metadata = md)
  NormalizeData(obj, assay = "Xenium")
})
slide_cell_ids <- unlist(lapply(samples, function(obj) obj$analyzer_cell_id), use.names = FALSE)
if (anyDuplicated(slide_cell_ids)) {
  stop("Analyzer cell IDs overlap between the four input sections.", call. = FALSE)
}

slide <- merge(samples[[1L]], y = samples[-1L], add.cell.ids = sample_ids, merge.data = TRUE)
if (inherits(slide[["Xenium"]], "Assay5")) {
  slide <- JoinLayers(slide, assay = "Xenium")
}

set.seed(42)
banksy_args <- list(
  object = slide, lambda = 0.8, assay = "Xenium", dimx = "x_centroid",
  dimy = "y_centroid", features = "all", group = "sample_id", split.scale = TRUE,
  k_geom = 10
)
# SeuratWrappers versions differ in the name used for the assay data argument.
data_argument <- if ("layer" %in% names(formals(SeuratWrappers::RunBanksy))) "layer" else "slot"
banksy_args[[data_argument]] <- "data"
slide <- do.call(SeuratWrappers::RunBanksy, banksy_args)
slide <- RunPCA(slide, assay = "BANKSY", features = rownames(slide[["BANKSY"]]), npcs = 30)
slide <- RunUMAP(slide, dims = 1:30)
slide <- FindNeighbors(slide, dims = 1:30)
slide <- FindClusters(slide, resolution = 0.9)

DefaultAssay(slide) <- "Xenium"
markers <- FindAllMarkers(slide, assay = "Xenium", only.pos = TRUE, logfc.threshold = 0.1)
cluster_counts <- as.data.frame.matrix(table(slide$seurat_clusters, slide$sample_id))
cluster_counts <- data.frame(cluster = rownames(cluster_counts), cluster_counts, check.names = FALSE)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(output_dir)) stop("Could not create the output directory.", call. = FALSE)
saveRDS(slide, file.path(output_dir, "slide1_banksy_clustered.rds"))
write.csv(markers, file.path(output_dir, "slide1_banksy_markers.csv"), row.names = FALSE)
write.csv(cluster_counts, file.path(output_dir, "slide1_cluster_counts.csv"), row.names = FALSE)
