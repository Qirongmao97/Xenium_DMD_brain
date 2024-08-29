library(Seurat)
library(Banksy)
library(R.utils)
library(data.table)
library(dplyr)
library(tidyverse)
library(SeuratWrappers)
library(scater)
library(ggplot2)



###load the seurat objects for each of the samples on slide1 
g896 <- readRDS('/exports/archive/hg-exon-skip/alireza/Xenium/slide1/slide1_g896.rds')
g895 <- readRDS('/exports/archive/hg-exon-skip/alireza/Xenium/slide1/slide1_g895.rds')
g886 <- readRDS('/exports/archive/hg-exon-skip/alireza/Xenium/slide1/slide1_g886.rds')
g706 <- readRDS('/exports/archive/hg-exon-skip/alireza/Xenium/slide1/slide1_g706.rds')


####Normalise each dataset separately before merging
g896 <- NormalizeData(g896)
g895 <- NormalizeData(g895)
g886 <- NormalizeData(g886)
g706 <- NormalizeData(g706)


####Merge all the samples together
slide1 <- merge(g896 , c(g895,g886,g706))






####import the metadata provided by xenium analyzer
cells <- fread("/exports/archive/hg-exon-skip/Qirong/Xenium/output-XETG00101__0019132__Region_1__20240123__092349/cells.csv.gz")
meta <- slide1@meta.data
cells <- cells %>% column_to_rownames(., var = "cell_id")



#Merge the metadata datasets together
meta_data <- merge(meta , cells , by = "row.names")
meta_data <- meta_data %>% column_to_rownames(., var = "Row.names")
columns_to_keep <- c("sample","x_centroid","y_centroid")
meta.data <- meta_data[,columns_to_keep]
names(meta.data)[names(meta.data) == "sample"] <- "sample_id"







#####add the metadata including the spatial coordinates to merged object
slide1 <- AddMetaData(slide1, metadata = meta.data)






#####run BANKSY algorithm on seurat object and optimise parameters, 
####use higher lambdas to identify regions and lower lambda to cluster cell types
slide1 = RunBanksy(slide1, lambda = 0.8, assay = 'Xenium', slot = 'data',
                   dimx = 'x_centroid', dimy = 'y_centroid', features = 'all',
                   group = 'sample_id', split.scale = TRUE, k_geom = 10)

####Run the common seurat pipeline on the merged object 
slide1 <- RunPCA(slide1, assay = 'BANKSY', features = rownames(slide1), npcs = 30)
slide1 <- RunUMAP(slide1, dims = 1:30)
slide1 <- FindNeighbors(slide1, dims = 1:30)
slide1 <- FindClusters(slide1, resolution = 0.9)



####check the marker genes fr the identified clusters
DefaultAssay(slide1) <- "Xenium"
markers <- FindAllMarkers(slide1 , only.pos = T , logfc.threshold = 0.1)
markers




####chech the number of cells assigned to each cluster for each sample
clusters_slide1 <- slide1@meta.data %>%
  group_by(sample, seurat_clusters) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0)


####save the results for the clustering and marker genes
saveRDS(slide1 , "/exports/archive/hg-exon-skip/alireza/Xenium/slide1/slide1-886-lambda1.rds")
write.csv(arranged.markers , "/exports/archive/hg-exon-skip/alireza/Xenium/slide1/markers_arranged-886.csv")









