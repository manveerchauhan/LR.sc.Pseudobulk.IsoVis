# script to perform pseudo-bulk DEGs for LR sc data and export counts for IsoVis
setwd("/data/gpfs/projects/punim0646/manveer/IsoVis_pseudobulk")

library(Seurat)
library(tidyverse)
library(DoubletFinder)
library(grid)
library(gridExtra)
library(reticulate)
library(openxlsx)
library(HGNChelper)


# load scType's gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load scType's cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# Import our feature_ID_converter python script
source_python("/data/gpfs/projects/punim0646/manveer/Scripts/scRNAseq-Analysis-pipeline/feature_ID_converter.py")
gtfFilePath <- "/data/gpfs/projects/punim0646/manveer/gencode.v43.chr_patch_hapl_scaff.annotation.gtf"


set.seed(48)
### Functions-----
## Function to filter out low quality cells from seurat objects based on feature counts per cell
RunSeuratQC <- function(seu.obj){
  ## calculate max and min feature threshold points------------------------------------
  max.features <- round(mean(seu.obj$nFeature_RNA) + (1.5 * sd(seu.obj$nFeature_RNA)))
  min.features <- round(mean(seu.obj$nFeature_RNA) - (1.5 * sd(seu.obj$nFeature_RNA)))
  
  ## Filter seurat object with calculated min/max feature and count thresholds-------
  seu.filtered <- subset(seu.obj, subset = nFeature_RNA > min.features & 
                           nFeature_RNA < max.features &
                           nCount_RNA > 800)
  
  VlnPlot(seu.obj, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
  VlnPlot(seu.filtered, features = c('nFeature_RNA', 'nCount_RNA', ncol = 2))
  
  return(seu.filtered)
}

## Function to run the standard seurat workflow processing steps
seuratWorkflowSteps <- function(seu.obj, npc = 20, 
                                cluster_resolution = 0.4, geneLvl = TRUE,
                                seu.obj.Qced.geneLvl = NULL) {
  cluster_resolution.figs = list()
  
  if(geneLvl){
    seu.filtered <- RunSeuratQC(seu.obj)
    print("Low quality cells removed from gene level object")
  } else{
    seu.filtered <- subset(seu.obj, cells = row.names(seu.obj.Qced.geneLvl@meta.data))
    if(identical(rownames(seu.filtered@meta.data), rownames(seu.obj.Qced.geneLvl@meta.data))){
      print("Barcodes in isoform level object match the QCed gene level object.")
    }
  }
  
  ## Run standard seurat workflow steps----
  seu.filtered <- NormalizeData(seu.filtered)
  seu.filtered <- FindVariableFeatures(seu.filtered, 
                                       selection.method = 'vst',
                                       nfeatures = 2000)
  seu.filtered <- ScaleData(seu.filtered, rownames(seu.filtered))
  seu.filtered <- RunPCA(seu.filtered,
                         features = VariableFeatures(object = seu.filtered))
  print(ElbowPlot(seu.filtered))
  
  seu.filtered <- FindNeighbors(seu.filtered, dims = 1:npc)
  seu.filtered <- FindClusters(seu.filtered, resolution = cluster_resolution)
  seu.filtered <- RunUMAP(seu.filtered, dims = 1:npc)
  
  return(seu.filtered)
}

## Function to use scType to automatically annotate cell types within a seurat object
find_cellTypes <- function(seurat.obj){
  # We need to convert gene ids to their names for this to work, this needs to be done manually if working with FLAMES
  gene_ids <- rownames(seurat.obj[["RNA"]]@scale.data)
  gene_names <- ListofGeneIDstoNames(gene_ids, gtfFilePath)
  
  # Set the gene names as the new row names in the RNA assay's scale.data matrix
  rownames(seurat.obj[["RNA"]]@scale.data) <- gene_names
  
  # Load the desired marker database file
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  tissue = "Brain"
  
  # prepare gene sets
  gs_list = gene_sets_prepare(db_, tissue)
  gs_list
  
  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = seurat.obj[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
  
  # merge by cluster-----------------------------------------------------
  cL_results = do.call("rbind", lapply(unique(seurat.obj@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(seurat.obj@meta.data[seurat.obj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat.obj@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[,1:3])
  
  # View UMAP with cell type labels
  seurat.obj@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]; 
    seurat.obj@meta.data$customclassif[seurat.obj@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  
  labeledUMAP <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif') + 
    labs(title = "Cell-type Annotations", color = "Labels \n(from scType)") +
    theme(plot.title = element_text(size = 20))
  
  baseUMAP <- DimPlot(seurat.obj, reduction = "umap", label = TRUE, repel = TRUE) + 
    labs(title = "Unsupervised Clustering", color = "cluster \n(from PCA)") +
    theme(plot.title = element_text(size = 20))
  
  labeledUMAP | baseUMAP
  
  ### Create pie charts displaying sctypes confidence scores for each unsupervised cluster's cell type annotation----
  # Code used here is adapted from : https://github.com/IanevskiAleksandr/sc-type/blob/master/README.md
  
  # Prepare the edges
  cL_results = cL_results[order(cL_results$cluster),]
  edges = cL_results
  edges$type = paste0(edges$type, "_", edges$cluster)
  edges$cluster = paste0("cluster ", edges$cluster)
  edges = edges[, c("cluster", "type")]
  colnames(edges) = c("from", "to")
  rownames(edges) = NULL
  
  # Prepare nodes
  nodes_lvl1 = sctype_scores[, c("cluster", "ncells")]
  nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster)
  nodes_lvl1$Colour = "#f1f1ef"
  nodes_lvl1$ord = 1
  nodes_lvl1$realname = nodes_lvl1$cluster
  nodes_lvl1 = as.data.frame(nodes_lvl1)
  nodes_lvl2 = data.frame(cluster = character(),
                          ncells = numeric(),
                          Colour = character(),
                          ord = numeric(),
                          realname = character(),
                          stringsAsFactors = FALSE)
  
  ccolss = c("#5f75ae", "#92bbb8", "#64a841", "#e5486e", "#de8e06", "#eccf5a", "#b5aa0f", "#e4b680", "#7ba39d", "#b15928", "#ffff99", "#6a3d9a", "#cab2d6", "#ff7f00", "#fdbf6f", "#e31a1c", "#fb9a99", "#33a02c", "#b2df8a", "#1f78b4", "#a6cee3")
  
  for (i in 1:length(unique(cL_results$cluster))) {
    dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]
    nodes_tmp = data.frame(cluster = paste0(dt_tmp$type, "_", dt_tmp$cluster),
                           ncells = dt_tmp$scores,
                           Colour = ccolss[i],
                           ord = 2,
                           realname = dt_tmp$type,
                           stringsAsFactors = FALSE)
    nodes_lvl2 = rbind(nodes_lvl2, nodes_tmp)
  }
  
  nodes = rbind(nodes_lvl1, nodes_lvl2)
  nodes$ncells[nodes$ncells < 1] = 1
  files_db = openxlsx::read.xlsx(db_)[, c("cellName", "shortName")]
  files_db = unique(files_db)
  nodes = merge(nodes, files_db, all.x = TRUE, all.y = FALSE, by.x = "realname", by.y = "cellName", sort = FALSE)
  nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]
  nodes = nodes[, c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
  
  # Identify duplicate rows based on the "cluster" column
  duplicated_rows <- duplicated(nodes$cluster)
  # Filter the dataframe to keep only the first occurrence of each unique value
  nodes <- nodes[!duplicated_rows, ]
  
  # Create a new dataframe by filtering out rows with ord equal to 2, then remove unwanted cols
  filtered_nodes <- subset(nodes, ord != 1)
  filtered_nodes <- subset(filtered_nodes, select = -c(4, 5))
  
  # Extract the number after "_" and replace the string with the extracted number in the cluster column
  filtered_nodes$cluster <- str_replace(filtered_nodes$cluster, ".*_([0-9]+)", "Cluster \\1")
  
  # Group the dataframe by 'cluster' column
  grouped_df <- aggregate(ncells ~ realname + cluster, data = filtered_nodes, FUN = sum)
  
  # Calculate proportions within each cluster group
  grouped_df <- filtered_nodes %>%
    group_by(cluster) %>%
    mutate(prop = ncells / sum(ncells)) %>%
    ungroup() %>% 
    arrange(cluster)
  
  # Creat a list to store each pie chart generated
  pieChartList <- list()
  # Create a separate pie chart for each unique cluster group
  for (cluster in unique(grouped_df$cluster)) {
    # Subset the data for the current cluster
    cluster_data <- grouped_df[grouped_df$cluster == cluster, ]
    
    # Create the pie chart
    pie_chart <- ggplot(cluster_data, aes(x = "", y = prop, fill = realname)) +
      geom_bar(stat = "identity", width = 1) +
      coord_polar("y", start = 0) +
      labs(title = paste(cluster)) +
      theme_void() +
      theme(legend.position = "right") +
      labs(fill = "Labels Considered")  # Rename the legend
    
    # Add the percentages to the pie chart with white color, rounded to 1 decimal place
    pie_chart <- pie_chart + 
      geom_text(aes(label = paste0(round(prop * 100, 1), "%"), x = 1.5), 
                position = position_stack(vjust = 0.5), size = 4)
    
    # Append the pie chart to the list
    pieChartList <- c(pieChartList, list(pie_chart))
  }
  pieChartList
  print(length(pieChartList))
  
  ### Organize what the function returns------------------
  outputs <- list(baseUMAP, labeledUMAP, seurat.obj, pieChartList)
  
  return(outputs)
}

## Function to export sctype annotation confidence pie charts into a pdf
generate_PDF <- function(piechartlist){
  
  # Open a PDF file for exporting pie charts and UMAPs
  fileName <- paste0("piechartAnnotationConfidence.pdf")
  pdf(fileName, width = 15, height = 10)  # Adjust width and height as needed
  
  # Add pie charts for scTypeDB labels
  piechart.layout <- grid.arrange(grobs = piechartlist, nrow = 3, ncol = 2,
                                  top = textGrob("Pie Charts denoting cell label confidence scores for each cluster\n"),
                                  gp = gpar(fontsize = 20))
  grid.draw(piechart.layout)
  
  # Close the PDF device
  dev.off()
}

## Function that returns a list of all differentially expressed features across cell types in a seurat object
find.cellType.DE.markers <- function(seu.obj, cellType.list, 
                                     num.top.isoforms = 10){
  # Find all the upregulated differentially expressed isoforms per cell type
  Idents(seu.obj) <- seu.obj$customclassif
  CellType.lr.all.markers <- FindAllMarkers(seu.obj,
                                            logfc.threshold = 0.25,
                                            min.pct = 0.1,
                                            only.pos = TRUE)
  
  top.markers.per.celltype <- data.frame()
  for(i in 1:length(cellType.list)) {
    # print the current celltype being viewed
    print(cellType.list[i])
    
    # Filter and extract top n upregulated isoforms for current cell type
    filtered.celltype.markers <- CellType.lr.all.markers %>%
      dplyr::filter(cluster == cellType.list[i] & p_val_adj < 0.05) %>%
      arrange(desc(abs(avg_log2FC)))
    
    
    top.markers.per.celltype <- bind_rows(top.markers.per.celltype,
                                          head(filtered.celltype.markers, num.top.isoforms))
  }
  View(top.markers.per.celltype)
  
  return(top.markers.per.celltype)
}


# Import raw count matrices from BLAZE paper -----
geneMatrixFilePath = "/data/gpfs/projects/punim0646/manveer/lr-cortdiff-day25_gene_count.csv"
isoMatrixFilePath = "/data/gpfs/projects/punim0646/manveer/FLAMES_Day-25+55_promethION_LSK110/lr-cortdiff-day25_transcript_count.csv"

## Read in the FLAMES transcript count csv files------
count_matrix.isoformLvl <- read.csv(isoMatrixFilePath, row.names = 1) # isoform level count matrix
count_matrix.geneLvl <- read.csv(geneMatrixFilePath, row.names = 1) # gene level count matrix (from BLAZE paper's python converter script)
count_matrix.geneLvl <- count_matrix.geneLvl[, -1] # remove first col of dataframe with transcript IDs

## Create our Seurat objects (one gene level and one isoform level)----
## (using features that appear in at least 3 cells, and cells that contain >= 1 feature)
seu.obj.geneLvl <- CreateSeuratObject(counts = count_matrix.geneLvl, project = "LR_D25_BLAZE",
                                 min.cells = 3, min.features = 1)
seu.obj.isoformLvl <- CreateSeuratObject(counts = count_matrix.isoformLvl, project = "LR_D25_BLAZE",
                                              min.cells = 3, min.features = 1)





# Run standard Seurat QC and preprocessing steps---------
seu.processed.geneLvl <- seuratWorkflowSteps(seu.obj.geneLvl, cluster_resolution = 0.6)
seu.processed.isoLvl <- seuratWorkflowSteps(seu.obj.isoformLvl, 
                                            cluster_resolution = 0.6,
                                            geneLvl = FALSE,
                                            seu.obj.Qced.geneLvl = seu.processed.geneLvl)
base.isoLvl.UMAP <- DimPlot(seu.processed.isoLvl, reduction = 'umap', label = TRUE) + 
  labs(title = "Unsupervised Clustering", color = "cluster \n(from PCA)") +
  theme(plot.title = element_text(size = 20))

# Automatically annotate cell types using sctype------------
seu.annotated.geneLvl.info <- find_cellTypes(seu.processed.geneLvl)
seu.annotated.geneLvl <- seu.annotated.geneLvl.info[[3]]
generate_PDF(seu.annotated.geneLvl.info[[4]])

## Based on scType confidence scores - relabel probable false positive annotations-----
View(seu.annotated.geneLvl@meta.data)
# Replace false positive annotations in the customclassif column
seu.annotated.geneLvl@meta.data <- seu.annotated.geneLvl@meta.data %>%
  mutate(customclassif = case_when(
    customclassif == "Oligodendrocyte precursor cells" ~ "Neural Progenitor cells",
    customclassif == "Myelinating Schwann cells" ~ "Neuroepithelial cells",
    TRUE ~ customclassif
  ))
# transfer cell annotations to isoform level seurat object
seu.processed.isoLvl$customclassif <- seu.annotated.geneLvl$customclassif

relabeled.seu.geneLvl.UMAP <- DimPlot(seu.annotated.geneLvl, reduction = "umap", 
                                                               label = TRUE, repel = TRUE, group.by = 'customclassif') + 
  labs(title = "Corrected Cell-type Annotations", color = "Labels \n(from scType)") +
  theme(plot.title = element_text(size = 20))
seu.annotated.geneLvl.info[[1]] | relabeled.seu.geneLvl.UMAP

# View unsupervised cluster UMAPS and annotated UMAPs-----------
annotated.isoLvl.UMAP <- DimPlot(seu.processed.isoLvl, reduction = "umap", 
                                 label = TRUE, repel = TRUE, group.by = 'customclassif') + 
  labs(title = "Cell-type Annotations", color = "Labels \n(from scType)") +
  theme(plot.title = element_text(size = 20))

base.isoLvl.UMAP | annotated.isoLvl.UMAP # isoform level object UMAPs
seu.annotated.geneLvl.info[[1]] | seu.annotated.geneLvl.info[[2]] # gene level object UMAPs
seu.annotated.geneLvl.info[[1]] | relabeled.seu.geneLvl.UMAP


## Aggregate counts across annotated cell types at the isoform and gene level--------
aggregated.counts.isoLvl <- AggregateExpression(seu.processed.isoLvl, 
                                     group.by = 'customclassif')
aggregated.counts.isoLvl <- as.data.frame(aggregated.counts.isoLvl$RNA)
View(aggregated.counts.isoLvl)

aggregated.counts.geneLvl <- AggregateExpression(seu.annotated.geneLvl,
                                                 group.by = 'customclassif')
aggregated.counts.geneLvl <- as.data.frame(aggregated.counts.geneLvl$RNA)


## Find top 10 differentially expressed genes (DEGs) per cell type at the single cell level (by p value)----
top.10.DEgenes.per.celltype <- find.cellType.DE.markers(seu.annotated.geneLvl,
                                                        cellType.list = colnames(aggregated.counts.isoLvl))
View(top.10.DEgenes.per.celltype)


## Filter our aggregated gene dataframe such that only the upregulated genes per cell type we identified remain
filtered.pseudobulk.geneLvl <- aggregated.counts.geneLvl[rownames(aggregated.counts.geneLvl) %in% rownames(top.10.DEgenes.per.celltype), ]
View(filtered.pseudobulk.geneLvl)


## Create a transcript id and gene id dictionary----------------------------------
# the genematrix file here is the output of the blaze paper's isoform count matrix converter python script
geneMatrixFilePath <- "/data/gpfs/projects/punim0646/manveer/lr-cortdiff-day25_gene_count.csv"
trans_gene_dict <- read.csv(geneMatrixFilePath, row.names = 1)[, 1, drop = FALSE]
trans_gene_dict$transcript_id <- gsub("(\\d)([BE])", "\\1,\\2", trans_gene_dict$transcript_id)
View(trans_gene_dict)

## Create a dictionary showing gene ids for each transcript id, then create a gene-id column for aggregated isoform pseudobulk matrix----
# Assuming trans_gene_dict has rownames as gene IDs and a column transcript_id with comma-separated values
trans_gene_expanded <- trans_gene_dict %>%
  rownames_to_column(var = "gene_id") %>%
  separate_rows(transcript_id, sep = ",") %>%
  distinct()

# Create a named vector for mapping
transcript_to_gene_map <- setNames(trans_gene_expanded$gene_id, trans_gene_expanded$transcript_id)

aggregated.counts.isoLvl$gene_id <- transcript_to_gene_map[rownames(aggregated.counts.isoLvl)]
View(aggregated.counts.isoLvl)

### Format and export pseudobulk isoforms for DEGs as csv file for IsoVis---------
pseudobulk.isoCounts.DEGs <- aggregated.counts.isoLvl[aggregated.counts.isoLvl$gene_id %in% rownames(filtered.pseudobulk.geneLvl), ]
pseudobulk.isoCounts.DEGs$transcript_id <- rownames(pseudobulk.isoCounts.DEGs)
rownames(pseudobulk.isoCounts.DEGs) <- NULL
# Move gene_id and transcript_id to the beginning of the dataframe
pseudobulk.isoCounts.DEGs <- pseudobulk.isoCounts.DEGs[c("gene_id", 
                                                           "transcript_id", 
                                                           setdiff(names(pseudobulk.isoCounts.DEGs), 
                                                                   c("gene_id", "transcript_id")))]

# Reorder the columns of the dataframe
new_column_order <- c("gene_id", "transcript_id", "Neuroepithelial cells", 
                      "Neural Progenitor cells", "Radial glial cells", "Immature neurons")
pseudobulk.isoCounts.DEGs <- pseudobulk.isoCounts.DEGs[, new_column_order]
View(pseudobulk.isoCounts.DEGs)


write.csv(pseudobulk.isoCounts.DEGs, 
          "/data/gpfs/projects/punim0646/manveer/IsoVis_pseudobulk/pseudobulkDEG-isoCounts.csv", 
          row.names = FALSE)





