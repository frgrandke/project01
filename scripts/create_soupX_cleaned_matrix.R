library(SoupX)
library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)
library(jsonlite)
library(RhpcBLASctl)

blas_p = blas_get_num_procs()
blas_set_num_threads(1)
omp_p = omp_get_num_procs()
omp_set_num_threads(1)

#save.image()
read_pisa = function(raw_dir){
  counts = Seurat::Read10X(raw_dir, gene.column = 1)
  CreateSeuratObject(counts=counts)
}

read_pisa_filtered = function(filtered_mtx) {
  counts = fread(filtered_mtx, sep='\t')
  mat = as.matrix(counts[, 2:ncol(counts)])
  rownames(mat) = counts$ID
  CreateSeuratObject(counts=mat)
}

create_adjusted_counts = function(raw, filtered) {
  filtered_counts = raw[, intersect(colnames(filtered), colnames(raw))]
  # filtering to remove cells that were incorrectly included when selecting the cell calling threshold
  # if no cell remains just output a matrix with one cell so that we discard this sample later
  if(all(filtered_counts$nFeature_RNA < snakemake@params$min_genes_per_cell)){
    return(list(out=filtered_counts@assays$RNA@counts[,1:2], contamination=-1, clusters=0,
                original_cells=ncol(filtered), common_cells=0))
  }
  
  filtered_counts = subset(filtered_counts, nFeature_RNA >= snakemake@params$min_genes_per_cell)
  if(ncol(filtered_counts) > 30000){
    stop(sprintf("More than 30k cells (%d) for samples %s", ncol(filtered_counts), snakemake@input$filtered))
  } else if(ncol(filtered_counts) < 300) {
    return(list(out=filtered_counts@assays$RNA@counts, contamination=-1, clusters=0,
                original_cells=ncol(filtered), common_cells=ncol(filtered_counts)))
  }
  
  filtered_counts = NormalizeData(filtered_counts, verbose = FALSE)
  filtered_counts = FindVariableFeatures(filtered_counts, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  filtered_counts = ScaleData(filtered_counts, verbose = FALSE)
  filtered_counts = RunPCA(filtered_counts, verbose=FALSE)
  
  #ElbowPlot(filtered_counts)
  
  # 12 PCs seem to be fine for our data
  filtered_counts = FindNeighbors(filtered_counts, dims = 1:12, verbose=FALSE)
  
  # rather do over-clustering than under clustering
  resolution = 0.6
  
  filtered_counts = FindClusters(filtered_counts, resolution = resolution, verbose = FALSE)
  count = 1
  
  # try to make at least 10 clusters
  while(length(levels(filtered_counts$seurat_clusters)) < 9 && count < 30){
    resolution = resolution + 0.05
    filtered_counts = FindClusters(filtered_counts, resolution = resolution, verbose = FALSE)
    count = count + 1
  }
  
  
  filtered_counts = RunUMAP(filtered_counts, dims = 1:12, verbose = FALSE)

  sc = SoupChannel(raw@assays$RNA@counts, filtered_counts@assays$RNA@counts)
  
  sc = setDR(sc, filtered_counts@reductions$umap@cell.embeddings)
  sc = setClusters(sc, filtered_counts$seurat_clusters)
  if(!is.null(snakemake@params$cont_fraction)) {
    sc = setContaminationFraction(sc, snakemake@params$cont_fraction, forceAccept = TRUE)
  } else {
    pdf(sprintf("%s.pdf", snakemake@output$counts))
    sc = autoEstCont(sc, forceAccept = TRUE)
    dev.off()
    if(sc$metaData$rho[1] > 0.5){
      stop()
    }
  }
  
  out = adjustCounts(sc)

  return(list(out=out, contamination=sc$metaData$rho[1], clusters=length(levels(sc$metaData$clusters)),
              original_cells=ncol(filtered), common_cells=ncol(filtered_counts)))
}

#raw_dir = "data/unfiltered_with_new/20200218-1"
raw_dir = snakemake@input$raw
#filtered_mtx = "data/filtered_with_new/20200218-1/counts.tsv.gz"
filtered_mtx = snakemake@input$filtered


raw = read_pisa(raw_dir)
filtered = read_pisa_filtered(filtered_mtx)

result = create_adjusted_counts(raw, filtered)
DropletUtils::write10xCounts(snakemake@output$counts, result$out)
result$out = NULL
write_json(result, snakemake@output$stats)

