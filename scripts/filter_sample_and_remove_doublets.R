library(Seurat)
library(data.table)
library(scDblFinder)

setDTthreads(snakemake@threads)

mtx = Read10X(snakemake@input$filtered)
seu = CreateSeuratObject(counts=mtx)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

stats = data.table(
  Sample=snakemake@wildcards$sample,
  Barcode=colnames(seu),
  high_mito=seu$percent.mt >= snakemake@params$max_mito_perc,
  low_features=seu$nFeature_RNA <  snakemake@params$min_genes_per_cell,
  high_features=seu$nFeature_RNA > snakemake@params$max_genes_per_cell
)

stats[, discard:=high_mito | low_features | high_features]
stats[, doublet:=NA]

if(sum(!stats$discard) < 10){
  fwrite(stats, snakemake@output$stats, sep='\t')
  saveRDS(seu[, 1], snakemake@output$filtered)
  quit(save="no")
}

seu = subset(seu, percent.mt < snakemake@params$max_mito_perc & nFeature_RNA >= snakemake@params$min_genes_per_cell & nFeature_RNA <= snakemake@params$max_genes_per_cell)

sce = as.SingleCellExperiment(seu)
# expected doublet ratio 1% per 1000 cells
expected_doublet_ratio = (0.01 * ncol(sce) / 1000)
sce = scDblFinder(sce, dbr = expected_doublet_ratio, verbose=FALSE)

stats[match(colnames(sce),stats$Barcode), doublet:=sce$scDblFinder.class != "singlet"]

fwrite(stats, snakemake@output$stats, sep='\t')

seu = seu[,colnames(sce)[sce$scDblFinder.class == "singlet"]]

saveRDS(seu, snakemake@output$filtered)

