library(Seurat)
library(data.table)
library(pbapply)

Sys.setenv(RETICUALTE_PYTHON=system("which python"))
library(reticulate)
snakemake@source("../src/save_as_h5ad.R")

#save.image()

metadata = fread(snakemake@input$metadata)
metadata = metadata[to_exclude_reason %in% c("", "Healthy patient that develops disease")][!Diagnosis %in% c("Lewy Body Disease", "Parkinson's Disease and Dementia")]

samples = snakemake@input$filtered[gsub(".rds", "", basename(snakemake@input$filtered)) %in% metadata$`Date Processed at CG`]
samples = samples[match(gsub(".rds", "", basename(samples)), metadata$`Date Processed at CG`)]

s1 = readRDS(samples[1])

nsamples = length(samples)
print("Parsing samples...")
other_samples = unlist(pblapply(samples[2:nsamples], readRDS, cl=snakemake@threads))

print("Merging samples...")
pbmc = merge(s1, other_samples,
              add.cell.ids=seq(1, nsamples))

####### add patient/sample information ########
print("Add patient/sample information")
pbmc_idx = as.numeric(gsub("_.*", "", colnames(pbmc)))
pbmc$Patient = metadata$SCMD[pbmc_idx]
pbmc$Visit = metadata$Visit[pbmc_idx]
pbmc$Sample = sprintf("%s_Y%s", pbmc$Patient, pbmc$Visit)
pbmc$Sample_batch = metadata$Sample_Batch[pbmc_idx]
pbmc$Processing_batch = metadata$Processing_Batch[pbmc_idx]
pbmc$Diagnosis = metadata$Diagnosis[pbmc_idx]
pbmc$Age = metadata$Age[pbmc_idx]
pbmc$Sex = metadata$Sex[pbmc_idx]
pbmc$Diagnosis_path = metadata$Diagnosis_path[pbmc_idx]

######### add gene IDs ############
print("Add gene IDs")
geneid2symbol = fread(snakemake@input$geneid2symbol)
geneid2symbol = geneid2symbol[!duplicated(GeneSymbol)]
geneid2symbol[grep("_", GeneSymbol), GeneSymbol:=gsub("_", "-", GeneSymbol)]

pbmc@assays$RNA@meta.features$ID_w_vers = geneid2symbol$Geneid[match(rownames(pbmc), geneid2symbol$GeneSymbol)]
pbmc@assays$RNA@meta.features$ID = gsub("\\.\\d+$", "", pbmc@assays$RNA@meta.features$ID_w_vers)
pbmc@assays$RNA@meta.features$Symbol = rownames(pbmc)

print("Saving combined object...")
saveRDS(pbmc, snakemake@output$rds_unfiltered, compress=FALSE)
seurat2anndata(pbmc, snakemake@output$h5ad_unfiltered)
#save.image()
print("Filter gene expression")
gene_expressed_in_cells = pblapply(unique(pbmc$Sample), function(s) {
  apply(pbmc@assays$RNA@counts[,pbmc$Sample == s] > 0, 1, sum)
})

gene_expressed_in_cells_df = setDT(gene_expressed_in_cells)
gene_expressed_in_cells_t = transpose(gene_expressed_in_cells_df)
colnames(gene_expressed_in_cells_t) = rownames(pbmc)
metadata[, merged_diagnosis:=Diagnosis]
metadata[Diagnosis == "Parkinson's Disease and Dementia", merged_diagnosis:="Parkinson's Disease with MCI"]

gene_expressed_in_cells_t[, biogroup:=metadata$merged_diagnosis[match(unique(pbmc$Sample), sprintf("%s_Y%s", metadata$SCMD, metadata$Visit))]]
min_cells_per_genes_per_sample = 3
gene_expressed_in_cells_per_biogroup = gene_expressed_in_cells_t[, pblapply(.SD, function(x) mean(x >= min_cells_per_genes_per_sample)), by=biogroup]
gene_expressed_in_at_least_ratio_biogroup = 0.5
genes_to_keep = colSums(gene_expressed_in_cells_per_biogroup[, 2:ncol(gene_expressed_in_cells_per_biogroup)] >= gene_expressed_in_at_least_ratio_biogroup) >= 1

pbmc = pbmc[names(genes_to_keep[genes_to_keep]),]

print("Saving combined object...")
saveRDS(pbmc, snakemake@output$rds, compress=FALSE)
seurat2anndata(pbmc, snakemake@output$h5ad)
