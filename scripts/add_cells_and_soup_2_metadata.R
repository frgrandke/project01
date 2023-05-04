library(pbapply)
library(data.table)

save.image()

setDTthreads(snakemake@threads)
metadata = fread(snakemake@input$metadata)

num_cells = unlist(pblapply(snakemake@input$filtered,
                            function(x) ncol(readRDS(x)), cl=snakemake@threads))
names(num_cells) = basename(dirname(snakemake@input$filtered))

contamination_stats = rbindlist(lapply(snakemake@input$soupx_stats, function(x) {
  res = as.data.table(jsonlite::fromJSON(x))
  res[, Sample:=gsub(".stats.json$", "", basename(x))]
  return(res)
}))

metadata[, soupx_contamination:=contamination_stats$contamination]
metadata[, num_cells:=num_cells]
metadata[num_cells < snakemake@params$min_cells_per_sample, to_exclude_reason:=paste(to_exclude_reason, "low cell count", sep=';')]
metadata[soupx_contamination == -1 | soupx_contamination > snakemake@params$max_soup, to_exclude_reason:=paste(to_exclude_reason, "high ambient RNA contamination", sep=';')]

fwrite(metadata, snakemake@output$metadata, sep='\t')

