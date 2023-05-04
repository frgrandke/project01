colors_full = c("Alzheimer's disease"="#ff7f00",
            "Healthy Control"="#33a02c",
            "Mild Cognitive Impairment"="#e31a1c",
            "Parkinson's Disease with MCI"="#91003f",
            "Parkinson's Disease"="#6a3d9a",
            "All"="#ffffff",
            "Neurodegeneration"="#1f78b4",
            "Cognitive Impairment"="#a6cee3"
)

colors_short = c("AD"="#ff7f00",
                 "HC"="#33a02c",
                 "MCI"="#e31a1c",
                 "PD & MCI"="#91003f",
                 "PD"="#6a3d9a",
                 "All"="#ffffff",
                 "Neurodegeneration"="#1f78b4",
                 "Cognitive Impairment"="#a6cee3"
)

darken <- function(color, factor=1.2){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  names(col) = names(color)
  col
}


diagnosis_order = c("HC", "MCI", "AD", "PD-MCI", "PD")
celltype_order = c("CD4+ T-Helper Cell","Naive CD4+ T cell","CD4+ Memory T cell",
                   "Gamma delta T cell","Mucosal associated invariant T cell",
                   "CD8+ Memory T cell","Naive CD8+ T cell","Proliferating CD8+ T cell",
                   "Memory B cell","Double negative B cell","Naive B cell",
                   "B cell","Transitional B cell","Conventional Dendritic cell",
                   "CD16+ Monocyte","Megacaryocytes", "Proliferating Monocyte",
                   "CD14+ Monocyte","Plasmacytoid Dendritic cell","HSPC",
                   "Red Blood Cells","NKT-like cell","Proliferating CD4+ T cell",
                   "Plasmablasts","CD56-Dim NK cell", "CD56-Bright NK cell")
                   
CellCluster_order <- c("CD4+ T cell", "CD8+ T cell", "NK cell", "NKT-like cell", "Gamma delta T cell", "Mucosal associated invariant T cell", "B cell", "Plasma cell", "Monocyte", "Megacaryocytes", "Dendritic cell")


celltype_cluster_order = c("Naive CD4+ T cells", "CD4+ regulatory T cell", "Follicular helper CD4 T cell",
"CD4+ T-Helper 1 Cell", "CD4+ T-Helper 2 Cell", "CD4+ T-Helper 17 Cell", 
"CD4+ central memory T cell", "CD4+ transitional memory T cell", "CD4+ effector memory T cell",
"Terminal effector CD4+ T cell", "Proliferating CD4+ T cell ", "Naive CD8+ T cell", 
"CD8+ stem cell memory T cell", "CD8+ central memory T cell", "CD8+ effector memory T cell", 
"CD8+ Terminal effector T cell", "Proliferating CD8+ T cell", "Gamma delta T cell",
"Mucosal associated invariant T cell", "NKT-like cell", "CD56- Natural Killer cell",
"CD56+ Natural Killer cell", "Naive B cell", "Classical Memory B cell", "IgM Memory B cell",
"B cell", "Transitional B cell", "Double negative 1 B cell", "Double negative 2 B cell",
"Double negative 3 B cell", "Double negative 4 B cell", "Plasmablasts", "Plasmacytoid Dendritic cell", 
"Conventional Dendritic cell Type 1", "Conventional Dendritic cell Type 2", 
"CD14+ Monocyte", "CD16+ Monocyte", "Proliferating Monocytes", "Megacaryocytes")

name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease only"="PD",
               "Parkinson's Disease"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")

get_celltype_colors = function(celltypes) {
  colors = rep("", length(celltypes))
  colors[grepl("^NK", celltypes)] = "#D62728FF"
  colors[grepl("T cell", celltypes)] = "#FF7F0EFF"
  
  colors[grepl("CD4 Memory", celltypes)] = "#fdd0a2"
  colors[grepl("CD4 T Helper2", celltypes)] = "#fd8d3c"
  colors[grepl("CD4 T Reg", celltypes)] = "#d94801"
  colors[grepl("CD4 Naive T", celltypes)] = "#7f2704"
  colors[grepl("CD8 Naive Cytotoxic", celltypes)] = "#980043"
  colors[grepl("CD8 Cytotoxic T", celltypes)] = "#ce1256"
  
  colors[grepl("CD4 T cell", celltypes)] = "#993404"
  colors[grepl("CD8 T cell", celltypes)] = "#bd0026"
  colors[grepl("CD4 CD3D", celltypes)] = "#993404"
  colors[grepl("CD8 CD3D", celltypes)] = "#bd0026"
  
  colors[grepl("Plasma", celltypes)] = "#1F77B4FF"
  
  colors[grepl("^B", celltypes)] = "#9467BDFF"
  colors[grepl("Mature B", celltypes)] = "#603a83"
  colors[grepl("Immature B", celltypes)] = "#ac8acc"
  
  colors[grepl("Monocyte", celltypes)] = "#a6d96a"
  colors[grepl("Mature Monocyte", celltypes)] = "#3d5d18"
  colors[grepl("CD16+ Monocyte", celltypes)] = "#d5edba"
  
  colors[grepl("Granulocyte", celltypes)] = "#1b6e4b"
  
  colors[grepl("Dendritic", celltypes)] = "#1a9850"
  
  colors[grepl("Myeloid cell", celltypes)] = "#2CA02CFF"
  names(colors) = celltypes
  return(colors)
}



make_dot_plot = function(pbmc, biogroup, celltype, genes, zscore=TRUE) {
  plot_df = rbindlist(lapply(unique(pbmc@meta.data[[biogroup]]), function(d) {
    sub_expr = pbmc@assays$RNA@data[genes, pbmc@meta.data[[biogroup]] == d]
    celltypes = pbmc@meta.data[[celltype]][pbmc@meta.data[[biogroup]] == d]
    result = data.frame(Gene=character(), Celltype=character(), avg_expr=numeric(), perc_expressed=numeric())
    for(c in unique(pbmc@meta.data[[celltype]])){
      percent_expressed = rowMeans(sub_expr[genes, celltypes==c] > 0) * 100
      avg = rowMeans(sub_expr[genes, celltypes==c])
      for(i in 1:length(percent_expressed)){
        g = names(percent_expressed)[i]
        result = rbind(result, data.frame(Gene=g, Celltype=c, avg_expr=avg[i], perc_expressed=percent_expressed[i]))
      }
    }
    result$Diagnosis = d
    return(as.data.table(result))
  }))
  #print(plot_df)
  
  if(zscore){
    plot_df[, scaled_expr:=(avg_expr-mean(avg_expr))/sd(avg_expr), by=c("Gene", "Diagnosis")]
  } else {
    plot_df[, scaled_expr:=avg_expr]
  }
  if(any(diagnosis_order %in% plot_df$Diagnosis)) {
    plot_df[, Diagnosis:=factor(Diagnosis, diagnosis_order)]
  }
  if(any(diagnosis_order %in% plot_df$Celltype)) {
    plot_df[, Celltype:=factor(Celltype, diagnosis_order)]
  }
  if(any(celltype_order %in% plot_df$Celltype)) {
    plot_df[, Celltype:=factor(Celltype, celltype_order)]
  }
  if(any(celltype_order %in% plot_df$Diagnosis)) {
    plot_df[, Diagnosis:=factor(Diagnosis, celltype_order)]
  }
  
  
  # specific plot
  color_range = c(-max(abs(min(plot_df$scaled_expr)), abs(max(plot_df$scaled_expr))))
  color_range[2] = -color_range[1]
  
  expr_range = c(0, max(plot_df$perc_expressed))
  
  p1 = ggplot(plot_df[order(Diagnosis)], aes(x=Celltype, y=sprintf("%13s", Gene), size=perc_expressed, color=scaled_expr)) +
    facet_grid(rows=vars(Diagnosis), scales = "free_y", space="free", switch="y") +
    geom_point() +
    scale_color_distiller(palette="Blues", name="Z-score", limits=color_range) + scale_size(name="Expressed (%)", limits=expr_range, range=c(0.5,4)) + 
    scale_x_discrete(guide = guide_axis(angle = 90)) +
    theme_cowplot(8) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          legend.position = "none",
          legend.box="vertical",
          legend.margin=margin(),
          strip.text.y.left=element_text(angle=0, size=8),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.placement = "outside",
          panel.spacing.y = unit(0.25, "lines"))
  
  color_strips = function(p, colors) {
    g <- ggplot_gtable(ggplot_build(p))
    
    strips <- which(grepl('strip-', g$layout$name))
    
    for (i in seq_along(strips)) {
      k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
      l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
      g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- colors[i]
      g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white"
    }
    
    return(g)
  }
  
  p1 = color_strips(p1, color_v[levels(plot_df$Diagnosis)])

  plot_grid(p1)
}
