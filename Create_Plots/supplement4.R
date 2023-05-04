library(ggplot2)
library(cowplot)
library(data.table)
library(ggsci)
library(ggforce)
library(readxl)
library(ggraph)
library(igraph)
library(gghalves)
library(Seurat)
library(reshape2)
library(ggsignif)
library(ggrastr)
library(UpSetR)
library(ComplexUpset)
library(pbapply)
library(ggpubr)
library(ggplotify)
library(ggforce)
library(concaveman)


source("scripts/helper.R")

pbmc = readRDS("results/annotated_celltypes/CompleteObjectAnnotated.rds")


tbl = fread("results/metadata.filtered.csv")
tbl[Diagnosis == "Parkinson's Disease only", Diagnosis:="Parkinson's Disease"]

colors = fread("data/colors.csv")
color_v = colors$Color
names(color_v) = colors$ID


labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID


name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease"="PD", "Parkinson's Disease only"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")


pbmc = SetIdent(pbmc, value="celltype")
pbmc$biogroup_short = name2short[as.character(pbmc$Diagnosis)]
pbmc$biogroup_short = factor(pbmc@meta.data$biogroup_short, levels=diagnosis_order)

sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}


labeller <- function(variable,value){
  return(labels_cell[value])
}

fill_title = function(p, palette){
  g <- ggplot_gtable(ggplot_build(p))
  
  strips <- which(grepl('strip-', g$layout$name))
  
  for (i in seq_along(strips)) {
    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    label <- g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$label
    if (label %in% labels_cell) {g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette[names(labels_cell)[labels_cell == label]]}
    else {g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette[label]}
    # g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "gray"
    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- c( "black", "white")[ 1+(sum( col2rgb(g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill) *c(299, 587,114))/1000 < 123) ]
  }
  return(g)
}

get_cells_per_sample_perc <- function(seur_obj){
  CpS = table(seur_obj$Sample, seur_obj$celltype_cluster)
  CpS_perc = prop.table(CpS, margin=1)
  CpS_perc_df = melt(CpS_perc, varnames = c("Sample", "Celltype"))
  #cells_per_sample_perc_df = cells_per_sample_perc_df[cells_per_sample_perc_df$Celltype != "Unknown/doublets",]
  CpS_perc_df$Diagnosis = factor(seur_obj$biogroup_short[match(CpS_perc_df$Sample, seur_obj$Sample)], levels=diagnosis_order)
  CpS_perc_df$Celltype = factor(CpS_perc_df$Celltype, levels=celltype_cluster_order)
  return (CpS_perc_df)
}

all_pairwise_diag_comparisons = list(c("HC","AD"),c("HC", "PD"), c("HC", "MCI"), c("HC", "PD-MCI"))
print(all_pairwise_diag_comparisons)
#all_pairwise_diag_comparisons = combn(diagnosis_order, 2, simplify = F)
print(all_pairwise_diag_comparisons)

getCelltypeProportion_Boxplot <- function(input_data, all_pairwise_diag_comparisons){
  p_1 = ggplot(input_data, aes(x=Diagnosis, y=value, fill=Diagnosis, color=Diagnosis)) +
    geom_boxplot(outlier.size=0.15, lwd=0.2) + facet_wrap(~Celltype, nrow=1,  labeller=labeller, scales = "free_y") +
    scale_fill_manual(values=color_v) +
    scale_color_manual(values=darken(color_v)) +
    theme_classic() + 
    scale_y_continuous(labels=scales::percent_format()) + xlab("") +
    theme_cowplot(10) + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                              strip.text=element_text(size=7), text = element_text(size=7), axis.text.y = element_text(size = 7), aspect.ratio = 1.2) + ylab("PBMCs")
  p_1_w_sig =  p_1 + 
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="greater"), color="black", tip_length=0) +
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="less"), color="black", tip_length=0) 
  
  return(p_1_w_sig)
}

pbmc_fem <- subset(pbmc, subset = Sex == "female")
cells_per_sample_perc_df_female <- get_cells_per_sample_perc(pbmc_fem)

celltype_sel <- unique(cells_per_sample_perc_df_female$Celltype)

p1 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_female[cells_per_sample_perc_df_female$Celltype %in% celltype_cluster_order[1:10],], all_pairwise_diag_comparisons)
p2 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_female[cells_per_sample_perc_df_female$Celltype %in% celltype_cluster_order[11:20],], all_pairwise_diag_comparisons)
p3 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_female[cells_per_sample_perc_df_female$Celltype %in% celltype_cluster_order[21:30],], all_pairwise_diag_comparisons)
p4 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_female[cells_per_sample_perc_df_female$Celltype %in% celltype_cluster_order[31:39],], all_pairwise_diag_comparisons)


p1 <-  as.ggplot(fill_title(p1, color_v))
p2 <-  as.ggplot(fill_title(p2, color_v))
p3 <-  as.ggplot(fill_title(p3, color_v))
p4 <-  as.ggplot(fill_title(p4, color_v))

p <- plot_grid(p1, p2, p3, p4, nrow = 4)

save_plot("results/figures/supp3_celltype_proportions_per_group_and_sample_female.svg",p, base_height = 140, base_width=280, units="mm")

pbmc_male <- subset(pbmc, subset = Sex == "male")
cells_per_sample_perc_df_male <- get_cells_per_sample_perc(pbmc_male)


p1 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_male[cells_per_sample_perc_df_male$Celltype %in% celltype_cluster_order[1:10],], all_pairwise_diag_comparisons)
p2 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_male[cells_per_sample_perc_df_male$Celltype %in% celltype_cluster_order[11:20],], all_pairwise_diag_comparisons)
p3 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_male[cells_per_sample_perc_df_male$Celltype %in% celltype_cluster_order[21:30],], all_pairwise_diag_comparisons)
p4 <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_male[cells_per_sample_perc_df_male$Celltype %in% celltype_cluster_order[31:39],], all_pairwise_diag_comparisons)

p1 <-  as.ggplot(fill_title(p1, color_v))
p2 <-  as.ggplot(fill_title(p2, color_v))
p3 <-  as.ggplot(fill_title(p3, color_v))
p4 <-  as.ggplot(fill_title(p4, color_v))

p <- plot_grid(p1, p2, p3, p4, nrow = 4)


#p <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_male, all_pairwise_diag_comparisons)
save_plot("results/figures/supp3_celltype_proportions_per_group_and_sample_male.svg", p, base_height = 140, base_width=280, units="mm")

