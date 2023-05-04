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
library(grid)
library(viridis)
library(aplot)
library(ggplotify)
library(svglite)
source("ADRC_theme.R")
source("scripts/helper.R")
source("scripts/plotting_functions.R")

fill_title = function(p, palette){
  g <- ggplot_gtable(ggplot_build(p))
  
  strips <- which(grepl('strip-', g$layout$name))
  
  for (i in seq_along(strips)) {
    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette[g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$label]
    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- "white"
  }
  return(g)
}



plot_density_umap = function(embedding, facet) {
  density_estim_per_sample = pblapply(unique(embedding$Sample), function(x, xrange, yrange, h) {
    MASS::kde2d(embedding[Sample == x]$UMAP_1, embedding[Sample == x]$UMAP_2, lims=c(xrange, yrange), h=h, n=200)
  }, xrange=range(embedding$UMAP_1), yrange=range(embedding$UMAP_2),
  h=c(MASS::bandwidth.nrd(embedding$UMAP_1), MASS::bandwidth.nrd(embedding$UMAP_2)))
  names(density_estim_per_sample) = unique(embedding$Sample)
  result_df = rbindlist(pblapply(levels(embedding[[facet]]), function(x){
    densities = density_estim_per_sample[as.character(unique(embedding$Sample[embedding[[facet]] == x]))]
    result = if (length(densities)<1) return(NULL) else densities[[1]]
    result$z = Reduce("+", lapply(densities, function(d) d$z)) / length(densities)
    rownames(result$z) = result$x
    colnames(result$z) = result$y
    result_l = melt(result$z)
    names(result_l) = c("V1", "V2", "z")
    result_l[[facet]] = x
    return(result_l)
  }))
  result_df$z2 = result_df$z
  result_df$z2[result_df$z2<1e-3] = NA
  result_df[[facet]] = factor(result_df[[facet]], levels=levels(embedding[[facet]]))
  ggplot(result_df, aes(V1, V2)) +
    facet_wrap(~biogroup, nrow=1) +
    geom_point_rast(size=0.05, alpha=0.25, data=embedding, aes(x=UMAP_1, y=UMAP_2)) + 
    geom_tile(aes(fill=z2)) +
    stat_contour(aes(z=z, color=..level..), bins=10) +
    scale_fill_viridis_c(na.value="transparent") +
    scale_color_viridis_c(na.value="transparent") +
    coord_cartesian(xlim=range(embedding$UMAP_1), ylim=range(embedding$UMAP_2)) +
    theme_adrc() + xlab("Dimension 1") + ylab("Dimension 2") + guides(colour=FALSE)
}



pbmc = readRDS("results/annotated_celltypes/CompleteObjectAnnotated.rds")
tbl = fread("results/metadata.filtered.csv")
tbl[Diagnosis == "Parkinson's Disease only", Diagnosis:="Parkinson's Disease"]

colors = fread("data/colors.csv", strip.white = F)
color_v = colors$Color
names(color_v) = colors$ID

name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease"="PD", "Parkinson's Disease only"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")

pbmc$biogroup_short = name2short[as.character(pbmc$Diagnosis)]
pbmc$biogroup_short = factor(pbmc@meta.data$biogroup_short, levels=diagnosis_order)



cells_per_sample = table(pbmc$Sample, pbmc$CellCluster)
cells_per_sample_perc = prop.table(cells_per_sample, margin=1)
cells_per_sample_perc_df = melt(cells_per_sample_perc, varnames = c("Sample", "Celltype"))
cells_per_sample_perc_df$Diagnosis = factor(pbmc$biogroup_short[match(cells_per_sample_perc_df$Sample, pbmc$Sample)], levels=diagnosis_order)
cells_per_sample_perc_df$Celltype = factor(cells_per_sample_perc_df$Celltype, levels=names(sort(table(pbmc$CellCluster), decreasing = TRUE)))
sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}

name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease"="PD", "Parkinson's Disease only"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")
tbl[, Diagnosis_short:=name2short[Diagnosis]]
tbl[, Diagnosis_short:=factor(Diagnosis_short, levels=diagnosis_order)]
tbl = tbl[!duplicated(Sample)]
tbl[Sex=="male", Sex:="Male"]
tbl[Sex=="female", Sex:="Female"]

pat_tbl = tbl[RVisit == 1]
plot_df2 = pat_tbl[, list(Frequency=.N), by=c("Diagnosis_short", "Sex")][order(Diagnosis_short, -Sex)]
plot_df2 = plot_df2[, list(Sex, Frequency, pos=cumsum(Frequency)-0.5*Frequency), by=Diagnosis_short]

# Figure 1 b) and c)
fig_1_b = ggplot(plot_df2, aes(x=factor(Diagnosis_short, levels = diagnosis_order), y = Frequency)) + geom_bar(aes(fill=Sex), stat="identity", alpha = 0.8) +
  geom_text(aes(label=Frequency, y=pos), size=2) +
  scale_y_continuous(expand = expansion(mult=c(0, 0))) +
  theme_adrc() + scale_fill_manual(values=color_v) + ylab("#Patients") + xlab("") +
  theme(legend.position="none", legend.direction = "vertical") #, strip.text=element_text(size=7))

fig_1_c = ggplot(pat_tbl, aes(x=factor(Diagnosis_short, levels = diagnosis_order), y = Age)) +
  geom_half_violin(aes(y=Age), data=pat_tbl[Sex == "Female"], side="l", color=darken(color_v["Female"]), fill=color_v["Female"], alpha = 0.8) + #+ geom_half_boxplot(aes(fill=Sex), side="l") + ggbeeswarm::geom_beeswarm(aes(fill=Sex)) +
  geom_half_boxplot(aes(y=Age), data=pat_tbl[Sex == "Female"], side="l", fill=color_v["Female"], width=0.2, alpha = 0.8) +
  geom_half_violin(aes(y=Age), data=pat_tbl[Sex == "Male"], side="r", color=darken(color_v["Male"]), fill=color_v["Male"], alpha = 0.8) +
  geom_half_boxplot(aes(y=Age), data=pat_tbl[Sex == "Male"], side="r", fill=color_v["Male"], width=0.2, alpha = 0.8) +
  theme_adrc() + ylab("Age (years)") + xlab("") + theme(legend.position="none")


figure1_bc = as.ggplot(plot_grid(fig_1_b + theme(plot.margin=margin(t=5)),
                    fig_1_c + theme(axis.text.x = element_blank(),
                                    plot.margin = margin(),
                                    axis.ticks.x = element_blank()), ncol=1, align="hv", labels = c("b", "c")))

pbmc = SetIdent(pbmc, value="celltype")
figure1_d = DimPlot(pbmc, reduction="umap", label=F) + scale_color_manual(values=color_v) + NoLegend() + scale_x_continuous(expand=expand_scale(mult = c(0.025,0.025)))+theme_void() + theme(aspect.ratio = 1)+ theme(legend.position = "none")
#save_plot("results/figures/figure1_umap_celltype_cluster.png", fill_title(celltype_umap, color_v), base_height = 60, base_width=60, units="mm")


# Figure 1e)
meta.data <- pbmc@meta.data

l1 <- data.frame(table(meta.data$Manual_Annotation))
l2 <- data.frame(table(meta.data$Group))
l3 <- data.frame(table(meta.data$CellCluster))
vertices <- data.frame(name = c(factor(""),l1$Var1, l2$Var1, l3$Var1), size = (c(length(meta.data$Manual_Annotation), l1$Freq, l2$Freq, l3$Freq)))
vertices <- vertices[!duplicated(vertices$name),]

e1 <- meta.data[, colnames(meta.data) %in% c("Group","Manual_Annotation")]
e1 <- e1[!duplicated(e1),]
e2 <- meta.data[, colnames(meta.data) %in% c("CellCluster", "Group")]
e2 <- e2[!duplicated(e2),]
edges <- data.frame(from = c(rep("", length(unique(meta.data$CellCluster))),e2$CellCluster, e1$Group), to = c(unique(meta.data$CellCluster),e2$Group,e1$Manual_Annotation)) ## , #
edges <- edges[!(edges$from == edges$to),]
edges <- edges[!duplicated(edges),]
mygraph <- graph_from_data_frame( edges, vertices=vertices )
edges$from[!edges$to %in% vertices$name]

figure1_e <- ggraph(mygraph, layout = 'circlepack', weight=size) +
  geom_node_circle() +
  theme_void()+
  geom_node_circle(aes(fill = depth)) +
  geom_node_text( aes(label=name, filter=leaf, fill=depth, size=size)) +
  scale_fill_viridis() + theme(legend.position="none", strip.text=element_text(size=7))

#ggsave("results/figures/Overview_CirclePlot.pdf", plot, height = 120, width=120, units="mm")


#Figure 1 f)

umap_embedding_df = as.data.table(pbmc@reductions$umap@cell.embeddings)
umap_embedding_df$biogroup = factor(pbmc$biogroup_short, levels=diagnosis_order)
print(head(umap_embedding_df))
umap_embedding_df$celltype = factor(pbmc$CellCluster, levels=celltype_order)
umap_embedding_df$Sample = pbmc$Sample

density_umap_per_biogroup = ggplot(umap_embedding_df, aes(x=UMAP_1, y=UMAP_2)) +
  geom_point_rast(size=0.05, alpha=0.25) +
  stat_density_2d(aes(fill=stat(nlevel)), geom="polygon", n=200, size=0.5) +
  facet_grid(. ~ biogroup) + scale_fill_viridis_c() +
  theme_adrc() + xlab("") + ylab("") +
  theme(legend.position = "none") + theme(legend.position="none",
                                          #strip.text=element_text(size=7),
                                          #text = element_text(size=7),
                                          axis.text.x = element_blank(), axis.text.y = element_blank(),aspect.ratio = 1, axis.line=element_blank(), axis.ticks=element_blank())

figure1_f <- as.ggplot(fill_title(density_umap_per_biogroup, color_v))

#save_plot("results/figures/figure1_density_umap_per_biogroup.pdf", fill_title(density_umap_per_biogroup, color_v), base_height = 50, base_width=120, units="mm")


save_plot("results/figures/FiguresForPaper/figure1_bc.svg", plot_grid(figure1_bc,labels = c("")) , base_width=210, base_height=230*0.5, unit = "mm")
save_plot("results/figures/FiguresForPaper/figure1_f.svg", plot_grid(figure1_f,labels = c("f")) , base_width=0.5*210, base_height=230*0.2, unit = "mm")
save_plot("results/figures/FiguresForPaper/figure1_d.svg", plot_grid(figure1_d,labels = c("d")) , base_width=0.6*210, base_height=230*0.5, unit = "mm")
save_plot("results/figures/FiguresForPaper/figure1_e.svg", plot_grid(figure1_e,labels = c("e")) , base_width=0.6*210, base_height=230*0.4, unit = "mm")


### Figure 2 ###################################################################

sigFunc = function(x){
  if(x < 0.001){"***"}
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}

all_pairwise_diag_comparisons = combn(diagnosis_order, 2, simplify = F)

color_v = c(color_v, "CD16+ Mono."=unname(color_v["CD16+ Monocyte"]), "CD14+ Mono."=unname(color_v["CD14+ Monocyte"]),
            "Imm. B cell"=unname(color_v["Immature B cell"]), "Mega."=unname(color_v["Megacaryocyte"]))

celltypes1 <- c("CD4+ T cell")
celltypes2 <-  c("CD8+ T cell", "NK cell", "B cell", "Monocyte", "Gamma delta T cell", "Mucosal associated invariant T cell")
celltypes3 <-  c("Dendritic cell", "Megacaryocytes", "NKT-like cell", "Plasma cell")

CellCluster_order <- c("CD4+ T cell", "CD8+ T cell", "NK cell", "B cell", "Monocyte", "Dendritic cell", "Gamma delta T cell", "Megacaryocytes", "Mucosal associated invariant T cell", "NKT-like cell", "Plasma cell")

get_cells_per_sample_perc <- function(seur_obj){
  CpS = table(seur_obj$Sample, seur_obj$CellCluster)
  CpS_perc = prop.table(CpS, margin=1)
  CpS_perc_df = melt(CpS_perc, varnames = c("Sample", "Celltype"))
  CpS_perc_df$Diagnosis = factor(seur_obj$biogroup_short[match(CpS_perc_df$Sample, seur_obj$Sample)], levels=diagnosis_order)
  CpS_perc_df$Celltype = factor(CpS_perc_df$Celltype, levels=CellCluster_order)
  return (CpS_perc_df)
}

pbmc_fem <- subset(pbmc, subset = Sex == "female")
cells_per_sample_perc_df_female <- get_cells_per_sample_perc(pbmc_fem)
figure2_a_fem <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_female, celltypes1, celltypes2, all_pairwise_diag_comparisons)

pbmc_male <- subset(pbmc, subset = Sex == "male")
cells_per_sample_perc_df_male <- get_cells_per_sample_perc(pbmc_male)
figure2_a_male <- getCelltypeProportion_Boxplot(cells_per_sample_perc_df_male, celltypes1, celltypes2, all_pairwise_diag_comparisons)

#### Figure 2b)

umap_data = data.table(pbmc@reductions$umap@cell.embeddings, biogroup=pbmc$biogroup_short, celltype=pbmc$celltype, celltype_cluster=pbmc$CellCluster, sample=pbmc$Sample, Sex = pbmc$Sex)
setnames(umap_data, c("UMAP_1", "UMAP_2"), c("V1", "V2"))

figure2_b_MCIvsHC_f = plot_density_difference(umap_data[biogroup=="MCI" & Sex == "female"], umap_data[biogroup=="HC" & Sex == "female"])
figure2_b_ADvsHC_f = plot_density_difference(umap_data[biogroup=="AD"& Sex == "female"], umap_data[biogroup=="HC" & Sex == "female"])
figure2_b_PDvsHC_f = plot_density_difference(umap_data[biogroup=="PD"& Sex == "female"], umap_data[biogroup=="HC" & Sex == "female"])
figure2_b_MCIvsHC_m = plot_density_difference(umap_data[biogroup=="MCI" & Sex == "male"], umap_data[biogroup=="HC" & Sex == "male"])
figure2_b_ADvsHC_m = plot_density_difference(umap_data[biogroup=="AD"& Sex == "male"], umap_data[biogroup=="HC" & Sex == "male"])
figure2_b_PDvsHC_m = plot_density_difference(umap_data[biogroup=="PD"& Sex == "male"], umap_data[biogroup=="HC" & Sex == "male"])

legend <- cowplot::get_legend(figure2_b_PDvsHC_m)
pdf("fig2_umap_density_legend.pdf")
grid.newpage()
grid.draw(legend)
dev.off()

# Figure 2 c)

meta.data <- pbmc@meta.data

labels = fread("data/CelltypeMapping.csv", strip.white = F)
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

# 
# c <- seq(nrow(meta.data[, colnames(meta.data) %in% c("celltype_cluster", "Diagnosis", "Sex")]))
# counts <- aggregate(c~., data=meta.data[,colnames(meta.data) %in% c("celltype_cluster", "Diagnosis", "Sex")], FUN=length)
# 
# 
# counts$CellAbbr <- labels_cell[counts$celltype_cluster]
# counts$CellAbbr <- factor(counts$CellAbbr, levels = labels_cell[celltype_cluster_order])
# counts$celltype_cluster <- factor(counts$celltype_cluster, levels = celltype_cluster_order)
# celltype_cluster_order[!celltype_cluster_order %in% unique(counts$celltype_cluster) ]
# counts <- counts[counts$Diagnosis != "",]
# 
# tmp = unlist(lapply(c(1:length(counts$c)), function(i) {
#   cellt <- counts$celltype_cluster[i]
#   gender <- counts$Sex[i]
#   diagn <- counts$Diagnosis[i]
#   (counts$c[i]/sum(counts$c[counts$Sex == gender& counts$Diagnosis == diagn]))
# }))
# counts$c <- tmp
# 
# tmp = unlist(lapply(c(1:length(counts$c)), function(i) {
#   cellt <- counts$celltype_cluster[i]
#   gender <- counts$Sex[i]
#   diagn <- counts$Diagnosis[i]
#   ((counts$c[i]) - (counts$c[counts$celltype_cluster == cellt & counts$Sex == gender & counts$Diagnosis == "Healthy Control"]))[1]
# }))
# counts$c <- tmp
# 
# counts[is.na(counts$celltype_cluster),]
# 
# counts <- counts[counts$Diagnosis != "Healthy Control",]
# print(head(counts$c))
# 
# print(max(counts$c))
# plot <- ggplot(counts[counts$Sex == "male",], aes(x = celltype_cluster, y = Diagnosis))+ 
#   geom_tile(aes(fill = c)) +  
#   scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", limits = c(-0.11, 0.11), name = "Abs. Change in the %Cell-type proportion") +
#   theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(Sex~.)+
#   scale_x_discrete(name = "") + 
#   scale_y_discrete(name = "") + labs(fill='Abs. Change in %Cells') + scale_x_discrete(name = "", labels = labels_cell[counts$celltype_cluster])
# 
# plot2 <- ggplot(counts[counts$Sex == "female",], aes(x = celltype_cluster, y = Diagnosis))+ 
#   geom_tile(aes(fill = c)) +  
#   scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", limits = c(-0.11, 0.11)) +
#   theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(Sex~.)+
#   scale_x_discrete(name = "") + 
#   scale_y_discrete(name = "") + labs(fill='Abs. Change in %Cells') +
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# 
# 
# 
# labels= ggplot(counts, aes(x=celltype_cluster, y=1, fill=celltype_cluster)) + geom_tile() +
#   scale_fill_manual(values = color_v) +
#   theme_void() + guides(fill="none")
# 
# p <- plot2 %>% insert_top(labels, height=.25) %>% insert_bottom(plot, height = 1.0)
# ggsave("results/figures/NumCells_Heatmap.pdf", p, width = 9, height = 3)



c <- seq(nrow(meta.data[, colnames(meta.data) %in% c("celltype_cluster", "Diagnosis", "Sex")]))
counts <- aggregate(c~., data=meta.data[,colnames(meta.data) %in% c("celltype_cluster", "Diagnosis", "Sex")], FUN=length)


counts$CellAbbr <- labels_cell[counts$celltype_cluster]
counts$CellAbbr <- factor(counts$CellAbbr, levels = labels_cell[celltype_cluster_order])
counts$celltype_cluster <- factor(counts$celltype_cluster, levels = celltype_cluster_order)
celltype_cluster_order[!celltype_cluster_order %in% unique(counts$celltype_cluster) ]
counts <- counts[counts$Diagnosis != "",]


tmp = unlist(lapply(c(1:length(counts$c)), function(i) {
  cellt <- counts$celltype_cluster[i]
  gender <- counts$Sex[i]
  diagn <- counts$Diagnosis[i]
  (counts$c[i]/sum(counts$c[counts$Sex == gender& counts$Diagnosis == diagn]))
}))
counts$c <- tmp

tmp = unlist(lapply(c(1:length(counts$c)), function(i) {
  cellt <- counts$celltype_cluster[i]
  gender <- counts$Sex[i]
  diagn <- counts$Diagnosis[i]
  print((counts$c[i]))
  print((counts$c[counts$celltype_cluster == cellt & counts$Sex == gender & counts$Diagnosis == "Healthy Control"]))
  ((counts$c[i]) / (counts$c[counts$celltype_cluster == cellt & counts$Sex == gender & counts$Diagnosis == "Healthy Control"]))[1]
}))
counts$c <- tmp

counts[is.na(counts$celltype_cluster),]

counts <- counts[counts$Diagnosis != "Healthy Control",]

plot <- ggplot(counts[counts$Sex == "male",], aes(x = celltype_cluster, y = Diagnosis))+ 
  geom_tile(aes(fill = c)) +  
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",midpoint = 1,  name = "Rel. Change in the %Cell-type proportion") +
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(Sex~.)+
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "") + labs(fill='Abs. Change in %Cells') + scale_x_discrete(name = "", labels = labels_cell[counts$celltype_cluster]) #+

plot2 <- ggplot(counts[counts$Sex == "female",], aes(x = celltype_cluster, y = Diagnosis))+ 
  geom_tile(aes(fill = c)) +  
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027",midpoint = 1) + 
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(Sex~.)+
  scale_x_discrete(name = "") +
  scale_y_discrete(name = "") + labs(fill='Abs. Change in %Cells') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), legend.position = "bottom")

labels= ggplot(counts, aes(x=celltype_cluster, y=1, fill=celltype_cluster)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill="none")

figure2_c <- plot2 %>% insert_top(labels, height=.25) %>% insert_bottom(plot, height = 1.0)

save_plot("results/figures/figure1_celltype_proportions_per_group_and_sample_female.pdf", figure2_a_fem, base_height = 40, base_width=220, units="mm")
save_plot("results/figures/figure1_celltype_proportions_per_group_and_sample_male.pdf", figure2_a_male, base_height = 40, base_width=220, units="mm")

save_plot("results/figures/fig2_umap_density_diff_hc_vs_mci_female.pdf", figure2_b_MCIvsHC_f, base_height = 50, base_width = 50, unit = "mm")
save_plot("results/figures/fig2_umap_density_diff_hc_vs_ad_female.pdf", figure2_b_ADvsHC_f, base_height = 50, base_width = 50, unit = "mm")
save_plot("results/figures/fig2_umap_density_diff_hc_vs_pd_female.pdf", figure2_b_PDvsHC_f, base_height = 50, base_width = 50, unit = "mm")
save_plot("results/figures/fig2_umap_density_diff_hc_vs_mci_male.pdf", figure2_b_MCIvsHC_m, base_height = 50, base_width = 50, unit = "mm")
save_plot("results/figures/fig2_umap_density_diff_hc_vs_ad_male.pdf", figure2_b_ADvsHC_m, base_height = 50, base_width = 50, unit = "mm")
save_plot("results/figures/fig2_umap_density_diff_hc_vs_pd_male.pdf", figure2_b_PDvsHC_m, base_height = 50, base_width = 50, unit = "mm")

pdf("results/figures/NumCells_Heatmap_2.pdf")
plot
dev.off()

