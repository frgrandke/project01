library(ggplot2)
library(dplyr)
library(data.table)
library(scales)
library(viridis)
library(cowplot)
library(ggpubr)
library(aplot)
library(ggplot2)
library(dplyr)
library(ggh4x)
library(ggplotify)
library(svglite)

colors = fread("data/colors.csv")
color_v = colors$Color
names(color_v) = colors$ID

labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

source("ADRC_theme.R")
source("scripts/helper.R")

color_strips = function(plot, palette, text_color="white") {
  g = ggplot_gtable(ggplot_build(plot))
  
  strips = which(grepl('strip-', g$layout$name))
  
  for (i in seq_along(strips)) {
    k = which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    l = which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    if(length(k) == 0){
      next
    }
    label = g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$label
    g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill = palette[label]
    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col = text_color
  }
  plot_grid(g)
}

fill_title = function(p, palette){
  g <- ggplot_gtable(ggplot_build(p))
  
  strips <- which(grepl('strip-', g$layout$name))
  
  for (i in seq_along(strips)) {
    k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
    label <- g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$label
    if (label == intToUtf8(9792)) {g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette["female"]}
    else if (label == intToUtf8(9794)) {g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette["male"]}
    else if (label %in% labels_cell) {g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette[names(labels_cell)[labels_cell == label]]}
    else {g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- palette[label]}
    g$grobs[[strips[i]]]$grobs[[1]]$children[[l]]$children[[1]]$gp$col <- c( "black", "white")[ 1+(sum( col2rgb(g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill) *c(299, 587,114))/1000 < 123) ]
  }
  return(g)
}


# Figure 4 a)

correlations <- 'results/correlation/celltype_and_celltype_cluster_correlations.csv'
correlations <- fread(correlations, sep=',')
correlations <- correlations[correlations$Diagnosis != "all"]

correlations_cluster = correlations[grepl("celltype_cluster", V1)]
correlations = correlations[grepl("celltype_cluster", V1)]

colors = fread("data/colors.csv", strip.white = F)
color_v = colors$Color
names(color_v) = colors$ID

#correlations_p = correlations[V2 %in% c("Quanterix.CSF___pTau181", "Quanterix3plex.CSF___Ab42",
#                                        "Quanterix4plex.CSF___GFAP",
#                                        "Quanterix4plex.CSF___Tau","Quanterix3plex.CSF___Ab40","Quanterix3plex.CSF___Tau",
#                                        "Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1",
#                                        "Quanterix3plex.CSF___Ab42/Ab40")]

correlations_p = correlations[V2 %in% c(
  "Quanterix4plex.CSF___GFAP",
  "Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1",
  "CSF P-Tau181", "CSF Total Tau", "CSF AB42", "CSF AB40", "CSF AB42/AB40")]

correlations_p = correlations_p[Diagnosis %in% c("Alzheimer's disease", "Healthy Control")]
correlations_p$V1 = unlist(lapply(correlations_p$V1, function(x) substr(x, 18, nchar(x))))
correlations_p[V2 == "Quanterix3plex.CSF___Tau", V2:="Tau (1)"]
correlations_p[V2 == "Quanterix4plex.CSF___Tau", V2:="Tau (2)"]
correlations_p[, V2:=gsub(".*___", "", V2)]
correlations_p[, V2:=gsub("CSF ", "", V2)]
order <- unlist(lapply(unique(correlations_p$V1), function(x) sum(correlations_p$Pearson[correlations_p$V1 == x & correlations_p$Diagnosis == "Alzheimer's disease"])))
correlations_p$V1 <- factor(correlations_p$V1, level = unique(correlations_p$V1)[order(order)])

plotn <- ggplot(correlations_p, aes(x=V1, y=V2 , fill= Pearson))+ 
  geom_tile() +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0, limits=c(-1,1),name="Pearson's\ncorrelation")+
  geom_point(aes(size = ifelse(correlations_p$pearson_p<0.05, "dot", "no_dot"))) +
  scale_size_manual(values = c(dot = 0.8, no_dot = NA), guide = "none")+
  theme_adrc() +
  xlab("") + ylab("") + 
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1))+
  facet_grid(~Diagnosis, scales = "free_y", space = "free_y")+
  scale_x_discrete(name = "", labels = labels_cell[correlations_p$V1])

figure4_a = fill_title(plotn, color_v)
save_plot("results/figures/sig_correlation_csf_marker_celltypes_heatmap_pearson_Alzheimer.pdf", figure4_a, base_width=200, base_height=48, unit="mm")

# Figure 4 b), c) and d)

metadata = fread("results/metadata_with_celltype_cluster_proportions_and_new_biomarker_data.csv")

data <- metadata[,c("Diagnosis", "Plasmacytoid Dendritic cell", "Quanterix4plex.CSF___GFAP",
  "Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1",
  "CSF P-Tau181", "CSF Total Tau", "CSF AB42", "CSF AB40", "CSF AB42/AB40")]
data <- metadata[metadata$Diagnosis %in% c("Healthy Control", "Alzheimer's disease", "Parkinson's Disease only", "Mild Cognitive Impairment"),]

p1 <- ggplot(data, aes(y = Quanterix4plex.CSF___GFAP, x = data[["CD4+ regulatory T cell"]], color = Diagnosis) )+ 
  geom_point()+ 
  scale_color_manual(values = color_v)+ 
  scale_x_continuous(limits = c(0,0.3), name = "Celltpye Proportion: Treg CD4+ T")+ 
  scale_y_continuous(name = "GFAP") +
  theme_adrc() + geom_smooth(method='lm', formula= y~x, se = FALSE)+ theme(legend.position = "none")

p2 <- ggplot(data, aes(y = Quanterix4plex.CSF___NFL, x = data$`Naive CD8+ T cell`, color = Diagnosis) )+ 
  geom_point()+ 
  scale_color_manual(values = color_v)+ 
  scale_y_continuous(limits = c(0,3000), name = "NFL")+ 
  scale_x_continuous(name = "Celltpye Proportion: Naive CD8+ T", limits = c(0,0.08))+ 
  theme_adrc() +
  geom_smooth(method='lm', formula= y~x, se = FALSE)+ theme(legend.position = "none")

p3 <- ggplot(data, aes(y = data$`CSF Total Tau`, x = data$`CD8+ central memory T cell`, color = Diagnosis) )+ 
  geom_point()+ 
  scale_color_manual(values = color_v)+ 
  scale_y_continuous(limits = c(0,2000), name = "Total Tau")+ 
  scale_x_continuous(name = "Celltpye Proportion: Tcm CD8+ T")+ 
  theme_adrc() +
  geom_smooth(method='lm', formula= y~x, se = FALSE)+ theme(legend.position = "none")

figure4_bcd <- plot_grid(p1, p2, p3, nrow = 1, labels = c( 'b', 'c', 'd'))

#ggplot(data, aes(x = Quanterix4plex.CSF___NFL, y = data[["Mucosal associated invariant T cell"]], color = Diagnosis) )+ geom_point()+ scale_color_manual(values = color_v)+ scale_y_continuous(limits = c(0,0.15))+ theme_classic()

#figure <- plot_grid(figure_a, figure_b, figure_c, nrow = 3, labels = c( 'a', 'b'), rel_heights = c(2,1,1))
#pdf("results/figures/Figure_4.pdf")
#figure
#dev.off()
#save_plot("results/figures/Figure_4.pdf", figure, base_width=230, base_height=200, unit="mm", family = "sans", device = cairo_pdf)                                                                               

#save_plot("results/figures/FiguresForPaper/Figure_4a.svg", plot_grid(figure_a, labels = c("a")), base_width=210, base_height=190*0.5, unit="mm")#, family = "sans", device = cairo_pdf)
save_plot("results/figures/FiguresForPaper/Figure_4a.svg", plot_grid(figure_a, labels = c("a")), base_width=210, base_height=190*0.25, unit="mm")#, family = "sans", device = cairo_pdf)
save_plot("results/figures/FiguresForPaper/Figure_4b.svg", figure4_bcd, base_width=210, base_height=190*0.25, unit="mm")#, family = "sans", device = cairo_pdf)                                                                                                                        
