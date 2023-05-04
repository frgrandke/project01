library(ggplot2)
library(ggrepel)
library(cowplot)
library(ggplot2)
library("ggpubr")
library(dplyr)
library(data.table)
library(aplot)
#library(ComplexHeatmap)
source("scripts/helper.R")
library(viridis)
library(ggplotify)

colors = fread("data/colors.csv")
color_v = colors$Color
names(color_v) = colors$ID

labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

source("ADRC_theme.R")

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


degs_male <- read.csv("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_male <- degs_male[degs_male$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"),]
degs_male$Sex = "male"
degs_female <- read.csv("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_female$Sex = "female"
degs_female <- degs_female[degs_female$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"),]

degs_male_t <- degs_male[(degs_male$gene %in% degs_female$gene) ,]
degs_female_t <- degs_female[(degs_female$gene %in% degs_male$gene),]

degs <- degs_male_t %>% right_join(degs_female_t, by= c("gene", "cluster_id", "contrast"))

Correlation <- NA
for (cellt in unique(degs$cluster_id)){
  for (cont in unique(degs$contrast[degs$cluster_id == cellt])){
    tryCatch(
      {
        corr <- cor.test(degs$logFC.x[degs$contrast == cont & degs$cluster_id == cellt], degs$logFC.y[degs$contrast == cont & degs$cluster_id == cellt], method=c("pearson", "kendall", "spearman"))
        res <- data.frame(Value = corr$estimate, p_val = corr$p.value, Celltype = cellt, Comparison = cont)
        if (any(is.na(Correlation))) {Correlation <- res} else {Correlation <- rbind(Correlation, res)}
      },
      error=function(cond) {
        print(cond)
      })
  }
}


# Figure 3 a)

degs_male <- degs_male[degs_male$p_adj.loc < 0.05 & abs(degs_male$logFC)>0.5,]
degs_female <- degs_female[degs_female$p_adj.loc < 0.05 & abs(degs_female$logFC)>0.5,]

degs_number <- rbind(degs_male, degs_female)

head(as.data.frame(table(comb = do.call(paste, degs_number[colnames(degs_number) %in% c("cluster_id", "contrast"),])), responseName = "n"))

c <- seq(nrow(degs_number[, colnames(degs_number) %in% c("cluster_id", "contrast", "Sex")]))
counts <- aggregate(c~., data=degs_number[,colnames(degs_number) %in% c("cluster_id", "contrast", "Sex")], FUN=length)


counts$contrast <- factor(counts$contrast, levels = c( "PDvsHC","MCIvsHC","ADvsHC"))

counts$CellAbbr <- labels_cell[counts$cluster_id]
counts$CellAbbr <- factor(counts$CellAbbr, levels = labels_cell[celltype_cluster_order])
counts$cluster_id <- factor(counts$cluster_id, levels = celltype_cluster_order)

labels= ggplot(counts, aes(x=cluster_id, y=1, fill=cluster_id)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill="none")
counts$c[counts$c == 0] <- NA

plot <- ggplot(counts[counts$Sex == "male",], aes(x = cluster_id, y = contrast))+ 
  geom_tile(aes(fill = c)) +  
  scale_fill_gradient2(low = "white", high = "#5D7322", trans = "log", limits = c(1,1000), breaks = c(1,5,50,500)) +
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(Sex~.)+
  scale_x_discrete(name = "", labels = labels_cell[counts$cluster_id]) + scale_y_discrete(name = "") + labs(fill='#DEGs ') + theme(legend.position="bottom") + theme(legend.title=element_text(size=8), legend.text=element_text(size=6))#+

plot2 <- ggplot(counts[counts$Sex == "female",], aes(x = cluster_id, y = contrast))+ 
  geom_tile(aes(fill = c)) +  
  scale_fill_gradient2(low = "white", high = "#9C489A", trans = "log", limits = c(1,1000), breaks = c(1,5,50,500)) +
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ facet_grid(Sex~.)+
  scale_x_discrete(name = "", labels = labels_cell[counts$cluster_id]) + scale_y_discrete(name = "") + labs(fill='#DEGs ') +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + theme(legend.position="bottom")+ theme(legend.title=element_text(size=8), legend.text=element_text(size=6))

legend1 <- cowplot::get_legend(plot)
#ggsave("results/figures/NumDEGs_Male_Female_l1.pdf", legend1, width = 2, height = 1)
legend2 <- cowplot::get_legend(plot2)
#ggsave("results/figures/NumDEGs_Male_Female_l2.pdf", legend2, width = 2, height = 1)

plot <- plot + theme(legend.position="none")
plot2 <- plot2 + theme(legend.position="none")

p <- plot2 %>% insert_top(labels, height=.25) %>% insert_bottom(plot, height = 1.0)
l <- as.ggplot(plot_grid(legend1, legend2, nrow = 1))
l <- l + theme(plot.margin = unit(c(0, 0, 1, 2), "cm"))
figure3_a <- plot_grid(as.ggplot(p), l, nrow = 2, rel_heights = c(1,0.1))

#ggsave("results/figures/NumDEGs_Male_Female.pdf", figure3_a, width = 7, height = 3)

# Figure 3 b)

degs_male <- read.csv("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_male <- degs_male[degs_male$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"),]
degs_male$Sex = "male"
degs_female <- read.csv("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_female$Sex = "female"
degs_female <- degs_female[degs_female$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"),]

degs_male <- degs_male[degs_male$p_adj.loc < 0.05 & abs(degs_male$logFC)>0.5,]
degs_female <- degs_female[degs_female$p_adj.loc < 0.05 & abs(degs_female$logFC)>0.5,]

degs_dereg <- full_join(degs_male[,colnames(degs_male) %in% c("gene", "cluster_id", "contrast")], degs_female[,colnames(degs_female) %in% c("gene", "cluster_id", "contrast")])


Correlation$Comparison <- factor(Correlation$Comparison, levels = c( "PDvsHC","MCIvsHC","ADvsHC"))

Correlation$CellAbbr <- labels_cell[Correlation$Celltype]
Correlation$CellAbbr <- factor(Correlation$CellAbbr, levels = labels_cell[celltype_cluster_order])
Correlation$Celltype <- factor(Correlation$Celltype, levels = celltype_cluster_order)

c <- seq(nrow(degs_dereg[, colnames(degs_dereg) %in% c("cluster_id", "contrast")]))
counts <- aggregate(c~., data=degs_dereg[,colnames(degs_dereg) %in% c("cluster_id", "contrast")], FUN=length)

counts$Celltype <- counts$cluster_id
counts$Comparison <- counts$contrast

data <- merge(Correlation, counts, by = c("Celltype", "Comparison"))

data$p_val[data$p_val == 0] <- min(data$p_val[data$p_val>0])

plot <- ggplot(data, aes(x=c, y = Value, size = -log10(p_val), color = Comparison)) + 
  geom_point()+ 
  theme_adrc()+ 
  scale_size(range = c(0,3))+
  facet_grid(.~Comparison) + 
  scale_color_manual(values = color_v,guide="none")+
  geom_text_repel(data=subset(data, p_val<0.05 & abs(c*Value) >100 & labels_cell[Celltype] %in% c("Naive CD4+ T", "Prolif. CD8+ T", "CD56- NK")), aes(label = labels_cell[Celltype]), color="black", size=2.5)+ 
  labs(size = "-log10\np-value", Comparison = "")+theme(panel.grid = element_line(colour = "gray92"), panel.border = element_rect(color = "black"))+
  geom_hline(yintercept = 0, linetype = 'dashed')  + xlab("Number of DEGs") + ylab("Correlation") + theme(legend.title=element_text(size=8), legend.text=element_text(size=6), aspect.ratio = 1.3)  #+ xlim(0,1100)+ ylim(-0.75,0.75)

figure3_b <- as.ggplot(fill_title(plot, color_v))
#save_plot("results/figures/ScatterPlotNumDEGSandCorrelation.pdf", figure3_b, base_height = 70, base_width=180, units="mm")

## Figure 3 c)

labels= ggplot(Correlation, aes(x=Celltype, y=1, fill=Celltype)) + geom_tile() +
  scale_fill_manual(values = color_v) +
  theme_void() + guides(fill="none")

plot <- ggplot(Correlation, aes(x = Celltype, y = Comparison, fill = Value))+ 
  geom_tile() +  
  geom_point(aes(size = ifelse(Correlation$p_val<0.05, ifelse(Correlation$Value > 0, "pos", "neg"), "no_dot"), shape = ifelse(Correlation$p_val<0.05, ifelse(Correlation$Value > 0, "pos", "neg"), "no_dot"))) +
  scale_size_manual(values = c(pos = 0.5, neg = 0.5, no_dot = NA), guide = "none")+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-1,1))+ 
  scale_shape_manual(values=c(pos = 19, neg = 19, no_dot = NA), guide = "none")+
  theme_adrc() + theme(axis.text.x = element_text(angle = 90, vjust = 0.7, hjust=1))+
  scale_x_discrete(name = "", labels = labels_cell[Correlation$Celltype]) + scale_y_discrete(name = "") + labs(fill='Correlation') +
  theme(legend.title=element_text(size=8), legend.text=element_text(size=6)) + 
  theme(legend.key.size = unit(3, 'mm'))
 
figure3_c <- plot %>% insert_top(labels, height=.3)
#ggsave("results/figures/Correlation_DEGs_Male_Female.pdf", figure3_c, width = 6, height = 2.2)

## Figure 3 d) 

a <- read.csv("results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
a <- a[a$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"),]

degs_male <- read.csv("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_male <- degs_male[degs_male$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"),]
degs_male$Sex = "male"
degs_female <- read.csv("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_female$Sex = "female"
degs_female <- degs_female[degs_female$contrast %in% c("ADvsHC", "PDvsHC", "MCIvsHC"),]

merged <- merge(degs_male, degs_female, by = c("gene", "cluster_id", "contrast"))

dge <- merge(a, merged, by = c("gene", "cluster_id", "contrast"), all = TRUE)

dge$color <- "none"
dge$color[dge$p_adj.loc.x < 0.05] <- "male"
dge$color[dge$p_adj.loc.y < 0.05 ] <- "female"
dge$color[dge$p_adj.loc.y < 0.05  & dge$p_adj.loc.x < 0.05 ] <- "male and female"

m <- dge[dge$cluster_id %in% c("Naive CD4+ T cells","Proliferating CD8+ T cell"),] #, "CD4+ effector memory T cell", "Follicular helper CD4 T cell", "Proliferating CD8+ T cell"),] Terminal effector CD4+ T cell
m$cl_abbr <- m$cluster_id
m$cl_abbr[m$cluster_id == "Naive CD4+ T cells"] <- "Naive CD4+ T"
m$cl_abbr[m$cluster_id == "Proliferating CD8+ T cell"] <- "Prolif. CD8+ T"
m$cl_abbr[m$cluster_id == "CD4+ effector memory T cell"] <- "Tem CD4+ T"
m$cl_abbr[m$cluster_id == "Follicular helper CD4 T cell"] <- "Tfh CD4+ T"
m$cluster_id <- factor(m$cluster_id, levels = c("Naive CD4+ T cells","Proliferating CD8+ T cell", "CD4+ effector memory T cell", "Follicular helper CD4 T cell"))


p <- ggplot(m, aes(x = logFC.x, y = logFC.y)) + 
  geom_point(aes(color = color), alpha = 0.7)+ 
  facet_grid(cols=vars(contrast), rows = vars(cluster_id)) +  xlim(-2.5,2.5)+ ylim(-2.5,2.5)+ 
  scale_fill_gradient2(low = "#d6604d", mid = "white", high = "#4393c3") +
  scale_colour_manual(name = "significant in", values = c("male and female" = "black", "female"= "#9C489A", "none"="lightgray","male"="#5D7322", "gray" = "gray")) + stat_cor(method = "pearson", label.x = -2, label.y = 2, size = 2.5)+
  geom_smooth(method=lm, color="black")+ theme_adrc() + theme(aspect.ratio = 1, legend.position = "bottom")+
  geom_hline(yintercept = 0, color="lightgrey", linetype="dashed") +
  geom_vline(xintercept = 0, color="lightgrey", linetype="dashed")+
  geom_hline(yintercept = 0.5, color="lightgrey", linetype="dotted")+
  geom_hline(yintercept = -0.5, color="lightgrey", linetype="dotted")+
  geom_vline(xintercept = 0.5, color="lightgrey", linetype="dotted")+
geom_vline(xintercept = -0.5, color="lightgrey", linetype="dotted")+ xlab("log2 Fold-change in male")+ ylab("log2 Fold-change in female")


figure3_d <-  fill_title(p, color_v)
#save_plot("results/figures/ScatterPlotFCComp.pdf", figure3_d, base_height = 5, base_width=8)

## Figure 3 e)

degs <- read.csv("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")

degs <- degs[degs$contrast %in% c("PDvsHC", "ADvsHC", "MCIvsHC"),]
degs <- degs[degs$p_adj.loc<0.05 & abs(degs$logFC)>0.5,]

degs_upset <- list(PDvsHC = degs$gene[degs$contrast == "PDvsHC"], ADvsHC = degs$gene[degs$contrast == "ADvsHC"],  MCIvsHC = degs$gene[degs$contrast == "MCIvsHC"])

m1 = make_comb_mat(degs_upset)


top_annotation = upset_top_annotation(m1,
                                      axis_param = list(gp = gpar(fontsize = 6)), 
                                      gp = gpar(fontsize = 8, fill = "black"), annotation_name_gp = grid::gpar(fontsize = 8), height = unit(1.5, "cm")
)
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))


plot <- UpSet(m1, comb_order = order(comb_size(m1)),   top_annotation= top_annotation, right_annotation = right_annotation, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8))


figure3_e <- as.ggplot(plot)
#pdf("results/figures/UpsetDisease_Male.pdf", height = 3, width = 5)
#figure3_e
#dev.off()

# Figure 3 f)

degs <- read.csv("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")

degs <- degs[degs$contrast %in% c("PDvsHC", "ADvsHC", "MCIvsHC"),]
degs <- degs[degs$p_adj.loc<0.05 & abs(degs$logFC)>0.5,]

degs_upset <- list(PDvsHC = degs$gene[degs$contrast == "PDvsHC"], ADvsHC = degs$gene[degs$contrast == "ADvsHC"],  MCIvsHC = degs$gene[degs$contrast == "MCIvsHC"])

m1 = make_comb_mat(degs_upset)

top_annotation = upset_top_annotation(m1,
                                      axis_param = list(gp = gpar(fontsize = 6)), 
                                      gp = gpar(fontsize = 8, fill = "black"), annotation_name_gp = grid::gpar(fontsize = 8), height = unit(1.5, "cm")
                                      )
right_annotation = upset_right_annotation(m1, gp = gpar(fill = color_v[rownames(m1)], fontsize = 2), annotation_name_gp = grid::gpar(fontsize = 8))

#bottom_annotation = upset


plot <- UpSet(m1, comb_order = order(comb_size(m1)),   top_annotation= top_annotation, right_annotation = right_annotation, row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 8))
figure3_f <- as.ggplot(plot)

#pdf("results/figures/UpsetDisease_Female.pdf", height = 3, width = 5)
#figure3_f
#dev.off()

#figure_e <- plot_grid(figure_e1, figure_e2, nrow = 2, labels = c("e", "f"))

# Figure 3 g)

degs_male <- read.csv("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_male$gender <- "male"
degs_female <- read.csv("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv", sep="\t")
degs_female$gender <- "female"
degs <- rbind(degs_female, degs_male)

unique(degs$contrast)
degs_ad <- degs[degs$contrast %in% c("ADvsHC"),]#, "PDMCIvsHC", "PDvsMCI", "PDvsHC", "ADvsPD"),]
genes <- as.data.frame(degs_ad$gene)

Alzheimer_kegg <- read.csv("AlzheimerGenes.csv", header=FALSE)

degs_kegg_ad <- degs_ad[gsub("^MT-", "", degs_ad$gene) %in% Alzheimer_kegg$V1 & degs_ad$contrast == "ADvsHC",]



freq <- as.data.frame(table(degs_kegg_ad$gene[degs_kegg_ad$p_adj.loc<0.05]))
use <- freq$Var1[freq$Freq >1]
degs_kegg_ad <- degs_kegg_ad[(degs_kegg_ad$gene %in% use) ,]


Cellt_to_use <- unique(degs_kegg_ad$cluster_id[degs_kegg_ad$p_adj.loc<0.05])
degs_kegg_ad <- degs_kegg_ad[ (degs_kegg_ad$cluster_id %in% Cellt_to_use),]

degs_kegg_ad$logFC[degs_kegg_ad$logFC>1.5] <- 1.5
degs_kegg_ad$logFC[degs_kegg_ad$logFC< (-1.5)] <- -1.5



degs_kegg_ad$cluster_id <- factor(degs_kegg_ad$cluster_id, levels = celltype_cluster_order)

female = sprintf(intToUtf8(9792))
male = intToUtf8(9794)

label_names <- list(
  'male'=male,
  'female'=female
)

sex_labeller <- function(variable,value){
  variable$gender <- unlist(lapply(variable$gender, function(x) ifelse(x=="female", female, male)))
  print((labels_cell[variable$cluster_id]))
  variable$cluster_id <- unlist(lapply(variable$cluster_id, function(x) as.character(labels_cell[x])))
  print(variable)                   
  return(variable)
}

p <- ggplot() + 
  geom_point(data=degs_kegg_ad, aes(x=cluster_id , y=gene,  size=-log10(p_adj.loc), fill=logFC, color = p_adj.loc<0.05 ), alpha = 0.8, shape = 21) +
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-1.5,1.5))+ 
  scale_color_manual(values=c("white","black"), guide = "none")+
  facet_nested(contrast~cluster_id+ gender, scales = "free_x",  labeller=sex_labeller,  remove_labels = "x")  + 
  scale_x_discrete(name = "", labels = element_blank()) +
  xlab("") + ylab("") + 
  theme_adrc() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + labs(size = "log10\np-value")

figure3_g <- as.ggplot(fill_title(p, color_v))
#save_plot("results/figures/Alzheimer_Kegg_tmp.pdf", figure3_g, base_height = 5, base_width=13) 



save_plot("results/figures/FiguresForPaper/Figure_3a.svg", plot_grid(figure3_a, labels = c("a")), base_width=210*(1/2.3), base_height=210*0.25, unit="mm")
save_plot("results/figures/FiguresForPaper/Figure_3b.svg", plot_grid(figure3_b, labels = c("b")), base_width=210*(1.3/2.3)*0.8, base_height=210*0.5*(1.1/2.1)*0.8, unit="mm")
save_plot("results/figures/FiguresForPaper/Figure_3c.svg", plot_grid(as.ggplot(figure3_c), labels = c("c")), base_width=210*(1.3/2.3), base_height=210*0.5*(0.8/2.1), unit="mm")
save_plot("results/figures/FiguresForPaper/Figure_3d.svg", plot_grid(figure3_d, labels = c("d")), base_width=210*(2/3), base_height=210*0.5, unit="mm")
save_plot("results/figures/FiguresForPaper/Figure_3e.svg", plot_grid(figure3_e, labels = c("e")), base_width=210*(1/3), base_height=210*0.5, unit="mm")
save_plot("results/figures/FiguresForPaper/Figure_3e.svg", plot_grid(figure3_e, labels = c("e")), base_width=210*(1/3), base_height=210*0.5, unit="mm")
save_plot("results/figures/FiguresForPaper/Figure_3f.svg", plot_grid(figure3_f, labels = c("e")), base_width=210*(1/3), base_height=210*0.5, unit="mm")
save_plot("results/figures/FiguresForPaper/Figure_3g.svg", plot_grid(figure3_g, labels = c("g")), base_width=210, base_height=210*0.5, unit="mm")


