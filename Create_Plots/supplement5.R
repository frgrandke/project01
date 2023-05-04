library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(tidyverse)
library(cowplot)
library(data.table)
library(ggsignif)

colors = fread("data/colors.csv")
color_v = colors$Color
names(color_v) = colors$ID

labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

source("ADRC_theme.R")
source("scripts/helper.R")

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

metadata <- read.csv("results/metadata_with_celltype_cluster_proportions_and_new_biomarker_data.csv", sep=",")
metadata <- metadata[metadata$to_exclude_reason == "",]
metadata <- metadata[metadata$Diagnosis %in% c("Healthy Control", "Parkinson's Disease only", "Alzheimer's disease", "Parkinson's Disease with MCI", "Mild Cognitive Impairment"),]

sigFunc = function(x){
  if(x < 0.001){"***"} 
  else if(x < 0.01){"**"}
  else if(x < 0.05){"*"}
  else{NA}
}

colors = fread("data/colors.csv")
color_v = colors$Color
names(color_v) = colors$ID

name2short = c("All"="All", "Neurodegeneration"="", "Cognitive Impairment"="",
               "Healthy Control"="HC", "Parkinson's Disease"="PD", "Parkinson's Disease only"="PD",
               "Alzheimer's disease"="AD", "Parkinson's Disease with MCI"="PD-MCI",
               "Mild Cognitive Impairment"="MCI")

metadata$biogroup_short = name2short[as.character(metadata$Diagnosis)]

cols <- c("Estimated.Total.Intracranial.Volume","Left.hemisphere.cortical.gray.matter.volume","Right.hemisphere.cortical.gray.matter.volume","Total.cortical.gray.matter.volume","Left.hemisphere.cerebral.white.matter.volume","Right.hemisphere.cerebral.white.matter.volume","Total.cerebral.white.matter.volume","Subcortical.gray.matter.volume","Total.gray.matter.volume","Left.Lateral.Ventricle","Left.Inf.Lat.Vent","Left.Cerebellum.White.Matter","Left.Cerebellum.Cortex","Left.Thalamus.Proper","Left.Caudate","Left.Putamen","Left.Pallidum","X3rd.Ventricle","X4th.Ventricle","Brain.Stem","Left.Hippocampus","Left.Amygdala","CSF","Left.Accumbens.area","Left.VentralDC","Left.vessel","Left.choroid.plexus","Right.Lateral.Ventricle","Right.Inf.Lat.Vent","Right.Cerebellum.White.Matter","Right.Cerebellum.Cortex","Right.Thalamus.Proper","Right.Caudate","Right.Putamen","Right.Pallidum","Right.Hippocampus","Right.Amygdala","Right.Accumbens.area","Right.VentralDC","Right.vessel","Right.choroid.plexus","WM.hypointensities","non.WM.hypointensities","Optic.Chiasm","CC_Posterior","CC_Mid_Posterior","CC_Central","CC_Mid_Anterior","CC_Anterior","Left.Whole.Hippocampus","Right.Whole.Hippocampus","Left.Hippocampal.Tail","Right.Hippocampal.Tail","Left.Subiculum","Right.Subiculum","Left.Presubiculum","Right.Presubiculum","Left.Parasubiculum","Right.Parasubiculum","Left.CA1","Right.CA1","Left.CA3","Right.CA3","Left.CA4","Right.CA4","Left.Hippocampal.Fissure","Right.Hippocampal.Fissure","Left.Molecular.Layer.HP","Right.Molecular.Layer.HP","Left.GC.ML.DG","Right.GC.ML.DG","Left.Fimbria","Right.Fimbria","Left.HATA","Right.HATA","Age_at_CSF_measurement","Diagnosis_at_CSF_measurement","Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1","Quanterix4plex.CSF___Tau","Quanterix4plex.CSF___GFAP","Quanterix3plex.CSF___Ab40","Quanterix3plex.CSF___Ab42","Quanterix3plex.CSF___Tau","Quanterix.CSF___pTau181")
data <- metadata[,colnames(metadata) %in% c("ADRC study ID","SCMD", "PIDN","Sample_Batch","Processing_Batch","Sex",	"Age",	"Visit", "biogroup_short")]
data$info <- metadata[["Age.at.Scan"]]
data$info_col <- "Age.at.Scan"

for (c in cols){
  tmp <- metadata[,colnames(metadata) %in% c("ADRC study ID","SCMD", "PIDN","Sample_Batch","Processing_Batch","Sex",	"Age",	"Visit", "biogroup_short")]
  tmp$info <- as.numeric(metadata[[c]])
  tmp$info_col <- c
  data <- rbind(data, tmp)
}

data <- data[which(!is.na(data$info)),]

data$biogroup_short <- factor(data$biogroup_short, levels = c("HC", "MCI", "AD", "PD-MCI", "PD"))

celltype_proportion_sample_1 = ggplot(data, aes(x=biogroup_short, y=info, fill=biogroup_short, color=biogroup_short)) +
  geom_boxplot(outlier.size=0.15, lwd=0.2) + facet_wrap(~info_col, scales = "free_y") +
  scale_fill_manual(values=color_v) +
  scale_color_manual(values=(color_v)) +
  #scale_y_continuous(labels=scales::percent_format()) + xlab("") +
  theme_cowplot(10) #+ theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                            #strip.text=element_text(size=7)) + ylab("PBMCs")

all_pairwise_diag_comparisons = list(c("HC", "MCI"), c("HC", "AD"), c("HC", "PD"), c("HC", "PD-MCI")) #combn(diagnosis_order, 2, simplify = F)


celltype_proportion_sample_1_w_sig =  celltype_proportion_sample_1 + 
  geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="greater"), color="black", tip_length=0) +
  geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="less"), color="black", tip_length=0) 
celltype_proportion_sample_1_w_sig

cols <- c("Estimated.Total.Intracranial.Volume","Left.hemisphere.cortical.gray.matter.volume","Right.hemisphere.cortical.gray.matter.volume","Total.cortical.gray.matter.volume","Left.hemisphere.cerebral.white.matter.volume","Right.hemisphere.cerebral.white.matter.volume","Total.cerebral.white.matter.volume","Subcortical.gray.matter.volume","Total.gray.matter.volume")
cols <- c("Left.Lateral.Ventricle","Right.Lateral.Ventricle","Left.Inf.Lat.Vent","Right.Inf.Lat.Vent","Left.Cerebellum.White.Matter","Right.Cerebellum.White.Matter","Left.Cerebellum.Cortex","Right.Cerebellum.Cortex","Left.Thalamus.Proper","Right.Thalamus.Proper","Left.Caudate","Right.Caudate","Left.Putamen","Right.Putamen","Left.Pallidum","Right.Pallidum","Left.Hippocampus","Right.Hippocampus","Left.Amygdala","Right.Amygdala","Left.Accumbens.area","Right.Accumbens.area","Left.VentralDC","Right.VentralDC","Left.vessel","Right.vessel","Left.choroid.plexus","Right.choroid.plexus")
cols <- c("X3rd.Ventricle","X4th.Ventricle","Brain.Stem","CSF")
cols <- c("CC_Posterior","CC_Mid_Posterior","CC_Central","CC_Mid_Anterior","CC_Anterior")
cols <- c("Left.Whole.Hippocampus","Right.Whole.Hippocampus","Left.Hippocampal.Tail","Right.Hippocampal.Tail","Left.Subiculum","Right.Subiculum","Left.Presubiculum","Right.Presubiculum","Left.Parasubiculum","Right.Parasubiculum", "Left.CA1","Right.CA1","Left.CA3","Right.CA3","Left.CA4","Right.CA4","Left.Hippocampal.Fissure","Right.Hippocampal.Fissure","Left.Molecular.Layer.HP","Right.Molecular.Layer.HP","Left.GC.ML.DG","Right.GC.ML.DG","Left.Fimbria","Right.Fimbria", "Left.HATA","Right.HATA")
cols <- c("Quanterix.CSF___pTau181", "Quanterix3plex.CSF___Ab42",
          "Quanterix4plex.CSF___GFAP",
          "Quanterix4plex.CSF___Tau","Quanterix3plex.CSF___Ab40","Quanterix3plex.CSF___Tau",
          "Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1",
          "Quanterix3plex.CSF___Ab42/Ab40",
          "CSF P-Tau181", "CSF Total Tau", "CSF AB42", "CSF AB40", "CSF AB42/AB40")

cols <-  c("Quanterix.CSF___pTau181", "Quanterix3plex.CSF___Ab42",
           "Quanterix4plex.CSF___GFAP",
           "Quanterix4plex.CSF___Tau","Quanterix3plex.CSF___Ab40","Quanterix3plex.CSF___Tau",
           "Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1",
           "Quanterix3plex.CSF___Ab42/Ab40")

data <- metadata[,colnames(metadata) %in% c("ADRC study ID","SCMD", "PIDN","Sample_Batch","Processing_Batch","Sex",	"Age",	"Visit", "biogroup_short")]
data$info <- metadata[[cols[1]]]
data$info_col <- cols[1]

for (c in cols[-1]){
  tmp <- metadata[,colnames(metadata) %in% c("ADRC study ID","SCMD", "PIDN","Sample_Batch","Processing_Batch","Sex",	"Age",	"Visit", "biogroup_short")]
  tmp$info <- as.numeric(metadata[[c]])
  tmp$info_col <- c
  data <- rbind(data, tmp)
}

data <- data[which(!is.na(data$info)),]
data$info_col <- factor(data$info_col, levels = cols)

data$biogroup_short <- factor(data$biogroup_short, levels = c("HC", "MCI", "AD", "PD-MCI", "PD"))

#data$info_col <- unlist(lapply(as.character(data$info_col), function(x) substr(x, 16, nchar(x))))

labeller <- function(value){
  value <- str_replace(value, "___", "")
  return(value)
}


celltype_proportion_sample_1 = ggplot(data, aes(x=biogroup_short, y=info, fill=biogroup_short, color=biogroup_short)) +
  geom_boxplot(outlier.size=0.15, lwd=0.2) + facet_wrap(~info_col, scales = "free_y",nrow=1, labeller = labeller) +
  scale_fill_manual(values=color_v) +
  scale_color_manual(values=darken(color_v)) +
  #scale_y_continuous(labels=scales::percent_format()) + xlab("") +
  theme_adrc() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 1.2)#,
#strip.text=element_text(size=7)) + ylab("PBMCs")

all_pairwise_diag_comparisons = list(c("HC", "MCI"), c("HC", "AD"), c("HC", "PD"), c("HC", "PD-MCI")) #combn(diagnosis_order, 2, simplify = F)


celltype_proportion_sample_1_w_sig =  celltype_proportion_sample_1 + 
  geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="greater"), color="black", tip_length=0) +
  geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="less"), color="black", tip_length=0) 
celltype_proportion_sample_1_w_sig

save_plot("results/figures/Alzheimer_CSF_marker_significance.svg", celltype_proportion_sample_1_w_sig, base_height = 40, base_width = 220, unit = "mm")

#cols <- c("Left.Whole.Hippocampus","Right.Whole.Hippocampus","Left.Hippocampal.Tail","Right.Hippocampal.Tail","Left.Subiculum","Right.Subiculum","Left.Presubiculum","Right.Presubiculum","Left.Parasubiculum","Right.Parasubiculum", "Left.CA1","Right.CA1","Left.CA3","Right.CA3","Left.CA4","Right.CA4","Left.Hippocampal.Fissure","Right.Hippocampal.Fissure","Left.Molecular.Layer.HP","Right.Molecular.Layer.HP","Left.GC.ML.DG","Right.GC.ML.DG","Left.Fimbria","Right.Fimbria", "Left.HATA","Right.HATA")
cols <- c("Left.hemisphere.cortical.gray.matter.volume","Right.hemisphere.cortical.gray.matter.volume","Total.cortical.gray.matter.volume","Subcortical.gray.matter.volume","Total.gray.matter.volume","Estimated.Total.Intracranial.Volume","Left.hemisphere.cerebral.white.matter.volume","Right.hemisphere.cerebral.white.matter.volume","Total.cerebral.white.matter.volume")

labeller <- function(variable,value){
  print(value)
  value <- str_replace_all(value, "[:punct:]", " ")
  value <- str_replace(value, "gray", "\ngray")
  value <- str_replace(value, "white", "\nwhite")
  value <- str_replace(value, "Intracranial", "\nIntracranial")
  return(value)
}
cols <- c("Left.hemisphere.cortical.gray.matter.volume","Right.hemisphere.cortical.gray.matter.volume","Total.cortical.gray.matter.volume","Subcortical.gray.matter.volume","Total.gray.matter.volume","Estimated.Total.Intracranial.Volume","Left.hemisphere.cerebral.white.matter.volume","Right.hemisphere.cerebral.white.matter.volume","Total.cerebral.white.matter.volume")

#cols <- c("Total.cerebral.white.matter.volume","Total.cortical.gray.matter.volume","Left.hemisphere.cortical.gray.matter.volume","Right.hemisphere.cortical.gray.matter.volume","Left.hemisphere.cerebral.white.matter.volume","Right.hemisphere.cerebral.white.matter.volume","Subcortical.gray.matter.volume","Total.gray.matter.volume","Estimated.Total.Intracranial.Volume")

data <- metadata[,colnames(metadata) %in% c("ADRC study ID","SCMD", "PIDN","Sample_Batch","Processing_Batch","Sex",	"Age",	"Visit", "biogroup_short")]
data$info <- metadata[[cols[1]]]
data$info_col <- cols[1]

for (c in cols[-1]){
  tmp <- metadata[,colnames(metadata) %in% c("ADRC study ID","SCMD", "PIDN","Sample_Batch","Processing_Batch","Sex",	"Age",	"Visit", "biogroup_short")]
  tmp$info <- as.numeric(metadata[[c]])
  tmp$info_col <- c
  data <- rbind(data, tmp)
}

data <- data[which(!is.na(data$info)),]
data$info_col <- factor(data$info_col, levels = cols)

data$biogroup_short <- factor(data$biogroup_short, levels = c("HC", "MCI", "AD", "PD-MCI", "PD"))



celltype_proportion_sample_1 = ggplot(data, aes(x=biogroup_short, y=info, fill=biogroup_short, color=biogroup_short)) +
  geom_boxplot(outlier.size=0.15, lwd=0.2) + facet_wrap(~info_col, scales = "free_y",nrow = 2, labeller = labeller) +
  scale_fill_manual(values=color_v) +
  scale_color_manual(values=darken(color_v)) +
  xlab("") + ylab("") +
  #scale_y_continuous(labels=scales::percent_format()) + xlab("") +
  theme_adrc() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), aspect.ratio = 1.2)#,
#strip.text=element_text(size=7)) + ylab("PBMCs")

all_pairwise_diag_comparisons = list(c("HC", "MCI"), c("HC", "AD"), c("HC", "PD"), c("HC", "PD-MCI")) #combn(diagnosis_order, 2, simplify = F)


celltype_proportion_sample_1_w_sig =  celltype_proportion_sample_1 + 
  geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="greater"), color="black", tip_length=0) +
  geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="less"), color="black", tip_length=0) 
celltype_proportion_sample_1_w_sig

save_plot("results/figures/Alzheimer_CSF_marker_significance_2.svg", celltype_proportion_sample_1_w_sig, base_height = 80, base_width = 220, unit = "mm")





correlations <- 'results/correlation/celltype_and_celltype_cluster_correlations.csv'
correlations <- fread(correlations, sep=',')
correlations <- correlations[correlations$Diagnosis != "all"]

correlations_cluster = correlations[grepl("celltype_cluster", V1)]
correlations = correlations[grepl("celltype_cluster", V1)]

colors = fread("data/colors.csv", strip.white = F)
color_v = colors$Color
names(color_v) = colors$ID

correlations_p = correlations[V2 %in% c("Quanterix.CSF___pTau181", "Quanterix3plex.CSF___Ab42",
                                        "Quanterix4plex.CSF___GFAP",
                                        "Quanterix4plex.CSF___Tau","Quanterix3plex.CSF___Ab40","Quanterix3plex.CSF___Tau",
                                        "Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1",
                                        "Quanterix3plex.CSF___Ab42/Ab40")]

correlations_p = correlations[V2 %in% c(
  "Quanterix4plex.CSF___GFAP",
  "Quanterix4plex.CSF___NFL","Quanterix4plex.CSF___UCHL1",
  "CSF P-Tau181", "CSF Total Tau", "CSF AB42", "CSF AB40", "CSF AB42/AB40")]

correlations_p = correlations_p[Diagnosis %in% c("Parkinson's Disease only", "Healthy Control")]
correlations_p$V1 = unlist(lapply(correlations_p$V1, function(x) substr(x, 18, nchar(x))))
correlations_p[V2 == "Quanterix3plex.CSF___Tau", V2:="Tau (1)"]
correlations_p[V2 == "Quanterix4plex.CSF___Tau", V2:="Tau (2)"]
correlations_p[, V2:=gsub(".*___", "", V2)]
correlations_p[, V2:=gsub("CSF ", "", V2)]
order <- unlist(lapply(unique(correlations_p$V1), function(x) sum(correlations_p$Pearson[correlations_p$V1 == x & correlations_p$Diagnosis == "Parkinson's Disease only"])))
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

figure_b = fill_title(plotn, color_v)
save_plot("results/figures/sig_correlation_csf_marker_celltypes_heatmap_pearson_Parkinson.svg", figure_b, base_width=220, base_height=50, unit="mm")

correlations_p = correlations[correlations$Diagnosis %in% c("Parkinson's Disease only", "Healthy Control")]

names <- as.data.frame(table(correlations$V2[correlations$pearson_p < 0.05 & correlations$Diagnosis %in% c("Parkinson's Disease only")]))
sel_names <- names$Var1[order(names$Freq, decreasing = T)][1:15]

correlations_p = correlations[V2 %in% sel_names]

correlations_p = correlations_p[Diagnosis %in% c("Parkinson's Disease only", "Healthy Control")]
correlations_p$V1 = unlist(lapply(correlations_p$V1, function(x) substr(x, 18, nchar(x))))


List <- (lapply(unique(correlations_p$V1), function(x) {
    t <- correlations_p[correlations_p$V1 == x & correlations_p$Diagnosis == "Parkinson's Disease only"]; 
    t[[as.character(x)]] <- t$Pearson
    t <- t  %>% select("V2", as.character(x))
    as.data.frame(t)
}))
matrix <- Reduce(inner_join,List)
rownames(matrix) <- matrix$V2
matrix <- as.matrix(matrix[,-1])

order_V1 <- hclust(dist(matrix))$order
order_V2 <- hclust(dist(t(matrix)))$order


order <- unlist(lapply(unique(correlations_p$V1), function(x) sum(correlations_p$Pearson[correlations_p$V1 == x & correlations_p$Diagnosis == "Parkinson's Disease only" & correlations_p$pearson_p < 0.05])))
correlations_p$V1 <- factor(correlations_p$V1, level = colnames(matrix)[order_V2])

order <- unlist(lapply(unique(correlations_p$V2), function(x) sum(correlations_p$Pearson[correlations_p$V2 == x & correlations_p$Diagnosis == "Parkinson's Disease only" & correlations_p$pearson_p < 0.05])))
correlations_p$V2 <- factor(correlations_p$V2, level = rownames(matrix)[order_V1])

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

figure_b = fill_title(plotn, color_v)
save_plot("results/figures/sig_correlation_csf_marker_celltypes_heatmap_pearson_Parkinson_Brain.svg", figure_b, base_width=220, base_height=70, unit="mm")

################################################################################

#correlations_p = correlations[correlations$Diagnosis %in% c("Alzheimer's disease", "Healthy Control")]

names <- as.data.frame(table(correlations$V2[correlations$pearson_p < 0.05 & correlations$Diagnosis %in% c("Alzheimer's disease")]))
sel_names <- names$Var1[order(names$Freq, decreasing = T)][1:15]

correlations_p = correlations[V2 %in% sel_names]

correlations_p = correlations_p[Diagnosis %in% c("Alzheimer's disease", "Healthy Control")]
correlations_p$V1 = unlist(lapply(correlations_p$V1, function(x) substr(x, 18, nchar(x))))


List <- (lapply(unique(correlations_p$V1), function(x) {
  t <- correlations_p[correlations_p$V1 == x & correlations_p$Diagnosis == "Alzheimer's disease"]; 
  t[[as.character(x)]] <- t$Pearson
  t <- t  %>% select("V2", as.character(x))
  as.data.frame(t)
}))
matrix <- Reduce(full_join,List)
rownames(matrix) <- matrix$V2
matrix <- as.matrix(matrix[,-1])

order_V1 <- hclust(dist(matrix))$order
order_V2 <- hclust(dist(t(matrix)))$order


#order <- unlist(lapply(unique(correlations_p$V1), function(x) sum(correlations_p$Pearson[correlations_p$V1 == x & correlations_p$Diagnosis == "Parkinson's Disease only" & correlations_p$pearson_p < 0.05])))
correlations_p$V1 <- factor(correlations_p$V1, level = colnames(matrix)[order_V2])

#order <- unlist(lapply(unique(correlations_p$V2), function(x) sum(correlations_p$Pearson[correlations_p$V2 == x & correlations_p$Diagnosis == "Parkinson's Disease only" & correlations_p$pearson_p < 0.05])))
correlations_p$V2 <- factor(correlations_p$V2, level = rownames(matrix)[order_V1])

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

figure_b = fill_title(plotn, color_v)
save_plot("results/figures/sig_correlation_csf_marker_celltypes_heatmap_pearson_Alzheimer's disease_Brain.svg", figure_b, base_width=220, base_height=70, unit="mm")

