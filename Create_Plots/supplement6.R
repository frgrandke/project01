library(facet_nested)
library(ggplot2)
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


source("ADRC_theme.R")
source("scripts/helper.R")


colors = fread("data/colors.csv")
color_v = colors$Color
names(color_v) = colors$ID

labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

pathways_f <- read.csv("Pipeline_Female/results/genetrail3/Pathway_Overview_deregulated.csv")
pathways_f$Sex <- "female"
pathways_m <- read.csv("Pipeline_Male/results/genetrail3/Pathway_Overview_deregulated.csv")
pathways_m$Sex <- "male"

pathways <- rbind(pathways_f, pathways_m)

table <- as.data.frame(table(pathways$celltype))
pathways <- pathways[pathways$celltype %in% table$Var1[order(table$Freq, decreasing = T)][1:10],]

#pathways$X.Name[order(pathways$P.value,  decreasing = F)][1:20]

pathways_to_use <- c()

for (d in c("ADvsHC", "MCIvsHC", "PDvsHC", "PDMCIvsHC")){
  for (s in c("female", "male")){
    freq.pathways <- as.data.frame(table(pathways$X.Name[pathways$disease == d & pathways$Sex == s]))
    pathways_to_use <- c(pathways_to_use, as.character(freq.pathways$Var1[order(freq.pathways$Freq, decreasing = T)][1:5]))
  }
}

# 
# freq.pathways <- as.data.frame(table(pathways$X.Name[pathways$disease == "PDvsHC"]))
# pathways_to_use <- c(pathways_to_use, as.character(freq.pathways$Var1[order(freq.pathways$Freq, decreasing = T)][1:5]))
# 
# freq.pathways <- as.data.frame(table(pathways$X.Name[pathways$disease == "MCIvsHC"]))
# pathways_to_use <- c(pathways_to_use, as.character(freq.pathways$Var1[order(freq.pathways$Freq, decreasing = T)][1:5]))
# 
# freq.pathways <- as.data.frame(table(pathways$X.Name[pathways$disease == "PDMCIvsHC"]))
# pathways_to_use <- c(pathways_to_use, as.character(freq.pathways$Var1[order(freq.pathways$Freq, decreasing = T)][1:5]))

path_sel <- pathways[pathways$P.value < 0.05 & pathways$X.Name %in% pathways_to_use,]
path_sel$X.Name <- factor(path_sel$X.Name, levels = unique(pathways_to_use))
path_sel$celltype <- factor(path_sel$celltype, levels = celltype_cluster_order)

female = sprintf(intToUtf8(9792))
male = intToUtf8(9794)

label_names <- list(
  'male'=male,
  'female'=female
)

sex_labeller <- function(variable,value){
  variable$Sex <- unlist(lapply(variable$Sex, function(x) ifelse(x=="female", female, male)))
  print((labels_cell[variable$cluster_id]))
  variable$celltype <- unlist(lapply(variable$celltype, function(x) as.character(labels_cell[x])))
  print(variable)                   
  return(variable)
}


plot <- ggplot(path_sel, aes(x = Sex, y = X.Name, fill = Regulation_direction)) + geom_tile() + 
  facet_nested(disease ~ celltype+Sex, scales = "free", space = "free",  labeller=sex_labeller,  remove_labels = "x")+ 
  scale_fill_gradient(low = "#4575b4", high = "#d73027")+ xlab("")+ ylab("")+
  theme_adrc()+ 
  theme(legend.position = "none")+scale_x_discrete(name = "", labels = element_blank())

save_plot("results/figures/Supplement5.pdf", fill_title(plot, color_v), base_width = 280, base_height = 200, unit = "mm", family = "sans", device = cairo_pdf)

