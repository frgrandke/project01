library(ggplot2)
library(cowplot)
source("scripts/helper.R")
source("ADRC_theme.R")

read_brain_input <- function(path, celltypes, sex){
  degs_brain <- NA
  for (c in cell_types){
    if (sex == "m") {path <- paste("../../../../../../Downloads/dataset_latent_var/human_cortex_AD_ct_", c, "_M_dataset.csv", sep = "")} else {path <- paste("../../../../../../Downloads/dataset_latent_var/human_cortex_AD_ct_", c, "_F_dataset.csv", sep = "")}
    print(path)
    data <- read.csv(path, sep = ",")
    print(data)
    data$cluster_id <- c
    if (is.na(degs_brain)[1]) {degs_brain <- data} else {degs_brain <- rbind(degs_brain, data)}
  }
  degs_brain$tissue <- "brain"
  degs_brain$logFC <- degs_brain$avg_log2FC
  degs_brain$p_adj.loc <- degs_brain$p_val_adj
  degs_brain$gene <- degs_brain$X
  
  degs_brain <- degs_brain[, colnames(degs_brain) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]
  
  return(degs_brain)
}

read_pbmc_input <- function(path){
  
  degs_pbmc <- read.csv(path, sep = "\t")
  
  degs_pbmc <- degs_pbmc[degs_pbmc$contrast == "ADvsHC" ,] 
  degs_pbmc$tissue <- "pbmc"
  degs_pbmc <- degs_pbmc[, colnames(degs_pbmc) %in% c("tissue", "logFC", "p_adj.loc", "gene", "cluster_id")]
  return(degs_pbmc)
}

# plotComparison <- function(brain_data, pbmc_data, celltypes, celltype_to_order){
#   
#   degs <- rbind(brain_data, pbmc_data)
#   print(degs[degs$gene %in% degs$gene[degs$tissue == "brain"] & degs$tissue == "pbmc",])
#   degs <- degs[degs$gene %in% degs$gene[degs$tissue == "brain"],]
#   degs <- degs[degs$gene %in% degs$gene[degs$tissue == "pbmc"],]
#   degs <- degs[degs$cluster_id %in% celltypes | degs$tissue == "brain",]
#   degs <- degs[degs$gene %in% degs$gene[degs$cluster_id %in% celltypes],]
#   
#   #degs$gene <- factor(degs$gene, levels = unique(degs$gene[degs$tissue == "pbmc"][order(degs$logFC[degs$tissue == "pbmc"], decreasing = T)]))
#   gene_list_to_order <- unique(degs$gene)
#   fc_to_order <- lapply(gene_list_to_order, function(gene) {
#     if (gene %in% degs$genes[degs$cluster_id %in% celltype_to_order]) {degs$logFC[degs$cluster_id %in% celltype_to_order & degs$gene == gene]} else 
#     {degs$logFC[degs$gene == gene][1]}
#   })
#   levels <- unique(degs$gene[degs$cluster_id %in% celltype_to_order][order(degs$logFC[degs$cluster_id %in% celltype_to_order], decreasing = T)])
#   degs$gene <- factor(degs$gene, levels = levels)
#   
#   degs$tissue <- factor(degs$tissue, levels = c("pbmc", "brain"))
#   degs$cluster_id <- factor(degs$cluster_id, levels = c("Excitatory neuron", "Inhibitory neuron", "Pvalb neuron", "OPC","Astrocyte",  "Oligodendrocyte", "Microglia", "Vascular", celltype_cluster_order))
#   
#   
#   degs$logFC[degs$logFC>0.5] <- 0.5
#   degs$logFC[degs$logFC<(-0.5)] <- (-0.5)
#   degs$p_adj.loc[degs$p_adj.loc == 0] <- min(degs$p_adj.loc[degs$p_adj.loc != 0])
# 
#   degs$p_adj.loc_norm <- NA
#   degs$p_adj.loc_norm[degs$tissue == "pbmc"] <- -log(degs$p_adj.loc[degs$tissue == "pbmc"]) / max(-log(degs$p_adj.loc[degs$tissue == "pbmc"]))
#   degs$p_adj.loc_norm[degs$tissue == "brain"] <- -log(degs$p_adj.loc[degs$tissue == "brain"]) / max(-log(degs$p_adj.loc[degs$tissue == "brain"]))
# 
#   return(degs)
# }




plotComparison <- function(brain_data, pbmc_data, celltypes, celltype_to_order){
  degs <- rbind(brain_data, pbmc_data)
  degs <- degs[degs$gene %in% degs$gene[degs$tissue == "brain"],]
  degs <- degs[degs$cluster_id %in% celltypes | degs$tissue == "brain",]
  degs <- degs[degs$gene %in% degs$gene[degs$cluster_id %in% celltypes],]
  degs$gene <- factor(degs$gene, levels = unique(degs$gene[degs$cluster_id %in% celltype_to_order][order(degs$logFC[degs$cluster_id %in% celltype_to_order], decreasing = T)]))
  degs$tissue <- factor(degs$tissue, levels = c("pbmc", "brain"))
  degs$cluster_id <- factor(degs$cluster_id, levels = c("Excitatory neuron", "Inhibitory neuron", "Pvalb neuron", "OPC","Astrocyte",  "Oligodendrocyte", "Microglia", "Vascular", celltype_cluster_order))
  degs$logFC[degs$logFC>1] <- 1
  degs$logFC[degs$logFC<(-1)] <- (-1)
  degs$p_adj.loc[degs$p_adj.loc < 10^(-10)] <- 10^(-10)
  plot <- ggplot() +
    geom_point(data = degs, aes(y= cluster_id, x=gene, fill = logFC, size = -log(p_adj.loc), colour = p_adj.loc<0.05 ), shape = 21) +
    facet_grid(tissue~., scales = "free", space = "free")+
    scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-1,1))+
    scale_colour_manual(values=c("white","black"), guide = "none")+
    theme_adrc()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(plot)
}

female_genes <- read.csv("Alzheimer_female", header = F)$V1
male_genes <- read.csv("Alzheimer_male", header = F)$V1

cell_types <- c("Astrocyte", "Excitatory neuron", "Inhibitory neuron", "Microglia", "Oligodendrocyte", "OPC", "Pvalb neuron", "Vascular")

cellt_female <- c("Naive CD4+ T cells","CD56- Natural Killer cell","Naive CD8+ T cell","Proliferating CD8+ T cell", "CD4+ T-Helper 1 Cell")
cellt_male <- c("Naive CD4+ T cells","Naive CD8+ T cell")
celltype_to_order <- c("Naive CD4+ T cells")

degs_brain_male <- read_brain_input(path, cell_types, sex = "m")
degs_pbmc_male <- read_pbmc_input("Pipeline_Male/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv")
degs_brain_female <- read_brain_input(path, cell_types, sex = "f")
degs_pbmc_female <- read_pbmc_input("Pipeline_Female/results/de/volcano/ds_pb_limma_voom_all_col_celltype_cluster_with_additional_stats.csv")

#plot_male_degs <- plotComparison(degs_brain_male[degs_brain_male$gene %in% male_genes,], degs_pbmc_male, cellt_male, cellt_male)

plot_male <- ggplot() + 
  geom_tile(data = plot_male_degs, aes(y= cluster_id, x=gene, fill = logFC)) + 
  facet_grid(tissue~., scales = "free", space = "free")+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-0.5,0.5))+ 
  #scale_colour_manual(values=c("white","black"), guide = "none")+
  theme_adrc()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#plot_female_degs <- plotComparison(degs_brain_female[degs_brain_female$gene %in% female_genes,], degs_pbmc_female, cellt_female, cellt_female)

plot_female <- ggplot() + 
  geom_tile(data = plot_female_degs, aes(y= cluster_id, x=gene, fill = logFC)) + 
  facet_grid(tissue~., scales = "free", space = "free")+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-0.5,0.5))+ 
  #scale_colour_manual(values=c("white","black"), guide = "none")+
  theme_adrc()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

save_plot("results/Plot_Brain_male_filtered_New.svg", plot_male, base_width = 220, base_height = 50, unit = "mm")
save_plot("results/Plot_Brain_female_filtered_New.svg", plot_female, base_width = 230, base_height = 50, unit = "mm")


# genes <- unique(read.csv(paste("../../../../../../Downloads/results/Alzheimer_female_pathways/", cell_types[1], "_F.csv", sep = ""))$names)
# genes_selected <-c("CHCHD10", "CORC1A", "BTF3", "FAU", "GAPDH", "HOPX", "IFITM2", "ITGB7", "LSM2", "MDH2", "PPDPF", "RAC2", "RPL27", "RPL3", "RPS5", "UBC", "YDJC", "ECH1", "BOLA1", "ITGB7", 
#   "PPN1", "HSPA11", "BDP1", "IBA57") # , "TTC3", "PLGRKT", "HTATIP2", "CEP290", "MSH3"
# high <- unique(degs_brain_female$gene)[order(unlist(lapply(unique(degs_brain_female$gene), function(x) sum(degs_brain_female$p_adj.loc[degs_brain_female$gene == x]<0.05))), decreasing = T)][1:10]
# genes <- as.character(unique(degs_brain_female$gene[degs_brain_female$gene %in% genes[!genes %in% genes_selected]]))
# degs_brain_female <- degs_brain_female[degs_brain_female$gene %in% c(genes, genes_selected, high),]
# degs_brain_female$gene <- factor(degs_brain_female$gene, levels = c("CHCHD10", "CORC1A", "BTF3", "FAU", "GAPDH", "HOPX", "IFITM2", "ITGB7", "LSM2", "MDH2", "PPDPF", "RAC2", "RPL27", "RPL3", "RPS5", "UBC", "YDJC", "ECH1", "BOLA1", 
#                                                                     "PPN1", "HSPA11", "BDP1", "IBA57",genes , high[!high %in% c(genes_selected, genes)]))
plot_female_degs <- plotComparison(degs_brain_female[degs_brain_female$p_adj.loc < 0.05 & degs_brain_female$gene %in% female_genes & abs(degs_brain_female$logFC) > 0.5 ,], degs_pbmc_female[degs_pbmc_female$p_adj.loc < 0.05,], cellt_female, celltype_to_order)
unique(degs_brain_female$gene[degs_brain_female$gene %in% genes])

plot_female <- ggplot() +
  geom_point(data = plot_female_degs, aes(y= cluster_id, x=gene, fill = logFC, size = p_adj.loc_norm, colour = p_adj.loc<0.05 ), shape = 21) +
  facet_grid(tissue~., scales = "free", space = "free")+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-0.5,0.5))+
  scale_colour_manual(values=c("white","black"), guide = "none")+
  theme_adrc()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  scale_size(range = c(0, 5))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())


plot_male_degs <- plotComparison(degs_brain_male[degs_brain_male$p_adj.loc < 0.05  & degs_brain_male$gene %in% male_genes & abs(degs_brain_male$logFC) > 0.5 ,], degs_pbmc_male[degs_pbmc_male$p_adj.loc <0.05,], cellt_male, celltype_to_order)

plot_male <- ggplot() +
  geom_point(data = plot_male_degs, aes(y= cluster_id, x=gene, fill = logFC, size = p_adj.loc_norm, colour = p_adj.loc<0.05 ), shape = 21) +
  facet_grid(tissue~., scales = "free", space = "free")+
  scale_fill_gradient2(low="#4575b4", mid="white", high="#d73027", limits = c(-0.5,0.5))+
  scale_colour_manual(values=c("white","black"), guide = "none")+
  theme_adrc()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_size(range = c(0, 4.5))+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())


save_plot("results/Plot_Brain_male_new.svg", plot_male, base_width = 220, base_height = 100, unit = "mm")
save_plot("results/Plot_Brain_female_new.svg", plot_female, base_width = 220, base_height = 100, unit = "mm")

