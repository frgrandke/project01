library(Seurat)
library(viridisLite)
library(ggplot2)
library(cowplot)
#library(Biobase)
adrc <- "results/annotated_celltypes/CompleteObjectAnnotated.rds"
adrc_obj <- readRDS(adrc)

labels = fread("data/CelltypeMapping.csv")
labels_cell = labels$Abbr
names(labels_cell) =  labels$ID

adrc_obj$cellt_cl_tmp <- labels_cell[adrc_obj$celltype_cluster]

source("scripts/helper.R")
adrc_obj$celltype_cluster <- factor(adrc_obj$celltype_cluster, levels = celltype_cluster_order)
adrc_obj$CellCluster <- factor(adrc_obj$CellCluster, levels = CellCluster_order)
adrc_obj$cellt_cl_tmp <- factor(adrc_obj$cellt_cl_tmp, levels = labels_cell[celltype_cluster_order])
cellt_cl_tmp

plot1 <- DotPlot(adrc_obj, group.by = "CellCluster", feature = c("CD14", "FCGR3A", # Monocytes
                                                                 "IGHA1", "IGLC1", "MZB1", "JCHAIN",  # Plasma cells
                                                                 "PPBP", "PF4", "GP1BB", "GNG11", # Megacaryocytes, "AC000093.1", "OAZ1"
                                                                 "CD1C", "FCER1A", "LILRA4", "GPR183", # Dendritic cells, "GZMB", "ITM2C", "CADM1"
                                                                 "NCAM1", "NCR1", # NK cells #"FCGR3A"
                                                                 "CD3D", "CD3E", # T cells
                                                                 "CD19", "MS4A1"# B cells
                                                                 ), scale=T, assay = "RNA", dot.scale = 5) +
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ scale_size(limits = c(0,95), range = c(1, 6)) + scale_colour_gradient(low = "lightgray", high = "blue", limits = c(-2,2.5))

sub <- subset(adrc_obj,subset = CellCluster == "CD4+ T cell" & celltype_cluster != "Proliferating CD4+ T cell ")
plot2 <- DotPlot(sub, group.by = "cellt_cl_tmp", feature = c("MKI67", # Prolif
                                                                       "CXCR5", "CD40LG", "CD4", # T FH "IFNG",
                                                                        "ZEB2", "SAMD3", "TYROBP", # Terma #"CD28", "FAS", "KLRC2", 
 "IL1B", 
 "CCR7","SELL", #Naive "TCF7", "FOXP1", "IL7R", "PNISR", 
 "TRF", #  "PTPRC", 
 "IKZF2","FOXP3", "CTLA4",
  
 "GATA3", "IL4I1", "IL4R", "CCR4", "CCR6",    # TH 2
"TBX21", "PRF1", "GZMK", # TH1"LTA", "TNF", "CXCR3", 
"IL21R", "IRF4", "POU2AF1", # TH17"RORC", 
"CD28","CD27", "ITGAL"), scale=T,assay = "RNA", dot.scale = 5)+
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ scale_size(limits = c(0,90), range = c(1, 6)) + scale_colour_gradient(low = "lightgray", high = "blue", limits = c(-2,2.5))


sub <- subset(adrc_obj,subset = CellCluster == "CD8+ T cell")
plot3 <- DotPlot(sub,  group.by = "cellt_cl_tmp", feature = c("IL7R", "TCF7", # Naive
"MKI67", # Proliferating
"CCR7", "SELL", "CD28", # Memory"PTPRC", , "ITGAL"
"CXCR3", "FAS"), scale=T, assay = "RNA") +
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ scale_size(limits = c(0,90), range = c(1, 6)) + scale_colour_gradient(low = "lightgray", high = "blue", limits = c(-2,2.5))


#
sub <- subset(adrc_obj,subset = CellCluster == "B cell")
sub$celltype_cluster <- factor(sub$celltype_cluster, levels = c("Naive B cell", "Transitional B cell", "IgM Memory B cell", "Classical Memory B cell", "Double negative 1 B cell", 
                                                                "Double negative 2 B cell", "Double negative 3 B cell", "Double negative 4 B cell", "B cell"))


plot4 <- DotPlot(sub, group.by = "cellt_cl_tmp", feature = c( "IGHD",  "IGHM", # General "MS4A1", 
                                                                                                                     "IGLL5", "TCL1A",# Trans"VPREB3" , "PCDH9", 
                                                                                                                     "PLPP5", "FCER2", # Naive"IL4R",
                                                                                                                     "CD79A",  #neg in C-Mem  "CD27","CD19", "B2M", "CD74", 
                                                                       "MARCKS", # M-memory# ND4
                                                                       "JCHAIN", "IGHA2", #DN1
                                                                       "EMP3", "CIB1", "PSAP", # DN2, "CD72", "DAPP1",
                                                                       "FCRL5", # DN3 IGHM
                                                                       "LCN10"), scale=T, assay = "RNA") +
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ scale_size(limits = c(0,90), range = c(1, 6)) + scale_colour_gradient(low = "lightgray", high = "blue", limits = c(-2,2.5))

plot5 <- DotPlot(subset(adrc_obj,subset = CellCluster %in% c("NK cell", "NKT-like cell")), group.by = "cellt_cl_tmp", feature = c("FCGR3A", # CD56-Dim
"CD3D", # NKT-like cells
"NCAM1" # CD56-Bright
), scale=T, assay = "RNA") +
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ scale_size(limits = c(0,90), range = c(1, 6)) + scale_colour_gradient(low = "lightgray", high = "blue", limits = c(-2,2.5))

#
sub <- subset(adrc_obj,subset = CellCluster %in% c("Monocyte"))
plot6 <- DotPlot(sub, group.by = "cellt_cl_tmp", feature = c("FCGR3A", # CD16 Monocytes
"CD14",  # CD14 Monocytes
"S100A8" # Prolif. Monocytes
), scale=T, assay = "RNA") +
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ scale_size(limits = c(0,90), range = c(1, 6)) + scale_colour_gradient(low = "lightgray", high = "blue", limits = c(-2,2.5))

plot7 <- DotPlot(subset(adrc_obj,subset = CellCluster %in% c("Dendritic cell")), group.by = "cellt_cl_tmp", feature = c("CD1C", "FCER1A", # cDC2
"ITM2C", "LILRA4",  # pDC
"CADM1" # cDC1
), scale=T, assay = "RNA") +
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))+ scale_size(limits = c(0,90), range = c(1, 6)) + scale_colour_gradient(low = "lightgray", high = "blue", limits = c(-2,2.5))



# plots <- align_plots(p3, p1, align = 'v', axis = 'l')
# # then build the bottom row
# bottom_row <- plot_grid(plots[[2]], p2, labels = c('B', 'C'), label_size = 12)
# 
# # then combine with the top row for final plot
# plot_grid(plots[[1]], bottom_row, labels = c('A', ''), label_size = 12, ncol = 1)
leg <- get_legend(plot2)

plot1 <- plot1 + NoLegend()
plot2 <- plot2 + NoLegend()
plot3 <- plot3 + NoLegend()
plot4 <- plot4 + NoLegend()
plot5 <- plot5+ NoLegend()
plot6 <- plot6 + NoLegend()
plot7 <- plot7 + NoLegend()

plot_1 <- plot_grid(plot1, labels = c("A"))
plot_2 <- plot_grid(plot2, labels = c("B"))
#plot <- plot_grid(plot2, plot4, labels = c("A", "B"), nrow =2) 

plot_3 <- plot_grid(plot4,  plot3, labels = c( 'C', "D"), label_size = 12, nrow = 1, rel_widths = c(1.3, 1))
#plot5 <- plot5 + scale_colo
plot_4 <- plot_grid(plot5, plot6 , plot7, leg, labels = c( 'E', 'F', "G"), label_size = 12, nrow = 1, rel_widths = c(1,1.1, 1.1, 0.5))

fig_marker <- plot_grid(plot_1, plot_2, plot_3, plot_4, nrow = 4, rel_heights = c(1,1,1,0.5))

#plot2+ scale_x_discrete(name = "", labels = labels_cell[correlations_p$V1])

#fig_marker = plot_grid(plot1, plot2, plot3, plot4, plot5, plot6, plot7, ncol = 1)
save_plot("results/figures/Marker_genes.pdf", fig_marker, base_height = 300, base_width = 200, unit="mm")
head(adrc_obj@meta.data)
unique(adrc_obj$CellCluster)
