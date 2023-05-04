library(Seurat)
library(viridisLite)
library(cowplot)
library(ggpubr)
library(data.table)
library(pbapply)
source("scripts/helper.R")

pbmc = readRDS("results/annotated_celltypes/CompleteObjectAnnotated.rds")
pbmc = pbmc[,pbmc$celltype != "Unknown/doublets"]
pbmc$biogroup_short = factor(name2short[as.character(pbmc$Diagnosis)], levels=diagnosis_order)
  
colors = fread("data/colors.csv")
color_v = colors$Color
names(color_v) = colors$ID

fig_mito_perc = ggplot(pbmc@meta.data, aes(x="", y = percent.mt, fill=biogroup_short)) +
  #geom_violin(color="#4575b4", fill="#74add1") +
  geom_violin() +
  geom_boxplot(width=0.2, fill="#74add1") +
  theme_cowplot(10) + ylab("Mitochondrial genes") + xlab("") +
  scale_y_continuous(labels=scales::percent_format(scale = 1), breaks = scales::pretty_breaks(5)) +
  scale_fill_manual(values = color_v) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

fig_mito_perc_per_diag = ggplot(pbmc@meta.data, aes(x=biogroup_short, y = percent.mt, fill=biogroup_short)) +
  geom_violin() +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  #stat_compare_means(label.y=0.7, method="anova") +
  theme_cowplot(10) + ylab("Mitochondrial genes") + xlab("") +
  scale_y_continuous(labels=scales::percent_format(scale = 1), breaks = scales::pretty_breaks(5)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())
  
fig_genes = ggplot(pbmc@meta.data, aes(x="", y = nFeature_RNA)) +
  geom_violin(color="#4575b4", fill="#74add1") +
  geom_boxplot(width=0.2, fill="#74add1") +
  theme_cowplot(10) + ylab("Genes per cell") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

fig_genes_per_diag = ggplot(pbmc@meta.data, aes(x=biogroup_short, y = nFeature_RNA, fill=biogroup_short)) +
  geom_violin() +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  #stat_compare_means(label.y=0.7) +
  theme_cowplot(10) + ylab("Genes per cell") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

fig_umis = ggplot(pbmc@meta.data, aes(x="", y = nCount_RNA)) +
  geom_violin(color="#4575b4", fill="#74add1") +
  geom_boxplot(width=0.2, fill="#74add1") +
  theme_cowplot(10) + ylab("UMIs per cell") + xlab("") +
  scale_y_log10(labels=scales::label_number_si(), breaks = c(100,500,1000,2000,5000,10000,25000,50000)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

fig_umis_per_diag = ggplot(pbmc@meta.data, aes(x=biogroup_short, y = nCount_RNA, fill=biogroup_short)) +
  geom_violin() +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  #stat_compare_means(label.y=0.7) +
  theme_cowplot(10) + ylab("UMIs per cell") + xlab("") +
  scale_y_log10(labels=scales::label_number_si(), breaks = c(100,500,1000,2000,5000,10000,25000,50000)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

cells_per_sample = data.table(Cells=table(pbmc$Sample))
cells_per_sample[, biogroup:=pbmc$biogroup_short[match(Cells.V1, pbmc$Sample)]]
fig_cells_per_sample = ggplot(cells_per_sample, aes(x="", y = Cells.N)) +
  geom_violin(color="#4575b4", fill="#74add1") +
  geom_boxplot(width=0.2, fill="#74add1", outlier.shape = NA) +
  theme_cowplot(10) + ylab("Cells per sample") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

print("Cells per sample")
print(mean(cells_per_sample$Cells.N))
print(sd(cells_per_sample$Cells.N))

fig_cells_per_sample_per_diag = ggplot(cells_per_sample, aes(x=biogroup, y = Cells.N, fill=biogroup)) +
  geom_violin() +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  #stat_compare_means(label.y=0.7) +
  theme_cowplot(10) + ylab("Cells per sample") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

umis_per_sample = as.data.table(pbmc@meta.data[, c("nCount_RNA", "Sample")])[, list(umis=sum(nCount_RNA)), by="Sample"]
umis_per_sample[, biogroup:=pbmc$biogroup_short[match(Sample, pbmc$Sample)]]
fig_umis_per_sample = ggplot(umis_per_sample, aes(x="", y = umis)) +
  geom_violin(color="#4575b4", fill="#74add1") +
  geom_boxplot(width=0.2, fill="#74add1") +
  theme_cowplot(10) + ylab("UMIs per sample") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

print("UMIs per sample")
print(mean(umis_per_sample$umis))
print(sd(umis_per_sample$umis))


fig_umis_per_sample_per_diag = ggplot(umis_per_sample, aes(x=biogroup, y = umis, fill=biogroup)) +
  geom_violin() +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  #stat_compare_means(label.y=40e6, label.x = 3) +
  theme_cowplot(10) + ylab("UMIs per sample") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

genes_per_sample = pbsapply(unique(pbmc$Sample), function(s){ sum(Matrix::rowSums(pbmc@assays$RNA@counts[,pbmc$Sample == s]) > 0) }, cl=64)
genes_per_sample_df = data.table(Genes=genes_per_sample, Sample=unique(pbmc$Sample))
genes_per_sample_df[, biogroup:=pbmc$biogroup_short[match(Sample, pbmc$Sample)]]

fig_genes_per_sample = ggplot(genes_per_sample_df, aes(x="", y = Genes)) +
  geom_violin(color="#4575b4", fill="#74add1") +
  geom_boxplot(width=0.2, fill="#74add1") +
  theme_cowplot(10) + ylab("Genes per sample") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())

print("Genes per sample")
print(mean(genes_per_sample_df$Genes))
print(sd(genes_per_sample_df$Genes))



fig_genes_per_sample_per_diag = ggplot(genes_per_sample_df, aes(x=biogroup, y = Genes, fill=biogroup)) +
  geom_violin() +
  geom_boxplot(width=0.2, outlier.shape = NA) +
  #stat_compare_means(label.y=40e6, label.x = 3) +
  theme_cowplot(10) + ylab("Genes per sample") + xlab("") +
  scale_y_continuous(labels=scales::label_number_si(), breaks = scales::pretty_breaks(5)) +
  scale_fill_manual(values = color_v) +
  scale_color_manual(values = darken(color_v)) +
  theme(legend.position = "none", axis.ticks.x = element_blank())


fig_abcde = plot_grid(fig_umis_per_diag, fig_genes_per_diag, fig_mito_perc_per_diag, fig_cells_per_sample_per_diag, fig_umis_per_sample_per_diag, fig_genes_per_sample_per_diag, nrow = 2)
save_plot("results/figures/qc_figure.pdf", fig_abcde, base_height = 100, base_width = 200, unit="mm")

stop()

# CD4 marker, IL7R, CD4, low NKG7, high CD3D, low cd8
# Cd8 marker, CD8A
# nk marker, NKG7, GNLY high, CD3d low

# IL7R, CD4, NKG7, CD3D, CD8A
# CST3, LYZ, FCGR3A, CD74, PPBP
# JCHAIN, MZB1, IGLL5, LGALS1, HSPA5, MS4A1, IL32, TNFAIP3, TXNDC5, FOS (low -> immature B)

markers = c("CD3D","IL7R", "CD4", "CD8A", "NKG7", "CST3", "LYZ", "FCGR3A", "CD74", "PPBP",
            "JCHAIN", "MZB1", "IGLL5", "LGALS1", "HSPA5", "MS4A1", "IL32", "TNFAIP3", "TXNDC5", "FOS")

DotPlot(pbmc, features=c("CD3D","IL7R", "CD4", "CD8A", "NKG7", "CST3", "LYZ", "FCGR3A", "CD74", "PPBP",
                         "MZB1", "MS4A1", "FOS"))

Seurat::DoHeatmap(pbmc, features=names(celltype_marker))
celltype_marker = c("CD3D"="TNK cells","IL7R"="TNK cells", "CD4"= "TNK cells", "CD8A"="TNK cells", "NKG7"="TNK cells", "CST3"="Myeloid", "LYZ"="Myeloid", "FCGR3A"="Myeloid", "CD74"="Myeloid", "PPBP"="Myeloid",
"MZB1"="B cells", "MS4A1"="B cells", "FOS"="B cells")

celltype_marker_expr_df = as.data.table(as.matrix(Matrix::t(pbmc@assays$RNA@data[names(celltype_marker),])), keep.rownames=TRUE)
celltype_marker_expr_df_melt = as.data.table(melt(celltype_marker_expr_df))
celltype_marker_expr_df_melt[, celltype:=pbmc@meta.data[rn,]$celltype]

ggplot(celltype_marker_expr_df_melt[celltype %in% c("CD4 T cell", "CD8 T cell", "NK cell") & variable %in% names(celltype_marker)[celltype_marker == "TNK cells"]], aes(x=variable, y=value)) + facet_grid(cols=vars(celltype)) + geom_violin() + coord_flip()


marker_plot = DotPlot(pbmc, features=c("CD3D","IL7R", "CD4", "CD8A", "NKG7", "CST3", "LYZ", "FCGR3A", "CD74", "PPBP",
                         "MZB1", "MS4A1", "FOS")) +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
    scale_color_viridis_c(option="magma", name="Scaled average expression") +
    guides(size=guide_legend(title="Percent expressed", override.aes=list(shape=21, colour="black", fill="white"))) + ylab("") + xlab("")

Idents(pbmc) = factor(pbmc$celltype, levels=celltype_order)
marker_plot2 = DotPlot(pbmc, features=c("CD3D","IL7R", "CD4", "CD8A", "NKG7", 
                                       "MS4A1", "FOS", "MZB1",
                                       "FCGR3A", "CST3", "LYZ", "CD74", "PPBP"), scale=TRUE, dot.scale = 8) +
  guides(size=guide_legend(title="Percent\nexpressed", override.aes=list(shape=21, colour="black", fill="white")),
         color=guide_colorbar(title="Scaled\naverage\nexpression")) + ylab("") + xlab("") + theme_cowplot(10) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))


#fig_abcde = plot_grid(fig_umis, fig_genes, fig_mito_perc, fig_umis_per_sample, fig_cells_per_sample, nrow = 1)

fig_abcde = plot_grid(fig_umis_per_diag, fig_genes_per_diag, fig_mito_perc_per_diag, fig_cells_per_sample_per_diag, fig_umis_per_sample_per_diag, fig_genes_per_sample_per_diag, nrow = 2)
fig = plot_grid(fig_abcde, marker_plot2, nrow=2)

save_plot("results/figures/qc_figure.pdf", fig, base_height = 200, base_width = 200, unit="mm")


#0 -> CD4 (IL7R, CD4)
#1 -> CD8 (CD8A)
#2 -> CD4 (No CD8/NKG7)
#3 -> NK (NKG7 & GNLY high, low CD3D)
#4 -> CD4 (low NKG7, Cd8, high CD3D, CD52)
#5 -> CD4 (low NKG7, Cd8, high CD3D)
#6 -> CD4 (low NKG7, high CD3D) (maybe CD8?)


