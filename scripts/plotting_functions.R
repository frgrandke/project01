plot_density_difference = function(data1, data2) {
  xrng = range(data1$V1, data2$V1)
  yrng = range(data1$V2, data2$V2)
  d1 = MASS::kde2d(data1$V1, data1$V2, lims=c(xrng, yrng), n=200)
  d2 = MASS::kde2d(data2$V1, data2$V2, lims=c(xrng, yrng), n=200)
  
  identical(d1$x, d2$x) # TRUE
  identical(d1$y, d2$y) # TRUE
  
  # Calculate the difference between the 2d density estimates
  diff12 = d1
  diff12$z = d1$z - d2$z
  
  ## Melt data into long format
  # First, add row and column names (x and y grid values) to the z-value matrix
  rownames(diff12$z) = diff12$x
  colnames(diff12$z) = diff12$y
  
  # Now melt it to long format
  diff12.m = melt(diff12$z, id.var=rownames(diff12))
  names(diff12.m) = c("V1","V2","z")
  
  #limits = c(0, max(abs(diff12.m$z)))
  #limits[1] = -limits[2]
  limits = c(0, 0.015)
  limits[1] = -limits[2]
  
  # Plot difference between densities
  ggplot(diff12.m, aes(V1, V2, z=z, fill=z)) +
    geom_tile() +
    stat_contour(aes(colour=..level..), binwidth=0.0005) +
    scale_fill_gradient2(name=expression(Delta*Density), low="#4575b4",mid="white", high="#d73027", midpoint=0, limits=limits, breaks=scales::pretty_breaks(5)) +
    scale_colour_gradient2(low=scales::muted("#4575b4"), mid="white", high=scales::muted("#d73027"), midpoint=0) +
    coord_cartesian(xlim=xrng, ylim=yrng) +
    guides(colour=FALSE) + ylab("") + xlab("") + theme_adrc()+ theme(aspect.ratio = 1, legend.position="none")
}

getCelltypeProportion_Boxplot <- function(input_data, celltypes1, celltypes2, all_pairwise_diag_comparisons){
  p_1 = ggplot(input_data[input_data$Celltype %in% celltypes1,], aes(x=Diagnosis, y=value, fill=Diagnosis, color=Diagnosis)) +
    geom_boxplot(outlier.size=0.15, lwd=0.2) + facet_wrap(~Celltype, nrow=1) +
    scale_fill_manual(values=color_v) +
    scale_color_manual(values=darken(color_v)) +
    scale_y_continuous(labels=scales::percent_format()) + xlab("") +
    theme_adrc() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("PBMCs")
  p_1_w_sig =  p_1 +
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="greater"), color="black", tip_length=0) +
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="less"), color="black", tip_length=0)
  
  
  cell_p_other = input_data[input_data$Celltype %in% celltypes2,]
  cell_p_other_wo_large_outlier = cell_p_other[cell_p_other$value < 0.6,]
  
  p_2 = ggplot(cell_p_other_wo_large_outlier, aes(x=Diagnosis, y=value, fill=Diagnosis, color=Diagnosis)) +
    geom_boxplot(outlier.size=0.15, data=cell_p_other_wo_large_outlier, lwd=0.2) + facet_wrap(~Celltype, nrow=1) + #
    scale_fill_manual(values=color_v) +
    scale_color_manual(values=darken(color_v)) +
    scale_y_continuous(labels=scales::percent_format(), limits = c(0,0.6)) + xlab("") + #
    theme_adrc() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())
  p_2_w_sig =  p_2 +
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="greater"), color="black", tip_length=0) +
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="less"), color="black", tip_length=0)
  
  cell_p_other = input_data[input_data$Celltype %in% celltypes3,]
  cell_p_other_wo_large_outlier = cell_p_other[cell_p_other$value < 0.2,]
  
  p_3 = ggplot(cell_p_other_wo_large_outlier, aes(x=Diagnosis, y=value, fill=Diagnosis, color=Diagnosis)) +
    geom_boxplot(outlier.size=0.15, data=cell_p_other_wo_large_outlier, lwd=0.1) + facet_wrap(~Celltype, nrow=1) + #
    scale_fill_manual(values=color_v) +
    scale_color_manual(values=darken(color_v)) +
    scale_y_continuous(labels=scales::percent_format(), limits = c(0,0.2)) + xlab("") + #, limits = c(0,0.15)
    theme_adrc() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())
  p_3_w_sig =  p_3 +
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="greater"), color="black", tip_length=0) +
    geom_signif(comparisons=all_pairwise_diag_comparisons, map_signif_level = sigFunc, margin_top = -0.15, size=0.2, textsize=2.5, step_increase = 0.045, test.args=list(alternative="less"), color="black", tip_length=0)
  
  
  celltype_proportions_p = plot_grid(fill_title(p_1_w_sig, color_v),
                                     fill_title(p_2_w_sig, color_v),
                                     fill_title(p_3_w_sig, color_v),
                                     rel_widths = c(2.9, 10.1, 7.1), nrow=1) #
  
  return(celltype_proportions_p)
}