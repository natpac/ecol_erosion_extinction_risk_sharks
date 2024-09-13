show_final_pcoa_fun <- function(pcoa_object_with_arrows, sp_plot=F, axis1 = 1, axis2 = 2,
                           variable_scalar = 10, max.overlaps_label = 10) {
  
  given_pcoa = pcoa_object_with_arrows
  
  # Store axes names
  kept_axes = paste0("Axis.", c(axis1, axis2))
  
  # data.frame with values from species
  pcoa_df = as.data.frame(given_pcoa$vectors)
  
  # data.frame with variable covariances
  arrow_df = as.data.frame(given_pcoa$U/variable_scalar)
  arrow_df$variable = rownames(arrow_df) # Add variable names
  
  # Eigen values per axes
  eigen_values = given_pcoa$values$Relative_eig[c(axis1, axis2)]
  names(eigen_values) = kept_axes
  
  # Axes labels in plot
  axes_labs = lapply(kept_axes, function(x) {
    paste0(gsub(".", " ", x, fixed = TRUE),
           " (", round(eigen_values[x] * 100, 2), "%)")
  })
  
  # Compute plot limits to get symetrical limits on each axis (= square plot)
  # limits = apply(given_pcoa$vectors[, kept_axes], 2, range)
  limits = apply(arrow_df[, kept_axes], 2, range)
  
  x_lim = c(limits[1, 1] - diff(limits[,1])/10,
            limits[2, 1] + diff(limits[,1])/10)
  
  y_lim = c(limits[1, 2] - diff(limits[,2])/10,
            limits[2, 2] + diff(limits[,2])/10)
  
  
  
  # Final plot
  ggplot(pcoa_df, aes_string(x = kept_axes[1],
                             y = kept_axes[2])) +
    # Add dashed axes
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    # Add species points
    {if(sp_plot==T)geom_point()} +
    # Add arrows
    geom_segment(data = arrow_df, x = 0, y = 0, alpha = 0.7, size = 1.2,
                 arrow = arrow(length = unit(3, "mm")),
                 aes_string(xend = kept_axes[1], yend = kept_axes[2])
                 , col = 'blue4') +
    # Add arrow labels
    ggrepel::geom_label_repel(data = arrow_df, aes(label = variable),
                              size = 3, fontface = "bold",
                              lineheight = 0.7
                              , max.overlaps = max.overlaps_label
    ) +
    # Specified limits to get square plot
    # to see the effect of this you can replace this line by
    # coord_equal() +
    coord_fixed(xlim = x_lim, ylim = y_lim) +
    # Labels & Theme
    labs(x = axes_labs[[1]],
         y = axes_labs[[2]]) +
    theme_bw() +
    theme(panel.grid = element_blank())
}