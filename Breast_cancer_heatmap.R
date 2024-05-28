geneSetDat <- readRDS("geneSetDat.rds")

## plots heatmaps of correlation matrices for each gene set
library(ggplot2)
library(gridExtra)

pdf("heatmap_correlation_matrices.pdf", width = 12, height = 15)  # Adjust width and height as needed

plots <- list()
plot_count <- 1
for (ss in 1:length(P0)) {
  dat_tmp <- geneSetDat[[ss]]$dat
  
  for (group in 1:3) {
    corr_mat <- cor(dat_tmp[[group]])
    cor_df <- as.data.frame(as.table(corr_mat))
    names(cor_df) <- c("Var1", "Var2", "Correlation")
    
    heatmap_plot <- ggplot(data = cor_df, aes(Var1, Var2, fill = Correlation)) +
      geom_tile() +
      scale_fill_gradient2(low = "#084594", mid = "white", high = "#e15759", 
                           midpoint = 0, name = "Correlation", limits = c(-1, 1)) +
      labs(title = paste0(geneSetDat[[ss]]$GOID, ", Group ", group), 
           x = "Genes", y = "Genes") +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            aspect.ratio = 1)
    
    plots[[plot_count]] <- heatmap_plot
    plot_count <- plot_count + 1
  }
  
  # If we have 12 plots, arrange them in a 3x5 grid and start a new page
  if (plot_count > 15) {
    grid.arrange(grobs = plots, ncol = 3, nrow = 5)
    plots <- list()  # Reset the list
    plot_count <- 1  # Reset the plot counter
  }
}

# If there are any remaining plots, arrange them as well
if (length(plots) > 0) {
  grid.arrange(grobs = plots, ncol = 3, nrow = 5)
}
dev.off()