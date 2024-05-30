################################################################################
######## plot heatmaps of correlation matrices for each gene set
################################################################################

## install necessary packages
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("gridExtra", quietly = TRUE)) install.packages("gridExtra")

## load necessary packages
library(ggplot2)
library(gridExtra)

## load gene sets data obtained from "geneSetDat.rds" file
geneSetDat <- readRDS("geneSetDat.rds")

## plot heatmaps of correlation matrices for each gene set into a pdf file
pdf("heatmap_correlation_matrices.pdf", width = 12, height = 15)

plots <- list()
plot_count <- 1
for (ss in 1:length(geneSetDat)) {
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
  
  # If we have 15 plots, arrange them in a 3x5 grid and start a new page
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


#### check Condition 2 for each group
corr_summ <- matrix(NA, 3, length(geneSetDat))
for (ss in 1:length(geneSetDat)) {
  dat_tmp <- geneSetDat[[ss]]$dat
  p_tmp <- geneSetDat[[ss]]$p
  
  corr_mat <- cor(dat_tmp[[1]])
  corr_summ[1, ss] <- max(sapply(c(1:nrow(corr_mat)), 
                                 function(z) sum(abs(corr_mat[z, ][-z]))/p_tmp))
  
  corr_mat <- cor(dat_tmp[[2]])
  corr_summ[2, ss] <- max(sapply(c(1:nrow(corr_mat)), 
                                 function(z) sum(abs(corr_mat[z, ][-z]))/p_tmp))
  
  corr_mat <- cor(dat_tmp[[3]])
  corr_summ[3, ss] <- max(sapply(c(1:nrow(corr_mat)), 
                                 function(z) sum(abs(corr_mat[z, ][-z]))/p_tmp))
}

## draw boxplots for each group
library(tidyr)
library(ggplot2)
df <- as.data.frame(t(corr_summ))
colnames(df) <- c("Group1", "Group2", "Group3")
long_df <- pivot_longer(df, cols = everything(), names_to = "Group", values_to = "Value")

pdf(file = "real_data_condition_2_boxplot.pdf", width = 10, height = 8)
# Draw the boxplot
ggplot(long_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = 21, outlier.fill = "red", 
               outlier.color = "black", outlier.size = 3) +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette from RColorBrewer
  labs(title = "Boxplot", x = "Group", y = "Value") +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(size = 12),
    legend.position = "none"
  ) +
  geom_jitter(width = 0.2, size = 2, color = "black", alpha = 0.7) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'dashed', color = 'grey'))
dev.off()