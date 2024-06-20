################################################################################
######## This file includes R code for real data analysis
################################################################################


################################################################################
######## Part 1: obtain the 149 gene sets
################################################################################


## Install packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("Biobase", quietly = TRUE)) BiocManager::install("Biobase")
if (!require("breastCancerMAINZ", quietly = TRUE)) 
  BiocManager::install("breastCancerMAINZ")
if (!require("limma", quietly = TRUE)) BiocManager::install("limma")
if (!require("hgu95av2.db", quietly = TRUE)) BiocManager::install("hgu95av2.db")

## load packages
library(breastCancerMAINZ)
library(Biobase)
library(limma)
library(hgu95av2.db)

################################################################################
######## search for gene sets corresponding to different GO terms
######## and contains at least 25 genes
################################################################################

## load the mainz expression data from package "breastCancerMAINZ"
data(mainz)
dat <- Biobase::exprs(mainz)
## get the group id of each observation 
id <- Biobase::pData(mainz)[, 9]

## normalize the data using 'normalizeQuantiles' in the 'limma' package
dat <- limma::normalizeQuantiles(dat, ties=TRUE)
dat0 <- t(dat)

#### functions to identify genes sets corresponding to different GO terms.
#### Use 'hgu95av2GO2ALLPROBES' in 'hgu95av2.db' package to maps between 
#### manufacturer IDs and Gene Ontology (GO) IDs.
GOTFfun <- function(GOID) {
  x <- hgu95av2.db::hgu95av2GO2ALLPROBES[[GOID]]
  unique(x)
}

GOTFDat <- function(GOID, dat0, id) {
  TFs <- unique(unlist(lapply(GOID, GOTFfun)))
  
  inSel <- match(TFs, colnames(dat0), nomatch=0)
  
  if (sum(inSel!=0) > 1) {
    dat_G1 <- dat0[which(id==1), inSel]
    dat_G2 <- dat0[which(id==2), inSel]
    dat_G3 <- dat0[which(id==3), inSel]
  } else if (sum(inSel!=0) == 1) {
    dat_G1 <- matrix(dat0[which(id==1), inSel], sum(id==1), 1)
    dat_G2 <- matrix(dat0[which(id==2), inSel], sum(id==2), 1)
    dat_G3 <- matrix(dat0[which(id==3), inSel], sum(id==3), 1)
  }
  return(list("dat_G1" = dat_G1, "dat_G2" = dat_G2, "dat_G3" = dat_G3))
}

#### obtain all GO IDs that link to manufacturer IDs in 'dat0'
## The manufacturer IDs are given by the column names of 'dat0'
geneIDs <- colnames(dat0)
tmp <- AnnotationDbi::select(hgu95av2.db, keys=geneIDs, column="GO")
GOIDs <- unique(tmp$GO)
GOIDs <- GOIDs[-which(is.na(GOIDs))]
# In total, we obtain 2576 candidate GO IDs that link to manufacturer IDs

#### obtain gene sets corresponding to each GO ID with at least 25 genes
geneSetDat <- vector("list", length(GOIDs))
i <- 1;  j <- 1
while (j <= length(GOIDs)) {
  GOID <- GOIDs[j]
  dat_tmp <- GOTFDat(GOID, dat0, id)
  # print(ncol(dat_tmp$dat_G1))
  if (ncol(dat_tmp$dat_G1) >= 25) {
    dat_G1 <- dat_tmp$dat_G1
    dat_G2 <- dat_tmp$dat_G2
    dat_G3 <- dat_tmp$dat_G3
    
    dat <- list(dat_G1, dat_G2, dat_G3)
    N <- c(nrow(dat_G1), nrow(dat_G2), nrow(dat_G3))
    p <- ncol(dat_G1)
    
    geneSetDat[[i]] = list("GOID"=GOID, "dat"=dat, "N"=N, "p"=p)
    i <- i + 1
  }
  j <- j + 1
}
geneSetDat1 <- vector("list", i-1)
for (ss in 1:(i-1)) {
  geneSetDat1[[ss]] <- geneSetDat[[ss]]
}
geneSetDat <- geneSetDat1
# In total, we have 149 gene sets with at least 25 genes.

saveRDS(geneSetDat, "geneSetDat.rds")

## get the number of genes in each gene set
P0 <- unlist(lapply(geneSetDat, function(x) ncol(x$dat[[1]])))
summary(P0)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 25.00   29.00   40.00   57.44   63.00  271.00


################################################################################
######## Part 2: check Condition 2 for each group
################################################################################

################################################################################
######## plot heatmaps of correlation matrices for each gene set
################################################################################

## install necessary packages
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("gridExtra", quietly = TRUE)) install.packages("gridExtra")

## load necessary packages
library(tidyr)
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


################################################################################
######## calculate metric to quantitatively examine Condition 2
################################################################################

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

## draw boxplots of the metric for each group
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


################################################################################
######## Part 3: MANOVA on gene sets data
################################################################################

## Install packages
if (!require("mvtnorm", quietly = TRUE)) install.packages("mvtnorm")
if (!require("Gmedian", quietly = TRUE)) install.packages("Gmedian")
if (!require("extraDistr", quietly = TRUE)) install.packages("extraDistr")
if (!require("doRNG", quietly = TRUE)) install.packages("doRNG")
if (!require("moments", quietly = TRUE)) install.packages("moments")

## load packages
library(mvtnorm)
library(Gmedian)
library(extraDistr)
library(doRNG)

## load functions for MANOVA
source("./HD_MANOVA_GM_functions.R")

## load gene sets data
geneSetDat <- readRDS("geneSetDat.rds")

# use parallel computing
no_cores <- parallel::detectCores() - 1
cl <- parallel::makeCluster(no_cores)
doParallel::registerDoParallel(cl)
R_seed <- 202405
set.seed(R_seed)

B <- 100000 ## number of bootstrap iteration
real_MANOVA <- foreach(ss = 1:length(geneSetDat),
                       .options.RNG = R_seed,
                       .packages=c("mvtnorm", "Gmedian", "extraDistr")) %dopar% {
                         dat <- geneSetDat[[ss]]$dat
                         N <- geneSetDat[[ss]]$N
                         p <- geneSetDat[[ss]]$p
                         K <- 3
                         
                         res_median <- MANOVA_Median(dat, K=3, N, p, n_boot=B)
                         res_Lin2021 <- MANOVA_Mean_Lin2021(dat, K=3, N, p, n_boot=B)
                         
                         list("res_median"=res_median,
                              "res_Lin2021"=res_Lin2021)
                       }
parallel::stopCluster(cl)
# saveRDS(real_MANOVA, "res_real_data_MANOVA.rds")

# real_MANOVA <- readRDS("./res_real_data_MANOVA.rds")

#### summarize the results to get the p-values of each method
res_summ <- data.frame(matrix(NA, length(real_MANOVA), 5))
names(res_summ) <- c("median_non_std", "median_std_T_s", "median_std_T_s_tilde",
                     "mean_Lin2021_non_std", "mean_Lin2021_std")

for (i in 1:length(real_MANOVA)) {
  res_median <- real_MANOVA[[i]]$res_median
  res_Lin2021 <- real_MANOVA[[i]]$res_Lin2021
  
  res_summ$median_non_std[i] <- mean(res_median$test_stat_median<res_median$res_median_boot)
  res_summ$median_std_T_s[i] <- mean(res_median$test_stat_median_std<res_median$res_median_std_boot_old)
  res_summ$median_std_T_s_tilde[i] <- mean(res_median$test_stat_median_std<res_median$res_median_std_boot)
  res_summ$mean_Lin2021_non_std[i] <- mean(res_Lin2021$test_stat_mean<res_Lin2021$res_mean_boot)
  res_summ$mean_Lin2021_std[i] <- mean(res_Lin2021$test_stat_mean_std<res_Lin2021$res_mean_std_boot)
}

## adjust the p-values using the Benjamini{Yekutieli procedure
res_summ_adjust <- matrix(0, nrow(res_summ), 5)
for (j in 1:5) {
  res_summ_adjust[, j] <- p.adjust(res_summ[,j], method="BY")
}
colnames(res_summ_adjust) <- c("median_non_std", "median_std_T_s",
                               "median_std_T_s_tilde",
                               "mean_Lin2021_non_std", "mean_Lin2021_std")
saveRDS(res_summ_adjust, "real_data_adjusted_p_values.rds")

#### the number of gene sets identified by each method
apply(res_summ_adjust, 2, function(x) sum(x<0.01))
# median_non_std median_std_T_s median_std_T_s_tilde mean_Lin2021_non_std mean_Lin2021_std 
#             69            145                  141                   64              140


################################################################################
####### plot the marginal kurtoses of all genes in the additional 5 gene sets 
####### identified by "median_non_std" but not by "mean_Lin2021_non_std"
################################################################################
library(breastCancerMAINZ)
library(Biobase)
library(moments)
res_summ_adjust <- readRDS("./real_data_adjusted_p_values.rds")

tmp1 <- which(res_summ_adjust[, 1]<0.01 & res_summ_adjust[, 4]>0.01)

cols <- colnames(geneSetDat[[tmp1[1]]]$dat[[1]])
for (i in 2:length(tmp1)) {
  cols <- c(cols, colnames(geneSetDat[[tmp1[i]]]$dat[[1]]))
}
cols <- unique(cols)

data(mainz)
dat <- Biobase::exprs(mainz)
id <- Biobase::pData(mainz)[, 9]
dat <- limma::normalizeQuantiles(dat, ties=TRUE)
tmp3 <- which(colnames(t(dat))%in%cols)
dat11 <- t(dat)[,tmp3]

# calculate the marginal kurtoses of all genes
tmp4 <- apply(dat11, 2, kurtosis)

## plot the histrogram of kurtoses
pdf(file="real_data_kurtosis_histogram.pdf", width=10, height=8)
par(mar=c(5, 5, 2.5, 2.5), oma=c(1, 1, 0, 0), lwd=4)
hist(tmp4, nclass=20, xaxt="n", yaxt="n", xlab="Kurtosis",
     cex.axis=3, cex.lab=2, cex.main=2.5, main="Histogram of kurtosis")
axis(1, cex.axis=2, lwd=2)
axis(2, cex.axis=2, lwd=2)
abline(v=3, col="red", lwd=2)
dev.off()


################################################################################
######## Part 4: Simultaneous confidence intervals for gene set with GO:0030855
################################################################################

## install package "HDMedian" from the source file "HDMedian_0.1.0.tar.gz"
install.packages("./HDMedian_0.1.0.tar.gz", repos = NULL, type = "source")

library(HDMedian)
library(mvtnorm)
library(Gmedian)
library(extraDistr)
library(ggplot2)
library(ggpubr)
R_seed <- 202405
set.seed(R_seed)

geneSetDat <- readRDS("./geneSetDat.rds")
res_summ_adjust <- readRDS("./real_data_adjusted_p_values.rds")

## get the gene set data corresponding to GO:0030855
i_ind <- 142
res_summ_adjust[i_ind, ]
# median_non_std median_std_T_s median_std_T_s_tilde mean_Lin2021_non_std mean_Lin2021_std 
#     0.08797315     0.00614755           0.00000000           0.15276700       0.07332816
dat_tmp <- geneSetDat[[i_ind]]$dat

#### calculate simultaneous confidence intervals for each pair of groups
## mean-based methd
# group 1 vs group 2
mean12 <- SCI_TwoSample_HD_Mean_Konietschke(dat_tmp[[1]], dat_tmp[[2]], 
                                            nrow(dat_tmp[[1]]),
                                            nrow(dat_tmp[[2]]), 
                                            ncol(dat_tmp[[1]]), level=0.99)
mean12 <- mean12$res_SCI_mean_std

# group 1 vs group 3
mean13 <- SCI_TwoSample_HD_Mean_Konietschke(dat_tmp[[1]], dat_tmp[[3]], 
                                            nrow(dat_tmp[[1]]),
                                            nrow(dat_tmp[[3]]), 
                                            ncol(dat_tmp[[1]]), level=0.99)
mean13 <- mean13$res_SCI_mean_std

# group 2 vs group 3
mean23 <- SCI_TwoSample_HD_Mean_Konietschke(dat_tmp[[2]], dat_tmp[[3]], 
                                            nrow(dat_tmp[[2]]),
                                            nrow(dat_tmp[[3]]), 
                                            ncol(dat_tmp[[2]]), level=0.99)
mean23 <- mean23$res_SCI_mean_std

## median-based method
# group 1 vs group 2
median12 <- SCI_TwoSample_HD_Median(dat_tmp[[1]], dat_tmp[[2]], 
                                    nrow(dat_tmp[[1]]),
                                    nrow(dat_tmp[[2]]), 
                                    ncol(dat_tmp[[1]]), level=0.99)
median12 <- median12$res_SCI_std_ver2

# group 1 vs group 3
median13 <- SCI_TwoSample_HD_Median(dat_tmp[[1]], dat_tmp[[3]], 
                                    nrow(dat_tmp[[1]]),
                                    nrow(dat_tmp[[3]]), 
                                    ncol(dat_tmp[[1]]), level=0.99)
median13 <- median13$res_SCI_std_ver2

# group 2 vs group 3
median23 <- SCI_TwoSample_HD_Median(dat_tmp[[2]], dat_tmp[[3]], 
                                    nrow(dat_tmp[[2]]),
                                    nrow(dat_tmp[[3]]), 
                                    ncol(dat_tmp[[2]]), level=0.99)
median23 <- median23$res_SCI_std_ver2

#### plot the simultaneous confidence intervals
df1 <- data.frame(gene_index=c(1:ncol(dat_tmp[[1]])), 
                  md=colMeans(dat_tmp[[1]])-colMeans(dat_tmp[[2]]),
                  names=rep("SCIs of mean: Group 1 and Group 2", ncol(dat_tmp[[1]])),
                  lower=mean12[1,], upper=mean12[2,])

gp1 <- ggplot(df1, aes(gene_index, md),
              colour=factor(names))+
  theme_grey(base_size = 18) +
  geom_point(size=1.5) +
  facet_wrap(~names, ncol=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=0.75)+
  labs(x=("Gene Index"),y=("Mean Difference"))+
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1)
gp1

df2 <- data.frame(gene_index=c(1:ncol(dat_tmp[[1]])), 
                  md=colMeans(dat_tmp[[1]])-colMeans(dat_tmp[[3]]),
                  names=rep("SCIs of mean: Group 1 and Group 3", ncol(dat_tmp[[1]])),
                  lower=mean13[1,], upper=mean13[2,])

gp2 <- ggplot(df2, aes(gene_index, md),
              colour=factor(names))+
  theme_grey(base_size = 18) +
  geom_point(size=1.5) +
  facet_wrap(~names, ncol=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=0.75)+
  labs(x=("Gene Index"),y=("Mean Difference"))+
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1)
gp2

df3 <- data.frame(gene_index=c(1:ncol(dat_tmp[[2]])), 
                  md=colMeans(dat_tmp[[2]])-colMeans(dat_tmp[[3]]),
                  names=rep("SCIs of mean: Group 2 and Group 3", ncol(dat_tmp[[2]])),
                  lower=mean23[1,], upper=mean23[2,])

gp3 <- ggplot(df3, aes(gene_index, md),
              colour=factor(names))+
  theme_grey(base_size = 18) +
  geom_point(size=1.5) +
  facet_wrap(~names, ncol=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=0.75)+
  labs(x=("Gene Index"),y=("Mean Difference"))+
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1)
gp3

df4 <- data.frame(gene_index=c(1:ncol(dat_tmp[[1]])),
                  md=as.numeric(Weiszfeld(dat_tmp[[1]])$median-Weiszfeld(dat_tmp[[2]])$median),
                  names=rep("SCIs of median: Group 1 and Group 2", ncol(dat_tmp[[1]])),
                  lower=median12[1,], upper=median12[2,])
tmp_upper_index4 <- which(df4$upper<0) # NULL
gp4 <- ggplot(df4, aes(gene_index, md),
              colour=factor(names))+
  theme_grey(base_size = 18) +
  geom_point(size=1.5) +
  facet_wrap(~names, ncol=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=0.75)+
  labs(x=("Gene Index"),y=("Median Difference"))+
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1)
gp4

df5 <- data.frame(gene_index=c(1:ncol(dat_tmp[[1]])),
                  md=as.numeric(Weiszfeld(dat_tmp[[1]])$median-Weiszfeld(dat_tmp[[3]])$median),
                  names=rep("SCIs of median: Group 1 and Group 3", ncol(dat_tmp[[1]])),
                  lower=median13[1,], upper=median13[2,])
upper_index5 <- which(df5$upper<0) # 12 13 14 15 21
gp5 <- ggplot(df5, aes(gene_index, md),
              colour=factor(names))+
  theme_grey(base_size = 18) +
  geom_point(size=1.5) +
  facet_wrap(~names, ncol=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=0.75)+
  labs(x=("Gene Index"),y=("Median Difference"))+
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1) +
  geom_segment(aes(x = upper_index5[1], y = lower[upper_index5[1]],
                   xend = upper_index5[1], yend = upper[upper_index5[1]]),
               color = "red") +
  geom_segment(aes(x = upper_index5[2], y = lower[upper_index5[2]],
                   xend = upper_index5[2], yend = upper[upper_index5[2]]),
               color = "red") +
  geom_segment(aes(x = upper_index5[3], y = lower[upper_index5[3]],
                   xend = upper_index5[3], yend = upper[upper_index5[3]]),
               color = "red") +
  geom_segment(aes(x = upper_index5[4], y = lower[upper_index5[4]],
                   xend = upper_index5[4], yend = upper[upper_index5[4]]),
               color = "red") +
  geom_segment(aes(x = upper_index5[5], y = lower[upper_index5[5]],
                   xend = upper_index5[5], yend = upper[upper_index5[5]]),
               color = "red")

gp5

df6 <- data.frame(gene_index=c(1:ncol(dat_tmp[[2]])),
                  md=as.numeric(Weiszfeld(dat_tmp[[2]])$median-Weiszfeld(dat_tmp[[3]])$median),
                  names=rep("SCIs of median: Group 2 and Group 3", ncol(dat_tmp[[2]])),
                  lower=median23[1,], upper=median23[2,])
tmp_upper_index6 <- which(df6$upper<0) # 12 14
gp6 <- ggplot(df6, aes(gene_index, md),
              colour=factor(names))+
  theme_grey(base_size = 18) +
  geom_point(size=1.5) +
  facet_wrap(~names, ncol=1)+
  geom_errorbar(aes(ymin = lower, ymax = upper), size=0.75)+
  labs(x=("Gene Index"),y=("Median Difference"))+
  geom_hline(yintercept=0, linetype="dashed", color = "blue", size=1) +
  geom_segment(aes(x = tmp_upper_index6[1], y = lower[tmp_upper_index6[1]],
                   xend = tmp_upper_index6[1], yend = upper[tmp_upper_index6[1]]),
               color = "red") +
  geom_segment(aes(x = tmp_upper_index6[2], y = lower[tmp_upper_index6[2]],
                   xend = tmp_upper_index6[2], yend = upper[tmp_upper_index6[2]]),
               color = "red")
gp6

## plot the figures into one pdf file
pdf(file = paste0("real_data_SCIs.pdf"),
    width = 12, # The width of the plot in inches
    height = 15) # The height of the plot in inches
ggarrange(gp1, gp4, gp2, gp5, gp3, gp6, ncol=2, nrow=3)
dev.off()

