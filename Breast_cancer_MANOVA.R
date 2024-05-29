##########################################################################
######## MANOVA on gene sets data
##########################################################################
# load functions for MANOVA
source("./HD_MANOVA_GM_functions.R")

# load necessary packages
library(mvtnorm)
library(Gmedian)
library(flare)
library(extraDistr)
library(doRNG)

# load gene sets data
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

## summarize the results to get the p-values of each test
res_summ <- data.frame(matrix(NA, length(real_MANOVA), 5))
names(res_summ) <- c("median_non_std", "median_std_T_s", "median_std_T_s_tilde",
                     "mean_Lin2021_non_std", "mean_Lin2021_std")

for (i in 1:length(real_MANOVA)) {
  res_median <- real_MANOVA[[i]]$res_median
  res_Lin2021 <- real_MANOVA[[i]]$res_Lin2021

  res_summ$median_non_std[i] <- mean(res_median$test_stat_median<res_median$res_median_boot)
  res_summ$median_std_T_s[i] <- mean(res_median$test_stat_median_std<res_median$res_median_std_boot)
  res_summ$median_std_T_s_tilde[i] <- mean(res_median$test_stat_median_std<res_median$res_median_std_boot_old)
  res_summ$mean_Lin2021_non_std[i] <- mean(res_Lin2021$test_stat_mean<res_Lin2021$res_mean_boot)
  res_summ$mean_Lin2021_std[i] <- mean(res_Lin2021$test_stat_mean_std<res_Lin2021$res_mean_std_boot)
}

## adjust the p-values using the Benjamini{Yekutieli procedure
res_summ_adjust <- matrix(0, nrow(res_summ), 5)
for (j in 1:5) {
  res_summ_adjust[, j] <- p.adjust(res_summ[,j], method="BY")
}
saveRDS(res_summ_adjust, "real_data_adjusted_p_values.rds")

## the number of gene sets identified by each method
apply(res_summ_adjust, 2, function(x) sum(x<0.01))
## 69 141 145  64 140


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

# plot the histrogram of kurtoses
pdf(file="real_data_kurtosis_histogram.pdf", width=10, height=8)
par(mar=c(5, 5, 2.5, 2.5), oma=c(1, 1, 0, 0), lwd=4)
hist(tmp4, nclass=20, xaxt="n", yaxt="n", xlab="Kurtosis",
     cex.axis=3, cex.lab=2, cex.main=2.5, main="Histogram of kurtosis")
axis(1, cex.axis=2, lwd=2)
axis(2, cex.axis=2, lwd=2)
abline(v=3, col="red", lwd=2)
dev.off()

