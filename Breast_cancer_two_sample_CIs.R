##########################################################################
######## Simultaneous confidence intervals for gene set with GO:0030855
##########################################################################

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
