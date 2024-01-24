### Simulate angular central Gaussian distribution vectors ---------------------
library(rotasym) # https://github.com/egarpor/rotasym
# Package citation for rotasym:
# García-Portugués et al. (2020) Journal of American Statistical Association
# URL: https://www.tandfonline.com/doi/abs/10.1080/01621459.2019.1665527

## Provide singular values / standard deviations 
simACG <- function(s_vec, 
                   n_obs,
                   n_boot) {
  # Assert that all entries of s_vec are positive
  assertthat::assert_that(all(s_vec > 0),
                          msg="All entries of s_vec must be positive.")
  
  # Provide messages
  message(date(), ": Dimension of angular vectors = ", length(s_vec))
  message(date(), ": No. observations, n, = ", n_obs)
  message(date(), ": No. bootstrap replicates, n_boot, = ", n_boot)
  
  # Generate n_obs spherical/angular vectors
  # Square singular values / standard deviations to get eigenvalues / variances
  x <- r_ACG(n_obs, Lambda=diag(s_vec^2)) 
  
  # Compute empirical mean and variance of simulated vectors
  comp_means <- colMeans(x)
  comp_vars <- apply(x,2,function(z) {mean(z^2)})
  comp_vector <- c(comp_means,comp_vars)
  
  # Get bootstrap errors 
  oneBootstrap <- function(df) {
    df_resampled <- df[sample(nrow(df),size=nrow(df),replace=TRUE),]
    resampled_comp_means <- colMeans(df_resampled)
    resampled_comp_vars <- apply(df_resampled,2,function(z) {mean(z^2)})
    return(c(resampled_comp_means,resampled_comp_vars))
  }
  
  booted_statistics <- replicate(n_boot, oneBootstrap(x), simplify='array')
  
  # Get bootstrap CIs for each statistic
  boot_upper_vals <- c()
  boot_lower_vals <- c()
  
  for (j in 1:nrow(booted_statistics)) {
    # get vector of bootstrapped statistics
    bootstrapped_vector <- booted_statistics[j,]
    
    # sort from lowest to highest
    lower_and_upper_cutoffs <- quantile(bootstrapped_vector,probs=c(0.025,0.975))
    lower_cutoff <- lower_and_upper_cutoffs[1]; upper_cutoff <- lower_and_upper_cutoffs[2]
    
    # compute upper and lower bound, add to vector
    boot_upper <- 2*comp_vector[j] - lower_cutoff
    boot_lower <- 2*comp_vector[j] - upper_cutoff
    boot_upper_vals <- c(boot_upper_vals, boot_upper); names(boot_upper_vals) <- NULL
    boot_lower_vals <- c(boot_lower_vals, boot_lower); names(boot_lower_vals) <- NULL
  }
  
  # Return dataframes
  mean_df <- data.frame(Analytical=rep(0,length(s_vec)),
                        Empirical=comp_means,
                        Upper=boot_upper_vals[1:length(s_vec)],
                        Lower=boot_lower_vals[1:length(s_vec)])
  var_df <- data.frame(Analytical=s_vec/sum(s_vec),
                       Empirical=comp_vars,
                       Upper=boot_upper_vals[(length(s_vec)+1):(2*length(s_vec))],
                       Lower=boot_lower_vals[(length(s_vec)+1):(2*length(s_vec))])
  return(list(MEAN = mean_df,
              VARIANCE = var_df))
}

### Main Body ------------------------------------------------------------------
## Test on bivariate random vector
set.seed(2023)
test_s_vec <- c(2,1)
sim_result <- simACG(s_vec=test_s_vec,
                     n_obs=1e4,
                     n_boot=1e3)
message(date(), ": Finished.")
sim_result$MEAN
sim_result$VARIANCE

# > sim_result$MEAN
# ANALYTICAL    EMPIRICAL
# 1          0 0.0007460246
# 2          0 0.0006563815
# > sim_result$VARIANCE
# ANALYTICAL EMPIRICAL
# 1  0.6666667 0.6660469
# 2  0.3333333 0.3339531

melted_df <- sim_result$VARIANCE %>% 
  mutate(COMPONENT=paste0("S",1:nrow(sim_result$VARIANCE))) %>%
  melt(id.vars="COMPONENT")
ggplot(melted_df,aes(y=value,x=factor(COMPONENT,levels=paste0("S",1:nrow(sim_result$VARIANCE))))) +
  geom_point(aes(shape=variable)) +
  geom_line(aes(group=variable)) +
  theme_bw() +
  xlab('Subspace associated with singular vector') +
  ylab('Squared magnitude of cosine similarity\n of singular vector with random PRS')

## Test on 3-variate random vector
test_s_vec <- c(3,2,1)
sim_result <- simACG(s_vec=test_s_vec,
                     n_obs=1e6)
sim_result$MEAN
sim_result$VARIANCE

# > sim_result$MEAN
# THEORY     EMPIRICAL
# 1      0 -1.864713e-04
# 2      0 -7.080303e-04
# 3      0  6.434021e-05
# > sim_result$VARIANCE
# THEORY EMPIRICAL
# 1 0.5000000 0.5257665
# 2 0.3333333 0.3378968
# 3 0.1666667 0.1363367

melted_df <- sim_result$VARIANCE %>% 
  mutate(COMPONENT=paste0("S",1:nrow(sim_result$VARIANCE))) %>%
  melt(id.vars="COMPONENT")
ggplot(melted_df,aes(y=value,x=factor(COMPONENT,levels=paste0("S",1:nrow(sim_result$VARIANCE))))) +
  geom_point(aes(shape=variable)) +
  geom_line(aes(group=variable)) +
  theme_bw() +
  xlab('Subspace associated with singular vector') +
  ylab('Squared magnitude of cosine similarity\n of singular vector with random PRS')


## Test on Galinsky eigenvalues
# Data citation Galinsky eigenvalues:
# Galinsky et al. (2016) American Journal of Human Genetics
# URL: https://www.cell.com/ajhg/fulltext/S0002-9297(16)30395-0#supplementaryMaterial
set.seed(2023)
Galinsky2016TableS2Eigenvals <- c(20.99,9.35,7.76,5.18,
                                  5.13,4.62,4.61,4.59,
                                  4.59,4.57)
test_s_vec <- rep(1,10)
sim_result_null <- simACG(s_vec=test_s_vec,
                          n_obs=1e4,
                          n_boot=1e3)
sim_result_null$MEAN
sim_result_null$VARIANCE

input_s_vec <- sqrt(Galinsky2016TableS2Eigenvals)
sim_result <- simACG(s_vec=input_s_vec,
                     n_obs=1e4,
                     n_boot=1e3)
sim_result$MEAN
sim_result$VARIANCE

# > sim_result$MEAN
# ANALYTICAL     EMPIRICAL
# 1           0  7.086913e-04
# 2           0 -1.805030e-04
# 3           0  1.656087e-04
# 4           0 -1.820976e-04
# 5           0 -4.586059e-06
# 6           0  2.827231e-05
# 7           0 -2.444250e-04
# 8           0 -1.343253e-04
# 9           0 -1.692531e-04
# 10          0 -6.657755e-04
# > sim_result$VARIANCE
# ANALYTICAL  EMPIRICAL
# 1  0.17837217 0.23873777
# 2  0.11904926 0.12998854
# 3  0.10845554 0.11199261
# 4  0.08861062 0.08000290
# 5  0.08818192 0.07917168
# 6  0.08368389 0.07216144
# 7  0.08359328 0.07236097
# 8  0.08341175 0.07206844
# 9  0.08341175 0.07192567
# 10 0.08322983 0.07158999

null_melted_df <- sim_result_null$VARIANCE %>%
  mutate(COMPONENT=paste0("S",1:nrow(sim_result$VARIANCE))) %>%
  melt(id.vars="COMPONENT")

plot_A <- ggplot(null_melted_df,aes(y=value,x=factor(COMPONENT,levels=paste0("S",1:nrow(sim_result$VARIANCE))))) +
  geom_point(aes(shape=variable)) +
  geom_line(aes(group=variable)) +
  theme_bw() +
  xlab('Subspace associated with PC') +
  ylab('Mean squared magnitude of cosine similarity\n of PC with random PGS') +
  theme(legend.pos=c(0.6,0.8)) +
  scale_shape_discrete(name="") +
  #ylim(c(0.0996,0.1004)) +
  ggtitle('A. All Eigenvalues/Singular Values\n    Equal')

test_data_null <- sim_result_null$VARIANCE %>%
  mutate(COMPONENT=paste0("S",1:nrow(sim_result$VARIANCE)))

new_plot_A <- ggplot(test_data_null, aes(x = factor(COMPONENT,levels=paste0("S",1:nrow(sim_result$VARIANCE))))) +
  geom_line(aes(y = Analytical), 
            color = "grey", 
            lty='dashed',
            group = 1) +
  geom_point(aes(y = Analytical), 
             color = "black", 
             shape=1) +
  geom_point(aes(y = Empirical), color = "red", size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  xlab('Subspace associated with PC') +
  ylab('Mean squared magnitude of cosine similarity\n of PC with random PGS') +
  theme_bw() +
  ggtitle('A. All Eigenvalues/Singular Values\n    Equal') +
  theme(plot.title=element_text(face='bold'))

melted_df <- sim_result$VARIANCE %>% 
  mutate(COMPONENT=paste0("S",1:nrow(sim_result$VARIANCE))) %>%
  melt(id.vars="COMPONENT")

plot_B <- ggplot(melted_df,aes(y=value,x=factor(COMPONENT,levels=paste0("S",1:nrow(sim_result$VARIANCE))))) +
  geom_point(aes(shape=variable)) +
  geom_line(aes(group=variable)) +
  theme_bw() +
  xlab('Subspace associated with PC') +
  ylab('Mean squared magnitude of cosine similarity\n of PC with random PGS') +
  ylim(c(0.05,0.25)) +
  theme(legend.pos=c(0.6,0.8)) +
  scale_shape_discrete(name="") +
  ggtitle('B. Eigenvalues from UKB PCA\n    (Galinsky et al., 2016)') 

test_data <- sim_result$VARIANCE %>%
  mutate(COMPONENT=paste0("S",1:nrow(sim_result$VARIANCE)))

new_plot_B <- ggplot(test_data, aes(x = factor(COMPONENT,levels=paste0("S",1:nrow(sim_result$VARIANCE))))) +
  geom_line(aes(y = Analytical), 
            color = "grey", 
            lty='dashed',
            group = 1) +
  geom_point(aes(y = Analytical), 
             color = "black", 
             shape=1) +
  geom_point(aes(y = Empirical), color = "red", size = 2) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  xlab('Subspace associated with PC') +
  ylab('Mean squared magnitude of cosine similarity\n of PC with random PGS') +
  theme_bw() +
  ggtitle('B. Eigenvalues from UKB PCA\n    (Galinsky et al., 2016)') +
  theme(plot.title=element_text(face='bold'))

combined_fig <-gridExtra::grid.arrange(new_plot_A+xlab(''),new_plot_B,
                                       nrow=2)

# Optional: Save file
ggsave(combined_fig,
       filename = 'new_AngularCentralGaussian_010524.jpg',
       width = 4.5, height = 9,
       dpi = 300)
