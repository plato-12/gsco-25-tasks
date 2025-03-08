library(dirichletprocess)
library(ggplot2)
library(coda)

# Generate data from lognormal mixture
set.seed(123)
n <- 500
component1 <- rlnorm(n * 0.7, meanlog = 0, sdlog = 0.5)  # 70% from component1
component2 <- rlnorm(n * 0.3, meanlog = 1.5, sdlog = 0.3)  # 30% from component2
lnorm_data <- c(component1, component2)

# Histogram
hist(lnorm_data, breaks = 30, freq = FALSE,
     main = "Lognormal Mixture Data", xlab = "Value")
lines(density(lnorm_data), col = "blue", lwd = 2)

# Scaling
lnorm_data_scaled <- scale(lnorm_data)

# Different models with different alpha priors
# 1. Default prior: alphaPriors = c(2, 4)
dp_lnorm1 <- DirichletProcessGaussian(lnorm_data_scaled)
dp_lnorm1 <- Fit(dp_lnorm1, 2000)

# 2. Prior encouraging fewer clusters
dp_lnorm2 <- DirichletProcessGaussian(lnorm_data_scaled, alphaPriors = c(1, 2))
dp_lnorm2 <- Fit(dp_lnorm2, 2000)

# 3. Prior encouraging more clusters
dp_lnorm3 <- DirichletProcessGaussian(lnorm_data_scaled, alphaPriors = c(5, 1))
dp_lnorm3 <- Fit(dp_lnorm3, 2000)

# Posterior density estimates
x_grid <- seq(min(lnorm_data_scaled), max(lnorm_data_scaled), length.out = 1000)
posterior_frame <- PosteriorFrame(dp_lnorm1, x_grid, ndraws=1000)

# Generating actual posterior samples
n_samples <- 5000
cluster_probs <- dp_lnorm1$pointsPerCluster / sum(dp_lnorm1$pointsPerCluster)
selected_clusters <- sample(1:dp_lnorm1$numberClusters, n_samples, replace=TRUE, prob=cluster_probs)

# For each selected cluster, generating a sample from the Gaussian
# with that cluster's parameters
posterior_samples <- numeric(n_samples)
for (i in 1:n_samples) {
  cluster <- selected_clusters[i]
  mu <- dp_lnorm1$clusterParameters[[1]][,,cluster]
  sigma <- dp_lnorm1$clusterParameters[[2]][,,cluster]
  posterior_samples[i] <- rnorm(1, mu, sigma)
}

# Transform back to original scale
posterior_samples_original <- posterior_samples * sd(lnorm_data) + mean(lnorm_data)

# Calculating quantiles from the generated samples
quantiles <- quantile(posterior_samples_original, c(0.05, 0.95))
cat("5% and 95% quantiles of posterior samples:", quantiles[1], "and", quantiles[2], "\n")

# Plot the posterior samples and density
ggplot() +
  # Original data density
  geom_density(data = data.frame(x = lnorm_data), aes(x = x),
               color = "blue", linetype = 2) +
  # Posterior samples density
  geom_density(data = data.frame(x = posterior_samples_original),
               aes(x = x), color = "red") +
  # Histogram of posterior samples
  geom_histogram(data = data.frame(x = posterior_samples_original),
                 aes(x = x, y = after_stat(density)), bins = 30,
                 alpha = 0.2, fill = "black") +
  labs(title = "Posterior Samples Distribution",
       x = "Value", y = "Density") +
  theme_minimal() +
  annotate("text", x = Inf, y = Inf,
           label = paste("5% quantile:", round(quantiles[1], 2),
                         "\n95% quantile:", round(quantiles[2], 2)),
           hjust = 1.1, vjust = 1.5, size = 4)

# Comparing number of clusters for each model
n_clusters1 <- sapply(dp_lnorm1$labelsChain, function(x) length(unique(x)))
n_clusters2 <- sapply(dp_lnorm2$labelsChain, function(x) length(unique(x)))
n_clusters3 <- sapply(dp_lnorm3$labelsChain, function(x) length(unique(x)))

cluster_df <- data.frame(
  Iteration = rep(1:length(n_clusters1), 3),
  Clusters = c(n_clusters1, n_clusters2, n_clusters3),
  Prior = rep(c("Default (2,4)", "Few clusters (1,2)", "Many clusters (5,1)"),
              each = length(n_clusters1))
)

# Plot number of clusters
ggplot(cluster_df, aes(x = Iteration, y = Clusters, color = Prior)) +
  geom_line() +
  labs(title = "Number of Clusters by Prior",
       x = "MCMC Iteration", y = "Number of Clusters") +
  theme_minimal()

# Plot alpha chains to check convergence
alpha_df <- data.frame(
  Iteration = rep(1:length(dp_lnorm1$alphaChain), 3),
  Alpha = c(dp_lnorm1$alphaChain, dp_lnorm2$alphaChain, dp_lnorm3$alphaChain),
  Prior = rep(c("Default (2,4)", "Few clusters (1,2)", "Many clusters (5,1)"),
              each = length(dp_lnorm1$alphaChain))
)

# Plot alpha
ggplot(alpha_df, aes(x = Iteration, y = Alpha, color = Prior)) +
  geom_line() +
  labs(title = "Alpha Parameter Convergence by Prior",
       x = "MCMC Iteration", y = "Alpha") +
  theme_minimal()

# Average number of clusters in the second half of the chain
avg_clusters1 <- mean(n_clusters1[1000:2000])
avg_clusters2 <- mean(n_clusters2[1000:2000])
avg_clusters3 <- mean(n_clusters3[1000:2000])

cat("Average number of clusters (second half of chain):\n")
cat("Default prior (2,4):", avg_clusters1, "\n")
cat("Prior for fewer clusters (1,2):", avg_clusters2, "\n")
cat("Prior for more clusters (5,1):", avg_clusters3, "\n")

# Calculate convergence metrics
burnin <- 1000
effective_samples1 <- coda::effectiveSize(dp_lnorm1$alphaChain[burnin:2000])
effective_samples2 <- coda::effectiveSize(dp_lnorm2$alphaChain[burnin:2000])
effective_samples3 <- coda::effectiveSize(dp_lnorm3$alphaChain[burnin:2000])

cat("Effective sample sizes after burn-in:\n")
cat("Default prior:", effective_samples1, "\n")
cat("Prior for fewer clusters:", effective_samples2, "\n")
cat("Prior for more clusters:", effective_samples3, "\n")

