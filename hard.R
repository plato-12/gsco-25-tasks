# Gamma mixing kernel with conjugate prior

# Likelihood function for Gamma distribution
Likelihood.gamma <- function(mdobj, x, theta) {
  shape <- theta[[1]]
  rate <- theta[[2]]
  return(as.numeric(dgamma(x, shape = shape, rate = rate)))
}

# Prior draw function for Gamma parameters
# Using a conjugate structure with:
# shape ~ transformed Beta
# rate ~ Gamma
PriorDraw.gamma <- function(mdobj, n) {
  #  for shape
  a0 <- mdobj$priorParameters[1]
  b0 <- mdobj$priorParameters[2]
  #  for rate
  c0 <- mdobj$priorParameters[3]
  d0 <- mdobj$priorParameters[4]

  # Draw shape parameter (using a transform of Beta)
  shape_scale <- mdobj$priorParameters[5]  # Scale factor for shape
  u <- rbeta(n, a0, b0)
  shape <- -shape_scale * log(u)  # Transform to get positive values with appropriate scale

  # Draw rate from Gamma
  rate <- rgamma(n, shape = c0, rate = d0)

  theta <- list(
    array(shape, dim = c(1, 1, n)),
    array(rate, dim = c(1, 1, n))
  )
  return(theta)
}

# Posterior draw function using sufficient statistics and conjugacy
PosteriorDraw.gamma <- function(mdobj, x, n = 1) {

  a0 <- mdobj$priorParameters[1]
  b0 <- mdobj$priorParameters[2]
  c0 <- mdobj$priorParameters[3]
  d0 <- mdobj$priorParameters[4]
  shape_scale <- mdobj$priorParameters[5]

  n_obs <- length(x)

  # Sufficient statistics
  sum_x <- sum(x)
  sum_log_x <- sum(log(x))

  # For the gamma shape parameter, we need to use MCMC since
  # there's no direct conjugate update
  shape_samples <- numeric(n)
  rate_samples <- numeric(n)

  for(i in 1:n) {
    # Initial values
    alpha_current <- 1  # Starting value for shape

    # Since we know the posterior mode for rate given shape
    # we can use Gibbs sampling
    n_iter <- 500
    n_burnin <- 300

    for(j in 1:n_iter) {
      # Update rate given shape (conjugate update)
      beta_posterior <- rgamma(1, shape = c0 + n_obs * alpha_current,
                               rate = d0 + sum_x)

      # Update shape using Metropolis-Hastings
      alpha_proposal <- exp(rnorm(1, log(alpha_current), 0.1))

      # Log posterior ratio for shape
      log_ratio <- (alpha_proposal - alpha_current) * sum_log_x -
        n_obs * (lgamma(alpha_proposal) - lgamma(alpha_current)) +
        n_obs * alpha_proposal * log(beta_posterior) -
        n_obs * alpha_current * log(beta_posterior) +
        (a0 - 1) * (log(1 - exp(-alpha_proposal/shape_scale)) -
                      log(1 - exp(-alpha_current/shape_scale))) +
        (alpha_current - alpha_proposal)/shape_scale

      # Accept/reject shape proposal
      if(log(runif(1)) < log_ratio) {
        alpha_current <- alpha_proposal
      }

      # Store samples after burn-in
      if(j > n_burnin) {
        shape_samples[i] <- alpha_current
        rate_samples[i] <- beta_posterior
      }
    }
  }

  theta <- list(
    array(shape_samples, dim = c(1, 1, n)),
    array(rate_samples, dim = c(1, 1, n))
  )

  return(theta)
}

# Predictive function for Gamma
Predictive.gamma <- function(mdobj, x) {
  pred <- numeric(length(x))
  n_samples <- 1000

  # Draw samples from the prior
  prior_samples <- PriorDraw.gamma(mdobj, n_samples)
  shape_samples <- prior_samples[[1]][1, 1, ]
  rate_samples <- prior_samples[[2]][1, 1, ]

  # For each x, compute the predictive density
  for(i in seq_along(x)) {
    if(x[i] <= 0) {
      pred[i] <- 0  # Gamma density is 0 for x ≤ 0
    } else {
      # Compute Gamma density for each prior sample
      densities <- dgamma(x[i], shape = shape_samples, rate = rate_samples)
      # Average over samples
      pred[i] <- mean(densities)
    }
  }

  return(pred)
}

# Creating the Gamma mixing distribution
gammaMd <- MixingDistribution(distribution = "gamma",
                              priorParameters = c(2, 2, 2, 0.5, 5),  # a0, b0, c0, d0, shape_scale
                              conjugate = "conjugate")

# simulated data from a mixture of two Gamma distributions
set.seed(123)
y <- c(rgamma(200, shape = 2, rate = 1),   # Shape=2, Mean=2
       rgamma(300, shape = 5, rate = 0.5)) # Shape=5, Mean=10

dp <- DirichletProcessCreate(y, gammaMd)
dp <- Initialise(dp)
dp <- Fit(dp, 1000)  # Run for 1000 iterations

# posterior predictive frame for plotting
x_range <- seq(0.1, 25, length.out = 100)
pf <- PosteriorFrame(dp, x_range, 1000)

# data frame for true distribution
trueFrame <- data.frame(
  x = x_range,
  y = 0.4 * dgamma(x_range, shape = 2, rate = 1) +
    0.6 * dgamma(x_range, shape = 5, rate = 0.5)
)

# Plot
ggplot() +
  geom_ribbon(data = pf,
              aes(x = x, ymin = X5., ymax = X95.),
              colour = NA,
              fill = "orange",
              alpha = 0.2) +  # credible intervals
  geom_line(data = pf, aes(x = x, y = Mean), colour = "orange") +  # mean
  geom_line(data = trueFrame, aes(x = x, y = y), linetype = "dashed") +  # true
  labs(title = "Gamma Mixture Model Fit",
       x = "Value",
       y = "Density") +
  theme_minimal()

# Histogram with fitted density
ggplot() +
  geom_histogram(aes(x = y, y = ..density..), bins = 30, fill = "lightgray", color = "black") +
  geom_line(data = pf, aes(x = x, y = Mean), colour = "orange", size = 1) +
  geom_line(data = trueFrame, aes(x = x, y = y), linetype = "dashed") +
  labs(title = "Gamma Mixture Model Fit with Data Histogram",
       x = "Value",
       y = "Density") +
  theme_minimal()
