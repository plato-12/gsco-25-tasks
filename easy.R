library(dirichletprocess)
library(ggplot2)

# PART 1: Univariate model with faithful dataset
data(faithful)

faithful_scaled <- scale(faithful$waiting)

set.seed(123)
dp_faithful <- DirichletProcessGaussian(faithful_scaled)
dp_faithful <- Fit(dp_faithful, 1000)

plot(dp_faithful)


# PART 2: Multivariate model with iris dataset
data(iris)

iris_scaled <- scale(iris[, 1:4])

set.seed(456)
dp_iris <- DirichletProcessMvnormal(iris_scaled)
dp_iris <- Fit(dp_iris, 1000)

plot(dp_iris)

# number of clusters
unique_clusters <- length(unique(dp_iris$clusterLabels))
cat("Number of clusters in iris data:", unique_clusters, "\n")

# cluster sizes
table(dp_iris$clusterLabels)
