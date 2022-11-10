######### R script: logistic_regn.R ##########

# For performing Bayesian logistic regression.

# Created: 28 SEP 2022
# Last changed: 28 SEP 2022

# Load libraries:

library(MASS)
library(magic)
library(lattice)
library(pracma)
library(ellipse)

# Required functions:

setwd("functions")

source("fpca_algs.r")
source("X_design.r")
source("ZOSull.r")
source("OmegaOSull.r")
source("vec.r")
source("vecInverse.r")
source("tr.r")
source("trapint.r")
source("cprod.r")
source("wait.r")
source("fourier_basis.r")
source("ise.r")
source("stan_mods.r")
source("vmp_functions.r")

setwd("..")

# Establish simulation variables:

n <- 100                            # number of interior knots
data_col <- "blue"                  # colour of the data in the plots
criterion <- 1e-6                   # convergence criterion

n_vmp <- 100                        # number of VMP iterations
n_g <- 1000                         # length of the plotting grid
vmp_col <- "black"                  # colour of the VMP plots

beta <- c(-2.2, 3.8, 1.3)           # beta coefficient vector
d <- length(beta)                   # dimension of spline vector

plot_dim <- c(1, d)

# Simulate the data:

x <- matrix(NA, n, d - 1)
for(j in 1:(d - 1)) {
	
	x[, j] <- rnorm(n, 0, 1/j)
}
X <- X_design(x)
X_beta <- as.vector(X %*% beta)
prob <- 1/(1 + exp(-X_beta))
y <- rbinom(n, 1, prob)

# Establish hyperparameters:

sigsq_beta <- 1e10

# VMP simulations:

eta <- vector("list", length = 4)
names(eta) <- c(
	"beta->p(Y|beta)", "p(Y|beta)->beta",
	"beta->p(beta)", "p(beta)->beta"
)

eta$"p(beta)->beta" <- gauss_prior_frag(rep(0, d), sigsq_beta*diag(d))
eta$"p(Y|beta)->beta" <- gauss_prior_frag(rep(0, d), diag(d)) - eta$"p(beta)->beta"

for(iter in 1:n_vmp) {
	
	cat("starting iteration", iter, "of", n_vmp, "\n")
	
	eta$"beta->p(Y|beta)" <- eta$"p(beta)->beta"
	eta$"beta->p(beta)" <- eta$"p(Y|beta)->beta"
	
	# Update p(Y|beta) fragment:
	
	eta_in <- list(
		eta$"beta->p(Y|beta)",
		eta$"p(Y|beta)->beta"
	)
	
	logistic_lik_fragment <- logistic_lik_frag(eta_in, y, X)
	
	eta$"p(Y|beta)->beta" <- logistic_lik_fragment$"eta"[[1]]
}

# Save the original q_beta:

eta_beta <- list(eta$"p(Y|beta)->beta", eta$"p(beta)->beta")
q_beta <- gauss_q(eta_beta, use_vech = TRUE)
E_q_beta <- q_beta[[1]]
Cov_q_beta <- q_beta[[2]]

par(mfrow = plot_dim)

for(j in 1:d) {
	
	x_g_min <- E_q_beta[j] + qnorm(0.9999)*Cov_q_beta[j, j]
	x_g_max <- E_q_beta[j] + qnorm(0.0001)*Cov_q_beta[j, j]
	x_g <- seq(x_g_min, x_g_max, length = n_g)
	y_g <- dnorm(x_g, mean = E_q_beta[j], sd = sqrt(Cov_q_beta[j, j]))
	plot(x_g, y_g, col = "blue", type = "l")
	abline(v = beta[j], col = "red")
}









