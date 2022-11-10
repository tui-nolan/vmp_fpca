######### R script: logistic_flr.R ##########

# For performing functional linear regression via
# vmp with an FPCA decomposition on the functional
# covariates

# Created: 24 SEP 2022
# Last changed: 06 OCT 2022

# Load libraries:

library(MASS)
library(magic)
library(rstan)
rstan_options(auto_write = TRUE)
library(lattice)
library(ellipse)

# Required functions:

setwd("functions")

source("X_design.r")
source("ZOSull.r")
source("OmegaOSull.r")
source("vec.r")
source("vecInverse.r")
source("tr.r")
source("trapint.r")
source("cprod.r")
source("wait.r")
source("stan_mods.r")
source("fourier_basis.r")
source("vmp_functions.r")
source("fpca_algs.r")

setwd("..")

# Establish simulation variables:

N <- 250                             # number of curves
T_vec <- sample(50:70, N, replace = TRUE)    # number of time observations for each curve
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 2                              # guess of the number of FPCA basis functions
data_col <- "black"                 # colour of the data in the plots
criterion <- 1e-6                   # convergence criterion
d <- (K+2)*(L+1)              # dimension of spline vector
alpha <- c(-2.2, 3.8, 1.3)              # logistic coefficients

n_vmp_fpca <- 200                        # number of FPCA VMP iterations
n_vmp_logistic <- 100                # number of logistic regression iterations
n_g <- 1000                         # length of the plotting grid
vmp_col <- "red"                    # colour of the VMP plots
vmp_lwd <- 1                        # line width for vmp plots

sigma_zeta_vec <- 1/(1:L)^2      # sd for the scores
sigma_eps_x <- 1                        # sd of covariate residuals

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
mu_beta <- rep(0, 2)
A <- 1e5
sigsq_zeta <- 1
sigma_zeta <- sqrt(sigsq_zeta)
Sigma_zeta <- sigsq_zeta*diag(L)

sigsq_alpha <- 1e10

plot_dim <- c(1, L + 1)

# Set the mean function and the FPCA basis functions:

mu <- function(t) return(3*sin(pi*t) - 3/2)
Psi <- fourier_basis(L)

coeff_func <- function(t, a, Psi) {
	
	L <- length(Psi)
	Psi <- sapply(Psi, do.call, list(t))
	
	ans <- as.vector(Psi %*% a)
	return(ans)
}

beta_0 <- alpha[1]
beta_1 <- function(t) {
	
	ans <- coeff_func(t, alpha[-1], Psi)
	return(ans)
}

# Generate functional data:

flr_dataset <- flr_data(
	T_vec, K, n_g, sigma_zeta_vec, beta_0,
	sigma_eps_x, mu, Psi, beta_1,
	likelihood = "logistic"
)

time_obs <- flr_dataset$"time_obs"
time_g <- flr_dataset$"time_g"
int_knots <- flr_dataset$"int_knots"
C <- flr_dataset$"C"
C_g <- flr_dataset$"C_g"
zeta <- flr_dataset$"zeta"
mu_g <- flr_dataset$"mu_g"
Psi_g <- flr_dataset$"Psi_g"
beta_1_g <- flr_dataset$"beta_1_g"
x <- flr_dataset$"x"
y <- flr_dataset$"y"

# VMP simulations:

eta_vec <- vmp_gauss_fpca(
	n_vmp_fpca, N, L, C, x, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, plot_elbo = TRUE
)

# Get the posterior estimates

eta_in <- list(
	eta_vec$"p(nu|Sigma_nu)->nu", eta_vec $"p(Y|nu,zeta,sigsq_eps)->nu",
	eta_vec $"p(zeta)->zeta", eta_vec $"p(Y|nu,zeta,sigsq_eps)->zeta"
)
fpc_rotns <- fpc_orthogonalization(eta_in, 1:N, time_g, C_g, Psi_g)

# Summarise the VMP results:

gbl_hat <- fpc_rotns$"gbl_curves"
mu_hat <- gbl_hat[, 1]
Psi_hat <- gbl_hat[, -1]
Zeta_hat <- fpc_rotns$"zeta"

# Logistic VMP simulations:

A <- X_design(Zeta_hat)

eta <- vector("list", length = 4)
names(eta) <- c(
	"alpha->p(Y|alpha)", "p(Y|alpha)->alpha",
	"alpha->p(alpha)", "p(alpha)->alpha"
)

eta$"p(alpha)->alpha" <- gauss_prior_frag(rep(0, L + 1), sigsq_alpha*diag(L + 1))
eta$"p(Y|alpha)->alpha" <- gauss_prior_frag(rep(0, L + 1), diag(L + 1)) - eta$"p(alpha)->alpha"

for(iter in 1:n_vmp_logistic) {
	
	cat("starting iteration", iter, "of", n_vmp_logistic, "\n")
	
	eta$"alpha->p(Y|alpha)" <- eta$"p(alpha)->alpha"
	eta$"alpha->p(alpha)" <- eta$"p(Y|alpha)->alpha"
	
	# Update p(Y|alpha) fragment:
	
	eta_in <- list(
		eta$"alpha->p(Y|alpha)",
		eta$"p(Y|alpha)->alpha"
	)
	
	logistic_lik_fragment <- logistic_lik_frag(eta_in, y, A)
	
	eta$"p(Y|alpha)->alpha" <- logistic_lik_fragment$"eta"[[1]]
}

# Save the original q_alpha:

eta_alpha <- list(eta$"p(Y|alpha)->alpha", eta$"p(alpha)->alpha")
q_alpha <- gauss_q(eta_alpha, use_vech = TRUE)
E_q_alpha <- q_alpha[[1]]
Cov_q_alpha <- q_alpha[[2]]

# Plot the alpha results:

par(mfrow = plot_dim)

for(j in 1:(L + 1)) {
	
	x_g_min <- E_q_alpha[j] + qnorm(0.9999)*Cov_q_alpha[j, j]
	x_g_max <- E_q_alpha[j] + qnorm(0.0001)*Cov_q_alpha[j, j]
	x_g <- seq(x_g_min, x_g_max, length = n_g)
	y_g <- dnorm(x_g, mean = E_q_alpha[j], sd = sqrt(Cov_q_alpha[j, j]))
	plot(x_g, y_g, col = "blue", type = "l")
	abline(v = alpha[j], col = "red")
}

wait()

# Determine the coefficient function:

par(mfrow = c(1, 1))

sd_vec <- sqrt(diag(tcrossprod(Psi_hat %*% Cov_q_alpha[-1, -1], Psi_hat)))
beta_1_hat <- as.vector(Psi_hat %*% E_q_alpha[-1])
beta_1_low <- beta_1_hat + qnorm(0.025)*sd_vec
beta_1_upp <- beta_1_hat + qnorm(0.975)*sd_vec

y_min <- min(beta_1_hat, beta_1_low, beta_1_upp)
y_max <- max(beta_1_hat, beta_1_low, beta_1_upp)
y_lim <- c(y_min, y_max)
plot(time_g, beta_1_hat, col = "blue", type = "l", ylim = y_lim)
lines(time_g, beta_1_g, col = "red")
lines(time_g, beta_1_low, col = "blue", lty = 2)
lines(time_g, beta_1_upp, col = "blue", lty = 2)

acc <- trapint(time_g, (beta_1_hat - beta_1_g)^2)/trapint(time_g, beta_1_g^2)
print(acc)

############ End of logistic_flr.R ############