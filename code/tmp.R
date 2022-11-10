######### R script: tmp.R ##########

# For performing a simple functional principal
# components analysis via MCMC to determine
# the FPCA basis functions. The FPCA basis functions
# are determined via semiparametric regression.

# Created: 17 AUG 2020
# Last changed: 19 AUG 2022

# Load libraries:

library(MASS)
library(magic)
library(rstan)
rstan_options(auto_write = TRUE)
library(lattice)
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
source("vmp_functions.r")
source("stan_mods.r")
source("ise.r")
source("fourier_basis.r")

setwd("..")

N_vec <- c(10, 50, 100, 250, 500)       # number of curves
n_sims <- 20

M_min <- 10
M_max <- 15

m_sample <- M_min                          # number of second level curves for the plots

M_sample <- 1:M_min   # specific second level curves for the plots

criterion  <- 1e-5                        # convergence criterion

n_int_knots <- 10                         # number of interior knots
K <- n_int_knots + 2                      # number of spline basis functions
L_1 <- 3                                  # number of first level eigenfunctions
L_2 <- 3                                  # number of second level eigenfunctions
L <- L_1 + L_2

n_vmp <- 200                              # number of VMP iterations
n_g <- 1000                               # length of the plotting grid
vmp_col <- "red"                 # colour of the VMP plots
d <- (K+2)*(L_1 + L_2 + 1)                # dimension of spline vector

n_burnin <- 1000                   # Length of burn-in.
n_mcmc <- 1000                     # Size of the kept sample.
n_thin <- 1                        # Thinning factor. 
tolerance <- 1e-10

sigma_zeta_1 <- 1/(1:L_1)                 # vector of st. dev.'s for the first level scores
sigma_zeta_2 <- 1/(1:L_2)                 # vector of st. dev.'s for the second level scores

sigma_eps <- 1                            # sd of the residuals
sigsq_eps <- sigma_eps^2

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
mu_beta <- rep(0, 2)
A <- 1e5
sigma_zeta <- 1
sigsq_zeta <- sigma_zeta^2

# Set the mean function and the FPCA basis functions:

mu <- function(t) return(10*sin(pi*t) - 5)
Psi_1 <- fourier_basis(L_1)
Psi_2 <- fourier_basis(L_1 + L_2)[(L_1 + 1):(L_1 + L_2)]

# Set up txt files:

psi_1_names <- rep(NA, L_1)
for(l in 1:L_1) {
	
	psi_1_names[l] <- paste("psi_1", l, sep="")
}

psi_2_names <- rep(NA, L_2)
for(l in 1:L_2) {
	
	psi_2_names[l] <- paste("psi_2", l, sep="")
}

col_names <- c("N", "sim", "mu", psi_1_names, psi_2_names)
n_col <- length(col_names)

write(col_names, "res/tmp_mlfpca.txt", ncol = n_col, append = FALSE)
write(col_names, "res/tmp_gauss_mlfpca.txt", ncol = n_col, append = FALSE)

for(N in N_vec) {
	
	n_sample <- N                             # number of curves for the plots
	M <- sample(M_min:M_max, N, replace = TRUE)     # number of observations on each subject
	N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots
	
	for(sim in 1:n_sims) {
		
		cat("N = ", N, "; simulation ", sim, "\n", sep = "")
		
		set.seed(sim)
		
		T_vec <- vector("list", length = N)       # number of time observations for each curve
		for(i in 1:N) {
			
			T_vec[[i]] <- round(runif(M[i], 20, 30))
		}
		
		mlfpca_data <- gauss_mlfpca_data(
			T_vec, K, n_g, sigma_zeta_1, sigma_zeta_2,
			sigma_eps, mu_func, Psi_1, Psi_2
		)
		
		time_obs <- mlfpca_data$time_obs
		time_g <- mlfpca_data$time_g
		int_knots <- mlfpca_data$int_knots
		X <- mlfpca_data$X
		Z <- mlfpca_data$Z
		C <- mlfpca_data$C
		C_g <- mlfpca_data$C_g
		zeta_1 <- mlfpca_data$zeta_1
		zeta_2 <- mlfpca_data$zeta_2
		mu_g <- mlfpca_data$mu_g
		Psi_1_g <- mlfpca_data$Psi_1_g
		Psi_2_g <- mlfpca_data$Psi_2_g
		Y <- mlfpca_data$Y
		
		# Run vmp_mlfpca algorithm:
		
		set.seed(1)
		
		mlfpca_alg <- vmp_mlfpca(
			n_vmp, N, M, L_1, L_2, C, Y,
			sigsq_beta, A, criterion, plot_elbo = FALSE
		)
		
		eta_vec <- mlfpca_alg$"eta_vec"
		
		eta_in <- list(
			eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(zeta)->zeta",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		)
		
		mlfpca_res <- mlfpc_rotation(
			eta_in, time_g, C_g, L_1, L_2,
			N_sample, M_sample, Psi_g = list(Psi_1_g, Psi_2_g)
		)
		
		gbl_hat <- mlfpca_res$"gbl_curves"
		mu_hat <- gbl_hat[, 1]
		Psi_1_hat <- gbl_hat[, 1:L_1 + 1]
		Psi_2_hat <- gbl_hat[, 1:L_2 + L_1 + 1]
		
		mu_acc <- ise(time_g, mu_g, mu_hat)
		
		Psi_1_acc <- rep(NA, L_1)
		for(l in 1:L_1) {
			
			Psi_1_acc[l] <- ise(time_g, Psi_1_g[, l], Psi_1_hat[, l])
		}
		
		Psi_2_acc <- rep(NA, L_2)
		for(l in 1:L_2) {
			
			Psi_2_acc[l] <- ise(time_g, Psi_2_g[, l], Psi_2_hat[, l])
		}
		
		log_acc_vec <- log(c(mu_acc, Psi_1_acc, Psi_2_acc))
		mlfpca_acc <- c(N, sim, log_acc_vec)
		write(mlfpca_acc, "res/tmp_mlfpca.txt", ncol = n_col, append = TRUE)
		
		# Run vmp_gauss_mlfpca algorithm:
		
		set.seed(1)
		
		eta_vec <- vmp_gauss_mlfpca(
			n_vmp, N, M, L_1, L_2, C, Y,
			sigsq_beta, A, criterion, plot_elbo = FALSE
		)
		
		eta_in <- list(
			eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(zeta)->zeta",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		)
		
		gauss_mlfpca_res <- mlfpc_rotation(
			eta_in, time_g, C_g, L_1, L_2,
			N_sample, M_sample, Psi_g = list(Psi_1_g, Psi_2_g)
		)
		
		gbl_hat <- gauss_mlfpca_res$"gbl_curves"
		mu_hat <- gbl_hat[, 1]
		Psi_1_hat <- gbl_hat[, 1:L_1 + 1]
		Psi_2_hat <- gbl_hat[, 1:L_2 + L_1 + 1]
		
		mu_acc <- ise(time_g, mu_g, mu_hat)
		
		Psi_1_acc <- rep(NA, L_1)
		for(l in 1:L_1) {
			
			Psi_1_acc[l] <- ise(time_g, Psi_1_g[, l], Psi_1_hat[, l])
		}
		
		Psi_2_acc <- rep(NA, L_2)
		for(l in 1:L_2) {
			
			Psi_2_acc[l] <- ise(time_g, Psi_2_g[, l], Psi_2_hat[, l])
		}
		
		log_acc_vec <- log(c(mu_acc, Psi_1_acc, Psi_2_acc))
		gauss_mlfpca_acc <- c(N, sim, log_acc_vec)
		write(gauss_mlfpca_acc, "res/tmp_gauss_mlfpca.txt", ncol = n_col, append = TRUE)
	}
}













     














