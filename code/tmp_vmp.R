######### R script: mlfpca_sim_st.R ##########

# For comparing a MlFPCA via MFVB and MCMC

# Created: 19 JUL 2022
# Last changed: 23 JUL 2022

# Load libraries:

library(MASS)
library(magic)
library(rstan)
rstan_options(auto_write = TRUE)
library(lattice)
library(pracma)

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

# Establish simulation variables:

N_vec <- c(10, 50, 100)                    # number of subjects
M_min <- 10                                # minimum second level observations
M_max <- 15                                # maximum second level observations

n_sims <- 100                               # number of simulations

criterion  <- 1e-5                        # convergence criterion
          
n_int_knots <- 10                         # number of interior knots
K <- n_int_knots + 2                      # number of spline basis functions
L_1 <- 3                                  # number of first level eigenfunctions
L_2 <- 3                                  # number of second level eigenfunctions
L <- L_1 + L_2
data_col <- "grey51"                      # colour of the data in the plots

n_burnin <- 500                   # Length of burn-in.
n_mcmc <- 500                     # Size of the kept sample.
n_thin <- 1                        # Thinning factor. 
tolerance <- 1e-10
mcmc_col <- "red"                  # colour of the MCMC lines in the plots
n_g <- 1000                               # length of the plotting grid

n_vmp <- 500                              # number of VMP iterations
n_g <- 1000                               # length of the plotting grid
vmp_col <- "deepskyblue2"                 # colour of the VMP plots
d <- (K+2)*(L_1 + L_2 + 1)                # dimension of spline vector

sigma_zeta_1 <- 1/(1:L_1)                   # vector of st. dev.'s for the first level scores
sigma_zeta_2 <- 1/(1:L_2)                 # vector of st. dev.'s for the second level scores

sigma_eps <- 1                            # sd of the residuals
sigsq_eps <- sigma_eps^2

plot_width <- 6
plot_height <- 4
bw_plot_dim <- c(L + 1, 1)              # c(ncol, nrow)
hist_plot_dim <- c(length(N_vec), 1)    # c(ncol, nrow)

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
mu_beta <- rep(0, 2)
A <- 1e5
sigsq_zeta <- 1
sigma_zeta <- sqrt(sigsq_zeta)
Sigma_zeta <- sigsq_zeta*diag(L)

# Set the mean function and the FPCA basis functions:

mu <- function(t) return(10*sin(pi*t) - 5)
Psi_1 <- fourier_basis(L_1)
Psi_2 <- fourier_basis(L_1 + L_2)[(L_1 + 1):(L_1 + L_2)]

# Begin the VMP simulations:

psi_1_names <- rep(NA, L_1)
for(l in 1:L_1) {
	
	psi_1_names[l] <- paste("psi_1", l, "_vmp", sep="")
}

psi_2_names <- rep(NA, L_2)
for(l in 1:L_2) {
	
	psi_2_names[l] <- paste("psi_2", l, "_vmp", sep="")
}

col_names <- c(
	"N", "sim", "mu_vmp",
	psi_1_names, psi_2_names,
	"zeta_1_vmp",
	"zeta_2_vmp"
)
n_col_acc <- length(col_names)
write(col_names, "res/tmp_mlfpca_acc.txt", ncol = n_col_acc, append = FALSE)

col_names <- c("N", "sim", "iter", "vmp")
n_col_speed <- length(col_names)
write(col_names, "res/tmp_mlfpca_speed.txt", ncol = n_col_speed, append = FALSE)

for(i_N in 1:length(N_vec)) {
	
	N <- N_vec[i_N]
	N_sample <- 1:N
	
	cat("Starting simulations with", N, "response curves \n")
	
	for(i_sim in 1:n_sims) {
		
		cat("Simulation", i_sim, "\n")
		
		set.seed(i_sim + 80)
		
		M <- sample(M_min:M_max, N, replace = TRUE)
		M_sample <- 1:M_min
		
		T_vec <- vector("list", length = N)       # number of time observations for each curve
		for(i in 1:N) {
			
			T_vec[[i]] <- round(runif(M[i], 20, 30))
		}
		
		ml_fpca_data <- gauss_mlfpca_data(
			T_vec, K, n_g, sigma_zeta_1, sigma_zeta_2,
			sigma_eps, mu, Psi_1, Psi_2
		)
		
		time_obs <- ml_fpca_data$time_obs
		time_g <- ml_fpca_data$time_g
		int_knots <- ml_fpca_data$int_knots
		X <- ml_fpca_data$X
		Z <- ml_fpca_data$Z
		C <- ml_fpca_data$C
		X_g <- ml_fpca_data$X_g
		Z_g <- ml_fpca_data$Z_g
		C_g <- ml_fpca_data$C_g
		zeta_1 <- ml_fpca_data$zeta_1
		zeta_2 <- ml_fpca_data$zeta_2
		mu_g <- ml_fpca_data$mu_g
		Psi_1_g <- ml_fpca_data$Psi_1_g
		Psi_2_g <- ml_fpca_data$Psi_2_g
		Y <- ml_fpca_data$Y
		
		####################################################
		#
		#  VMP  SIMULATIONS
		#
		####################################################
		
		start_time <- Sys.time()
		
		vmp_res <- vmp_mlfpca(
			n_vmp, N, M, L_1, L_2, C, Y,
			sigsq_beta, A, criterion, plot_elbo = FALSE
		)
		
		eta_vec <- vmp_res$"eta_vec"
		iter <- vmp_res$"iter"
		
		# Get the posterior estimates
		
		eta_in <- list(
			eta_vec$"p(nu|Sigma_nu)->nu",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec$"p(zeta)->zeta",
			eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
		)
		
		mlfpc_rotns <- mlfpc_rotation(
			eta_in, time_g, C_g, L_1, L_2,
			N_sample, M_sample, Psi_g = list(Psi_1_g, Psi_2_g)
		)
		
		end_time <- Sys.time()
		
		vmp_time <- difftime(end_time, start_time, units="secs")
		
		Zeta_1_hat_vmp <- mlfpc_rotns$"zeta_1"
		Zeta_2_hat_vmp <- mlfpc_rotns$"zeta_2"
		gbl_hat_vmp <- mlfpc_rotns$"gbl_curves"
		
		# Record accuracies:
		
		mu_hat_vmp <- gbl_hat_vmp[, 1]
		Psi_1_hat_vmp <- gbl_hat_vmp[, 1:L_1 + 1]
		Psi_2_hat_vmp <- gbl_hat_vmp[, 1:L_2 + L_1 + 1]
		
		mu_vmp_acc <- ise(time_g, mu_g, mu_hat_vmp)
		
		psi_1_vmp_acc <- rep(NA, L_1)
		for(l in 1:L_1) {
			
			psi_1_vmp_acc[l] <- ise(time_g, Psi_1_g[, l], Psi_1_hat_vmp[, l])
		}
		
		psi_2_vmp_acc <- rep(NA, L_2)
		for(l in 1:L_2) {
			
			psi_2_vmp_acc[l] <- ise(time_g, Psi_2_g[, l], Psi_2_hat_vmp[, l])
		}
		
		norm_diff <- apply(Zeta_1_hat_vmp - Reduce(rbind, zeta_1), 1, function(x) sqrt(cprod(x)))
		rmse_1_vmp <- sqrt(mean(norm_diff))
		
		norm_diff <- apply(
			Reduce(rbind, Zeta_2_hat_vmp) - Reduce(rbind, lapply(zeta_2, Reduce, f = rbind)),
			1, function(x) sqrt(cprod(x))
		)
		rmse_2_vmp <- sqrt(mean(norm_diff))
		
		####################################################
		#
		#  RESULTS
		#
		####################################################
		
		psi_1_acc <- psi_1_vmp_acc
		psi_2_acc <- psi_2_vmp_acc
		
		acc_res <- c(
			N, i_sim, mu_vmp_acc,
			psi_1_acc, psi_2_acc,
			rmse_1_vmp,
			rmse_2_vmp
		)
		write(acc_res, "res/tmp_mlfpca_acc.txt", ncol = n_col_acc, append = TRUE)
		
		speed_results <- c(N, i_sim, iter, vmp_time)
		write(speed_results, "res/tmp_mlfpca_speed.txt", ncol = n_col_speed, append=TRUE)
	}
}

############ End of mlfpca_sim_st.R ############