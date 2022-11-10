######### R script: fpca_L_sims.R ##########

# For performing Bayesian FPCA simulation study for
# determining the number of eigenfunctions (L)

# Created: 20 JUL 2022
# Last changed: 21 JUL 2022

# Load libraries:

library(MASS)
library(magic)
library(lattice)
library(pracma)
library(ellipse)

set.seed(0)

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

N_vec <- c(10, 50, 100)             # number of curves
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 4                              # number of FPCA basis functions
L_vec <- c(2, 4, 6)              # vector of guesses at L
data_col <- "red"                   # colour of the data in the plots
criterion <- 1e-6                   # convergence criterion

n_vmp <- 500                        # number of VMP iterations
n_g <- 1000                         # length of the plotting grid
vmp_col <- "black"                  # colour of the VMP plots
d <- (K+2)*(L+1)                    # dimension of spline vector

n_burnin <- 500                    # Length of burn-in.
n_mcmc <- 500                      # Size of the kept sample.
n_thin <- 1                        # Thinning factor. 
tolerance <- 1e-10
mcmc_col <- "red"                  # colour of the MCMC lines in the plots

sigma_zeta_vec <- 1/(1:L)         # sd for first and second scores
sigma_eps <- 1                      # sd of the residuals
sigsq_eps <- sigma_eps^2

n_sims <- 10                       # number of simulations

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
mu_beta <- rep(0, 2)
A <- 1e5
sigsq_zeta <- 1
sigma_zeta <- sqrt(sigsq_zeta)
Sigma_zeta <- sigsq_zeta*diag(L)

# Set the mean function and the FPCA basis functions:

mu <- function(t) return(3*sin(pi*t))
Psi <- fourier_basis(L)

# Begin the VMP simulations:

evec_names <- matrix(NA, 2, max(L_vec))
for(k in 1:max(L_vec)) {
	
	evec_names[1, k] <- paste("evec_", k, "_vmp", sep="")
	evec_names[2, k] <- paste("evec_", k, "_mcmc", sep="")
}
evec_names <- as.vector(evec_names)
col_names <- c("N", "sim", "L", evec_names)
n_col <- length(col_names)
write(col_names, "res/fpca_L.txt", ncol = n_col, append = FALSE)

for(i_N in 1:length(N_vec)) {
	
	N <- N_vec[i_N]
	N_sample <- 1:N
	
	cat("Starting simulations with", N, "response curves \n")
	
	for(i_sim in 1:n_sims) {
		
		cat("starting simulation", i_sim, "of", n_sims, "\n")
		
		set.seed(i_sim)
		
		T_vec <- round(runif(N, 20, 30))
		
		# Gather the data:
		
		fpca_data <- gauss_fpca_data(
			T_vec, N, K, n_g, sigma_zeta_vec, sigma_eps,
			mu, Psi
		)
		
		time_obs <- fpca_data$"time_obs"
		time_g <- fpca_data$"time_g"
		int_knots <- fpca_data$"int_knots"
		X <- fpca_data$"X"
		Z <- fpca_data$"Z"
		C <- fpca_data$"C"
		X_g <- fpca_data$"X_g"
		Z_g <- fpca_data$"Z_g"
		C_g <- fpca_data$"C_g"
		zeta <- fpca_data$"zeta"
		mu_t <- fpca_data$"mu_t"
		Psi_t <- fpca_data$"Psi_t"
		mu_g <- fpca_data$"mu_g"
		Psi_g <- fpca_data$"Psi_g"
		Y <- fpca_data$"Y"
		
		Y_vec <- Reduce(c, Y)
		
		for(L_k in L_vec) {
			
			Psi_g_k <- Psi_g[, 1:min(L, L_k)]
			
			if(L_k > L) {
				
				Zero <- matrix(0, n_g, L_k - L)
				Psi_g_k <- cbind(Psi_g, Zero)
			} else {
				
				Psi_g_k <- Psi_g[, 1:L_k]
			}
			
			####################################################
			#
			#  VMP  SIMULATIONS
			#
			####################################################
			
			vmp_res <- vmp_fpca(
				n_vmp, N, L_k, C, Y, sigma_zeta, mu_beta,
				Sigma_beta, A, time_g, C_g, Psi_g_k,
				criterion, plot_elbo = FALSE
			)
			
			eta_vec <- vmp_res$"eta_vec"
			
			# Get the posterior estimates
			
			eta_in <- list(
				eta_vec$"p(nu|Sigma_nu)->nu", eta_vec $"p(Y|nu,zeta,sigsq_eps)->nu",
				eta_vec $"p(zeta)->zeta", eta_vec $"p(Y|nu,zeta,sigsq_eps)->zeta"
			)
			fpc_rotns <- fpc_orthogonalization(eta_in, N_sample, time_g, C_g, Psi_g_k)
			
			Zeta_hat_vmp <- fpc_rotns$"zeta"
			gamma_hat_vmp <- apply(Zeta_hat_vmp, 2, var)
			cum_prop_var <- cumsum(gamma_hat_vmp/sum(gamma_hat_vmp))
			
			vmp_res <- rep(1, max(L_vec))
			vmp_res[1:L_k] <- cum_prop_var
			
			####################################################
			#
			#  MCMC  SIMULATIONS
			#
			####################################################
			
			all_data <- list(
				N=N, n_time_obs=sum(T_vec), K=K, L=L_k,
				sigma_beta=sqrt(sigsq_beta), A=A,
				Sigma_zeta=diag(L_k),
				X=do.call(rbind, X),
				Z=do.call(rbind, Z),
				T_vec=T_vec, Y=Y_vec
			)
			
			# Compile Stan code:
			
			compile_obj <- stan(
				model_code=fpca_model, data=all_data,
				iter=1, chains=1
			)
			
			# Obtain MCMC samples for each parameter using Stan:
			
			stan_obj <- stan(
				model_code=fpca_model, data=all_data, warmup=n_burnin,
				iter=(n_burnin+n_mcmc), chains=1, thin=n_thin,
				refresh=100, fit=compile_obj
			)
			
			# Summarise the MCMC samples:
			
			mcmc_summary <- summarise_mcmc(stan_obj, N_sample, C_g, Psi_g_k)
			
			Zeta_hat_mcmc <- mcmc_summary$"zeta"
			gamma_hat_mcmc <- apply(Zeta_hat_mcmc, 2, var)
			cum_prop_var <- cumsum(gamma_hat_mcmc/sum(gamma_hat_mcmc))
			
			mcmc_res <- rep(1, max(L_vec))
			mcmc_res[1:L_k] <- cum_prop_var
			
			####################################################
			#
			#  RESULTS
			#
			####################################################
			
			evec_res <- as.vector(rbind(vmp_res, mcmc_res))
			L_res <- c(N, i_sim, L_k, evec_res)
			write(L_res, "res/fpca_L.txt", ncol = n_col, append = TRUE)
		}
	}
}

exp_var <- fpca_exp_var_summary("res/fpca_L.txt")

############ End of gauss_fpca_sim_st_4.R ############