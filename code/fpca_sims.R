######### R script: fpca_sims.R ##########

# For running a speed comparison between the variational
# Bayesian approach and the MCMC approach.

# Created: 16 OCT 2020
# Last changed: 23 FEB 2020

# Load libraries:

library(MASS)
library(magic)
library(rstan)
rstan_options(auto_write = TRUE)
library(lattice)
library(ellipse)
library(matrixcalc)
library(pracma)

set.seed(0)

# Required functions:

source("X_design.r")
source("ZOSull.r")
source("OmegaOSull.r")
source("vec.r")
source("vecInverse.r")
source("tr.r")
source("trapint.r")
source("cprod.r")
source("wait.r")
source("ise.r")
source("logistic.r")
source("vmp_functions.r")
source("fpca_algs.r")

N_vec <- c(10, 50, 100)   # number of curves
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 3                              # number of FPCA basis functions
criterion <- 1e-5                   # convergence criterion
d <- (K+2)*(L+1)                    # dimension of spline vector

n_vmp <- 500                        # number of VMP iterations
n_mc <- 100                         # number of MC samples for MFVB CI
n_g <- 1000                         # length of the plotting grid

n_burnin <- 1000                    # Length of burn-in.
n_mcmc <- 1000                      # Size of the kept sample.
n_thin <- 1                         # Thinning factor. 
tolerance <- 1e-10

sigma_zeta_vec <- c(1, 0.5, 0)      # sd for the scores
sigma_eps <- 1                      # sd of the residuals
sigsq_eps <- sigma_eps^2

n_sims <- 100                       # number of simulations

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
psi_1 <- function(t) return(sqrt(2)*sin(2*pi*t))
psi_2 <- function(t) return(sqrt(2)*cos(2*pi*t))
psi_3 <- function(t) return(0)
Psi_func <- list(psi_1, psi_2, psi_3)

L_true <- 2                            # true number of basis functions

# Set up the FPCA model:

fpca_model <- "
	
	data {
		
		int<lower=1> N;                // number of curves
		int<lower=N> n_time_obs;       // total number of time observations
		int<lower=1> K;                // number of splines
		int<lower=1> L;                // number of basis functions
		real<lower=0> sigma_beta;      // fixed effects prior variance
		real<lower=0> A;               // cauchy hyperparameter
		real<lower=0> sigma_zeta;      // prior variance of the scores
		matrix[n_time_obs, 2] X;       // rbind of all design matrices
		matrix[n_time_obs, K] Z;       // rbind of all spline design matrices
		int<lower=1> T_vec[N];         // vector of time observations for
		                               // each curve
		vector[n_time_obs] Y;          // vector of all responses
	}
	
	parameters {
		
		matrix[N,L] zeta;
		
		real<lower=0> sigma_eps;
		
		vector[2] beta_mu;
		vector[K] u_mu;
		real<lower=0> sigma_mu;
		
		matrix[L, 2] beta_psi;
		matrix[L, K] u_psi;
		vector<lower=0>[L] sigma_psi;
	}
	
	transformed parameters {
		
		vector[n_time_obs] mu;
		matrix[L, n_time_obs] psi;
		
		mu = X*beta_mu + Z*u_mu;
		
		for(l in 1:L) {
			
			psi[l] = beta_psi[l]*X' + u_psi[l]*Z';
		}
	}
	
	model {
		
		int pos;
		pos = 1;
		
		for(i in 1:N) {
			
			vector[T_vec[i]] mu_i;
			matrix[L, T_vec[i]] psi_i;
			vector[T_vec[i]] Y_i_hat;
			
			mu_i = segment(mu, pos, T_vec[i]);
			psi_i = block(psi, 1, pos, L, T_vec[i]);
			Y_i_hat = mu_i + to_vector(zeta[i]*psi_i);
			
			segment(Y, pos, T_vec[i]) ~ normal(Y_i_hat, sigma_eps);
			
			pos = pos + T_vec[i];
			
			zeta[i] ~ normal(0, sigma_zeta);
		}
		
		sigma_eps ~ cauchy(0, A);
		
		beta_mu ~ normal(0, sigma_beta);
		u_mu ~ normal(0, sigma_mu);
		sigma_mu ~ cauchy(0, A);
		
		for(l in 1:L) {
			
			beta_psi[l] ~ normal(0, sigma_beta);
			u_psi[l] ~ normal(0, sigma_psi[l]);
			sigma_psi[l] ~ cauchy(0, A);
		}
	}
"

# Begin the speed simulations:

col_names <- c("N", "sim", "VMP", "MCMC")
n_col <- length(col_names)
write(col_names, "./res/comp_speed_res.txt", ncol=n_col, append=FALSE)

psi_names <- rep(NA, L_true)
for(l in 1: L_true) {
	
	psi_names[l] <- paste("psi_", l, sep="")
}
col_names <- c("N", "sim", "mu", psi_names)
n_col <- length(col_names)
write(col_names, "./res/gauss_fpca_acc.txt", ncol=n_col, append=FALSE)

write(seq(0, 1, length.out=n_g), "./res/mu.txt", ncol=n_g, append=FALSE)

for(l in 1:L) {
	
	txt_name <- paste("./res/psi_", l, ".txt", sep="")
	
	write(seq(0, 1, length.out=n_g), txt_name, ncol=n_g, append=FALSE)
}

write(seq(0, 1, length.out=n_g), "./res/mu_mcmc.txt", ncol=n_g, append=FALSE)

for(l in 1:L) {
	
	txt_name <- paste("./res/psi_mcmc_", l, ".txt", sep="")
	
	write(seq(0, 1, length.out=n_g), txt_name, ncol=n_g, append=FALSE)
}

for(i_N in 1:length(N_vec)) {
	
	N <- N_vec[i_N]
	
	cat("Starting simulations with", N, "response curves \n")
	
	for(i_sim in 1:n_sims) {
		
		set.seed(i_sim)
		
		T_vec <- round(runif(N, 20, 30))
		
		cat("starting simulation", i_sim, "of", n_sims, "\n")
		
		# Gather the data:
		
		fpca_data <- gauss_fpca_data(
			T_vec, N, K, n_g, sigma_zeta_vec, sigma_eps,
			mu, Psi_func
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
		
		# Run MCMC simulations:
		
		start_time <- Sys.time()
		
		all_data <- list(
			N=N, n_time_obs=sum(T_vec), K=K, L=L,
			sigma_beta=sqrt(sigsq_beta), A=A,
			sigma_zeta=sigma_zeta,
			X=do.call(rbind, X),
			Z=do.call(rbind, Z),
			T_vec=T_vec, Y=Y_vec
		)
		
		compile_obj <- stan(
			model_code=fpca_model, data=all_data,
			iter=1, chains=1
		)
		
		stan_obj <- stan(
			model_code=fpca_model, data=all_data, warmup=n_burnin,
			iter=(n_burnin+n_mcmc), chains=1, thin=n_thin,
			refresh=100, fit=compile_obj
		)
		
		mcmc_summary <- summarise_mcmc(stan_obj, C_g, Psi_g, use_logistic_mod=FALSE)
		
		end_time <- Sys.time()
		
		mcmc_time <- difftime(end_time, start_time, units="secs")
		
		# Run VMP simulations:
		
		start_time <- Sys.time()
		
		eta_vec <- vmp_gauss_fpca(
			n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
			Sigma_beta, A, time_g, C_g, Psi_g,
			criterion, n_mc=100, plot_elbo=FALSE
		)
		
		# Get the posterior estimates
		
		eta_in <- list(
			eta_vec$"p(nu|Sigma_nu)->nu", eta_vec $"p(Y|nu,zeta,sigsq_eps)->nu",
			eta_vec $"p(zeta)->zeta", eta_vec $"p(Y|nu,zeta,sigsq_eps)->zeta"
		)
		fpc_rotns <- fpc_rotation(eta_in, time_g, C_g, Psi_g)
		
		end_time <- Sys.time()
		
		vmp_time <- difftime(end_time, start_time, units="secs")
		
		gbl_estimates <- fpc_rotns$"gbl_summary"[,1:(L_true+1)]
		
		mu_q_mu <- gbl_estimates[,1]
		M_q_Psi <- gbl_estimates[,2:(L_true+1)]
		
		# Save the results:
		
		speed_results <- c(N, i_sim, vmp_time, mcmc_time)
		write(speed_results, "./res/comp_speed_res.txt", ncol=n_col, append=TRUE)
		
		mu_acc <- ise(time_g, mu_g, mu_q_mu)
		psi_acc <- rep(NA, L_true)
		for(l in 1:L_true) {
			
			psi_acc[l] <- ise(time_g, Psi_g[,l], M_q_Psi[,l])
		}
		
		results <- c(N, i_sim, mu_acc, psi_acc)
		write(results, "./res/gauss_fpca_acc.txt", ncol=n_col, append=TRUE)
		
		# Save the curves for the final value of N:
		
		if(i_N==length(N_vec)) {
			
			write(mu_q_mu, "./res/mu.txt", ncol=n_g, append=TRUE)
			
			for(l in 1:L_true) {
				
				txt_name <- paste("./res/psi_", l, ".txt", sep="")
				
				write(M_q_Psi[,l], txt_name, ncol=n_g, append=TRUE)
			}
			
			gbl_mcmc_summary <- mcmc_summary$"gbl_mcmc_summary"
			
			mu_g_mcmc <- gbl_mcmc_summary[,1]
			Psi_g_mcmc <- gbl_mcmc_summary[,2:(L_true+1)]
			
			write(mu_g_mcmc, "./res/mu_mcmc.txt", ncol=n_g, append=TRUE)
			
			for(l in 1:L_true) {
				
				txt_name <- paste("./res/psi_mcmc_", l, ".txt", sep="")
				
				write(Psi_g_mcmc[,l], txt_name, ncol=n_g, append=TRUE)
			}
		}
	}
}

############ End of fpca_sims.R ############