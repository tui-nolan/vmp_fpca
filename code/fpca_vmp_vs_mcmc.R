######### R script: fpca_vmp_vs_mcmc.R ##########

# For performing a simple functional principal
# components analysis via MCMC to determine
# the FPCA basis functions. The FPCA basis functions
# are determined via semiparametric regression.

# Created: 29 SEP 2020
# Last changed: 23 APR 2021

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
source("vmp_functions.r")
source("fpca_algs.r")

setwd("..")

# Establish simulation variables:

N <- 50                             # number of curves
n_sample <- 4                       # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots
T_vec <- round(runif(N, 20, 30))    # number of time observations for each curve
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 3                              # guess of the number of FPCA basis functions
data_col <- "black"                 # colour of the data in the plots
criterion <- 1e-5                   # convergence criterion
d <- (K+2)*(L+1)              # dimension of spline vector

n_vmp <- 200                        # number of VMP iterations
n_mc <- 100                         # number of MC samples for MFVB CI
n_g <- 1000                         # length of the plotting grid
vmp_col <- "red"                    # colour of the VMP plots
vmp_lwd <- 1                        # line width for vmp plots

n_burnin <- 1000                    # Length of burn-in.
n_mcmc <- 1000                      # Size of the kept sample.
n_thin <- 1                         # Thinning factor. 
tolerance <- 1e-10
mcmc_col <- "deepskyblue2"          # colour of the MCMC lines in the plots
mcmc_lwd <- 2                       # line width for mcmc plots

sigma_zeta_vec <- c(1, 0.5, 0)      # sd for the scores
sigma_eps <- 0.5                      # sd of the residuals
sigsq_eps <- sigma_eps^2

# Set up plotting variables:

plot_dim <- c(2, 2)                 # (ncol, nrow) for curve plots
plot_gbl_dim <- c(3, 1)             # (ncol, nrow) for the mean and basis function plots

n_g <- 1000                         # length of plotting grid

plot_width <- 2.9
plot_height <- 3.5

print_pdf <- FALSE

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

L_true <- 2                            # true number of basis functions

# Generate FPCA data:

fpca_data <- gauss_fpca_data(
	T_vec, N, K, n_g, sigma_zeta_vec, sigma_eps,
	mu, list(psi_1, psi_2, psi_3)
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

# Plot the data:

if(print_pdf) {
	
	pdf("./res/raw_gauss_data.pdf",width=plot_width, height=plot_height)
}

plot_fpca_data(N_sample, time_obs, Y, plot_dim, data_col)

if(print_pdf) {
	
	dev.off()
} else {
	
	wait()
}

####################################################
#
#  MCMC  SIMULATIONS
#
####################################################

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

# Set up Stan inputs:

all_data <- list(
	N=N, n_time_obs=sum(T_vec), K=K, L=L,
	sigma_beta=sqrt(sigsq_beta), A=A,
	sigma_zeta=sigma_zeta,
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

mcmc_summary <- summarise_mcmc(stan_obj, C_g, Psi_g, use_logistic_mod=FALSE)

Y_mcmc_summary <- mcmc_summary$"Y_g_mcmc_summary"
gbl_mcmc_summary <- mcmc_summary$"gbl_mcmc_summary"
zeta_mcmc_summary <- mcmc_summary$"zeta_mcmc_summary"

####################################################
#
#  VMP  SIMULATIONS
#
####################################################

# VMP simulations:

eta_vec <- vmp_gauss_fpca(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, plot_elbo=TRUE
)

# Get the posterior estimates

eta_in <- list(
	eta_vec$"p(nu|Sigma_nu)->nu", eta_vec $"p(Y|nu,zeta,sigsq_eps)->nu",
	eta_vec $"p(zeta)->zeta", eta_vec $"p(Y|nu,zeta,sigsq_eps)->zeta"
)
fpc_rotns <- fpc_rotation(eta_in, time_g, C_g, Psi_g)

# Summarise the VMP results:

Y_vmp_summary <- fpc_rotns$"Y_summary"
zeta_vmp_summary <- fpc_rotns$"zeta_summary"
gbl_vmp_summary <- fpc_rotns$"gbl_summary"

####################################################
#
#  PLOT  OF  COMPARISONS
#
####################################################

# Plot the fitted curves:

if(print_pdf) {
	
	pdf("./res/response_curves.pdf",width=plot_width, height=plot_height)
}

plot_fit_comparisons(
	N_sample, time_obs, time_g,
	Y, Y_vmp_summary, Y_mcmc_summary,
	plot_dim, vmp_col, mcmc_col, data_col
)

if(print_pdf) {
	
	dev.off()
} else {
	
	wait()
}

# Set up the plot of the scores:

if(print_pdf) {
	
	pdf("./res/scores.pdf",width=plot_width, height=plot_height)
}

plot_score_comparisons(
	N_sample, zeta, zeta_vmp_summary, zeta_mcmc_summary,
	data_col, vmp_col, mcmc_col, plot_dim,
	vmp_lwd = 1, mcmc_lwd = 2
)

if(print_pdf) {
	
	dev.off()
}

############ End of fpca_vmp_vs_mcmc.R ############