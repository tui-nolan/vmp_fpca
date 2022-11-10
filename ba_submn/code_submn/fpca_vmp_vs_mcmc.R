######### R script: fpca_vmp_vs_mcmc.R ##########

# For performing a simple functional principal
# components analysis via MCMC to determine
# the FPCA basis functions. The FPCA basis functions
# are determined via semiparametric regression.

# Created: 29 SEP 2020
# Last changed: 29 JUL 2022

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

set.seed(36)

# Establish simulation variables:

N <- 50                             # number of curves
n_sample <- 4                       # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots
T_vec <- sample(20:30, N, replace = TRUE)    # number of time observations for each curve
n_int_knots <- 10                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 4                              # guess of the number of FPCA basis functions
data_col <- "black"                 # colour of the data in the plots
criterion <- 1e-5                   # convergence criterion
d <- (K+2)*(L+1)              # dimension of spline vector

n_vmp <- 500                        # number of VMP iterations
n_g <- 200                         # length of the plotting grid
vmp_col <- "red"                    # colour of the VMP plots
vmp_lwd <- 1                        # line width for vmp plots

n_burnin <- 1000                    # Length of burn-in.
n_mcmc <- 1000                      # Size of the kept sample.
n_thin <- 1                         # Thinning factor. 
tolerance <- 1e-10
mcmc_col <- "deepskyblue2"          # colour of the MCMC lines in the plots
mcmc_lwd <- 2                       # line width for mcmc plots

sigma_zeta_vec <- 1/(1:L)      # sd for the scores
sigma_eps <- 1                      # sd of the residuals
sigsq_eps <- sigma_eps^2

# Set up plotting variables:

plot_dim <- c(2, 2)                 # (ncol, nrow) for curve plots

n_g <- 1000                         # length of plotting grid

plot_width <- 2
plot_height <- 4

print_pdf <- TRUE

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

# Generate FPCA data:

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

# Plot the data:

plot_fpca_data(N_sample, time_obs, Y, plot_dim, data_col)

####################################################
#
#  MCMC  SIMULATIONS
#
####################################################

# Set up Stan inputs:

all_data <- list(
	N=N, n_time_obs=sum(T_vec), K=K, L=L,
	sigma_beta=sqrt(sigsq_beta), A=A,
	Sigma_zeta=Sigma_zeta,
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

mcmc_summary <- summarise_mcmc(stan_obj, N_sample, C_g, Psi_g)

Y_mcmc_summary <- mcmc_summary$"Y_summary"
gbl_mcmc_hat <- mcmc_summary$"gbl_curves"
Zeta_mcmc_hat <- mcmc_summary$"zeta"

####################################################
#
#  VMP  SIMULATIONS
#
####################################################

# VMP simulations:

eta_vec <- vmp_gauss_fpca(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, plot_elbo = TRUE
)

# Get the posterior estimates

eta_in <- list(
	eta_vec$"p(nu|Sigma_nu)->nu", eta_vec $"p(Y|nu,zeta,sigsq_eps)->nu",
	eta_vec $"p(zeta)->zeta", eta_vec $"p(Y|nu,zeta,sigsq_eps)->zeta"
)
fpc_rotns <- fpc_orthogonalization(eta_in, N_sample, time_g, C_g, Psi_g)

# Summarise the VMP results:

Y_vmp_summary <- fpc_rotns$"Y_summary"
gbl_vmp_hat <- fpc_rotns$"gbl_curves"
Zeta_vmp_hat <- fpc_rotns$"zeta"

####################################################
#
#  PLOT  OF  COMPARISONS
#
####################################################

# Plot the fitted curves:

if(print_pdf) {
	
	pdf("./res/fpca_fits.pdf",width=plot_width, height=plot_height)
}

plot_fpca_fit_comparisons(
	N_sample, time_obs, time_g,
	Y, Y_vmp_summary, Y_mcmc_summary,
	plot_dim, vmp_col, mcmc_col, data_col
)

if(print_pdf) {
	
	dev.off()
}

############ End of fpca_vmp_vs_mcmc.R ############