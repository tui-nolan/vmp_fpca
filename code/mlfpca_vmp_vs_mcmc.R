######### R script: mlfpca_vmp_vs_mcmc.R ##########

# For performing a comparison of MlFPCA via MCMC
# and VMP

# Created: 12 MAY 2022
# Last changed: 12 MAY 2022

# Updates:
# 1. Constructing comparisons between VMP and MCMC.

# Load libraries:

library(MASS)
library(magic)
library(lattice)
library(pracma)
library(rstan)
rstan_options(auto_write = TRUE)
library(ellipse)
library(matrixcalc)

# Required functions:

setwd("functions")

source("X_design.r")
source("ZOSull.r")
source("OmegaOSull.r")
source("vec.r")
source("vecInverse.r")
source("tr.r")
source("cprod.r")
source("wait.r")
source("trapint.r")
source("stan_mods.r")
source("vmp_functions.r")
source("fourier_basis.r")
source("fpca_algs.r")

setwd("..")

#set.seed(1)

# Establish simulation variables:

M_min <- 10
M_max <- 15
N <- 50                                   # number of subjects

m_sample <- 3                                   # number of second level curves for the plots
n_sample <- 2                             # number of curves for the plots

M_sample <- round(seq(1, M_max, length.out = m_sample))   # specific second level curves for the plots
M <- sample(M_min:M_max, N, replace = TRUE)     # number of observations on each subject
while(sum(M == M_max) < n_sample) {
	
	M <- sample(M_min:M_max, N, replace = TRUE)
}

N_sample <- sort(sample(which(M == M_max), n_sample))   # specific curves for the plots

T_vec <- vector("list", length = N)       # number of time observations for each curve
for(i in 1:N) {
	
	T_vec[[i]] <- round(runif(M[i], 20, 30))
}

criterion  <- 1e-5                        # convergence criterion
          
n_int_knots <- 10                         # number of interior knots
K <- n_int_knots + 2                      # number of spline basis functions
L_1 <- 3                                  # number of first level eigenfunctions
L_2 <- 3                                  # number of second level eigenfunctions
L <- L_1 + L_2
data_col <- "grey51"                      # colour of the data in the plots

n_vmp <- 500                              # number of VMP iterations
n_g <- 1000                               # length of the plotting grid
vmp_col <- "red"                 # colour of the VMP plots
d <- (K+2)*(L_1 + L_2 + 1)                # dimension of spline vector

n_burnin <- 1000                   # Length of burn-in.
n_mcmc <- 1000                     # Size of the kept sample.
n_thin <- 1                        # Thinning factor. 
tolerance <- 1e-10
mcmc_col <- "deepskyblue2"                  # colour of the MCMC lines in the plots
mcmc_lwd <- 2                       # line width for mcmc plots

sigma_zeta_1 <- 1/(1:L_1)                   # vector of st. dev.'s for the first level scores
sigma_zeta_2 <- 1/(1:L_2)                 # vector of st. dev.'s for the second level scores
l_zeta <- L_1 + M*L_2                     # length of score vector

sigma_eps <- 1                            # sd of the residuals
sigsq_eps <- sigma_eps^2

# Set up plot-grid dimensions:

plot_dim <- c(m_sample, n_sample)         # (ncol, nrow) for curve plots

plot_width <- 5
plot_height <- 4

print_pdf <- FALSE              # save the plots in a PDF?

# Set the FPCA basis functions:

mu <- function(t) return(3*sin(pi*t) - 3/2)
Psi_1 <- fourier_basis(L_1)
Psi_2 <- fourier_basis(L_1 + L_2)[(L_1 + 1):(L_1 + L_2)]

# Generate the data:

#set.seed(1)

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

# Plot the data:

plot_mlfpca_data(N_sample, M_sample, time_obs, Y, plot_dim, data_col)

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
mu_beta <- rep(0, 2)
A <- 1e5
sigma_zeta <- 1
sigsq_zeta <- sigma_zeta^2

####################################################
#
#  VMP  SIMULATIONS
#
####################################################

# VMP simulations:

set.seed(1)

vmp_alg <- vmp_mlfpca(
	n_vmp, N, M, L_1, L_2, C, Y,
	sigsq_beta, A, criterion, plot_elbo = FALSE
)

eta_vec <- vmp_alg$"eta_vec"

# Orthogonalization and Rotation:

eta_in <- list(
	eta_vec$"p(nu|Sigma_nu)->nu",
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
	eta_vec$"p(zeta)->zeta",
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
)

vmp_res <- mlfpc_rotation(
	eta_in, time_g, C_g, L_1, L_2,
	N_sample, M_sample, Psi_g = list(Psi_1_g, Psi_2_g)
)

# Summarise the VMP results:

Y_vmp_summary <- vmp_res$"fits"
gbl_hat_vmp <- vmp_res$"gbl_curves"
Zeta_1_vmp <- vmp_res$"zeta_1"
Zeta_2_vmp <- vmp_res$"zeta_2"

####################################################
#
#  MCMC  SIMULATIONS
#
####################################################

# Set up Stan inputs:

n_time_obs <- sum(sapply(T_vec, sum))
X_in <- Reduce(rbind, lapply(X, function(X) Reduce(rbind, X)))
Z_in <- Reduce(rbind, lapply(Z, function(X) Reduce(rbind, X)))
T_vec <- Reduce(c, T_vec)
Y_in <- Reduce(c, lapply(Y, function(X) Reduce(c, X)))

all_data <- list(
	N = N, M = M, n_time_obs = n_time_obs, K = K, L_1 = L_1, L_2 = L_2,
	sigma_beta = sqrt(sigsq_beta), A = A,
	X = X_in, Z = Z_in, T_vec = T_vec, Y = Y_in
)

# Compile Stan code:

compile_obj <- stan(model_code = mlfpca_model, data = all_data, iter = 1, chains = 1)

# Obtain MCMC samples for each parameter using Stan:

stan_obj <- stan(
	model_code = mlfpca_model, data = all_data, warmup = n_burnin,
	iter = (n_burnin + n_mcmc), chains = 1, thin = n_thin,
	refresh = 100, fit = compile_obj
)

# Summarise the MCMC samples:

mcmc_res <- summarise_ml_mcmc(stan_obj, N, M, C_g, Psi_1_g, Psi_2_g, N_sample, M_sample)

Y_mcmc_summary <- mcmc_res$"fits"
gbl_hat_mcmc <- mcmc_res$"gbl_curves"
Zeta_1_mcmc <- mcmc_res$"zeta_1"
Zeta_2_mcmc <- mcmc_res$"zeta_2"

####################################################
#
#  PLOT  OF  COMPARISONS
#
####################################################

# Plot the fitted curves:

if(print_pdf) {
	
	pdf("./res/mlfpca_fits.pdf",width=plot_width, height=plot_height)
}

plot_mlfpca_fit_comparisons(
	N_sample, M_sample, time_obs, Y,
	time_g, Y_vmp_summary, Y_mcmc_summary,
	plot_dim, data_col, vmp_col, mcmc_col
)

if(print_pdf) {
	
	dev.off()
}

############ End of mlfpca_vmp_vs_mcmc_final.R ############