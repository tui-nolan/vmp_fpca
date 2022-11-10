######### R script: tmp_hpc.R ##########

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

# for HPC

task_id_string <- Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id <- as.numeric(task_id_string)
#task_id <- 1

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

n_sims <- 20

sim <- (task_id - 1) %% n_sims + 1
set.seed(1)

N_index <- (task_id - 1) %/% n_sims + 1
N <- c(10, 50, 100, 250, 500)[N_index]       # number of curves

M_min <- 10
M_max <- 15

m_sample <- M_min                                   # number of second level curves for the plots
n_sample <- N                             # number of curves for the plots

M_sample <- 1:M_min   # specific second level curves for the plots
M <- sample(M_min:M_max, N, replace = TRUE)     # number of observations on each subject

N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots

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

n_vmp <- 500                              # number of VMP iterations
n_g <- 1000                               # length of the plotting grid
vmp_col <- "red"                 # colour of the VMP plots
d <- (K+2)*(L_1 + L_2 + 1)                # dimension of spline vector

n_burnin <- 1000                   # Length of burn-in.
n_mcmc <- 1000                     # Size of the kept sample.
n_thin <- 1                        # Thinning factor. 
tolerance <- 1e-10

sigma_zeta_1 <- 1/(1:L_1)                   # vector of st. dev.'s for the first level scores
sigma_zeta_2 <- 1/(1:L_2)                 # vector of st. dev.'s for the second level scores
l_zeta <- L_1 + M*L_2                     # length of score vector

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

# Generate the data:

set.seed(sim)

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

vmp_res <- mlfpc_rotation(
	eta_in, time_g, C_g, L_1, L_2,
	N_sample, M_sample, Psi_g = list(Psi_1_g, Psi_2_g)
)

end_time <- Sys.time()

vmp_time <- difftime(end_time, start_time, units="secs")

gbl_hat_vmp <- vmp_res$"gbl_curves"
Zeta_1_hat_vmp <- vmp_res$"zeta_1"
Zeta_2_hat_vmp <- vmp_res$"zeta_2"

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
	"Zeta_1_hat_vmp",
	"Zeta_2_hat_vmp"
)
n_col_acc <- length(col_names)

file_name <- paste0("hpc_files/mlfpca_acc_res/tmp_", task_id, ".txt")
write(col_names, file_name, ncol = n_col_acc, append = FALSE)

acc_res <- c(
	N, task_id, mu_vmp_acc,
	psi_1_vmp_acc, psi_2_vmp_acc,
	rmse_1_vmp,
	rmse_2_vmp
)
write(acc_res, file_name, ncol = n_col_acc, append = TRUE)

