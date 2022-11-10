######### R script: mlfpca_hpc.R ##########

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

n_sims <- 100

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

mu <- function(t) return(3*sin(pi*t) - 3/2)
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
#  MCMC  SIMULATIONS
#
####################################################

start_time <- Sys.time()

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

end_time <- Sys.time()

mcmc_time <- difftime(end_time, start_time, units="secs")

gbl_hat_mcmc <- mcmc_res$"gbl_curves"
Zeta_1_hat_mcmc <- mcmc_res$"zeta_1"
Zeta_2_hat_mcmc <- mcmc_res$"zeta_2"

# Record accuracies:

mu_hat_mcmc <- gbl_hat_mcmc[, 1]
Psi_1_hat_mcmc <- gbl_hat_mcmc[, 1:L_1 + 1]
Psi_2_hat_mcmc <- gbl_hat_mcmc[, 1:L_2 + L_1 + 1]

mu_mcmc_acc <- ise(time_g, mu_g, mu_hat_mcmc)

psi_1_mcmc_acc <- rep(NA, L_1)
for(l in 1:L_1) {
	
	psi_1_mcmc_acc[l] <- ise(time_g, Psi_1_g[, l], Psi_1_hat_mcmc[, l])
}

psi_2_mcmc_acc <- rep(NA, L_2)
for(l in 1:L_2) {
	
	psi_2_mcmc_acc[l] <- ise(time_g, Psi_2_g[, l], Psi_2_hat_mcmc[, l])
}

norm_diff <- apply(Zeta_1_hat_mcmc - Reduce(rbind, zeta_1), 1, function(x) sqrt(cprod(x)))
rmse_1_mcmc <- sqrt(mean(norm_diff))

norm_diff <- apply(
	Reduce(rbind, Zeta_2_hat_mcmc) - Reduce(rbind, lapply(zeta_2, Reduce, f = rbind)),
	1, function(x) sqrt(cprod(x))
)
rmse_2_mcmc <- sqrt(mean(norm_diff))

####################################################
#
#  RESULTS
#
####################################################

psi_1_names <- matrix(NA, 2, L_1)
for(l in 1:L_1) {
	
	psi_1_names[1, l] <- paste("psi_1", l, "_vmp", sep="")
	psi_1_names[2, l] <- paste("psi_1", l, "_mcmc", sep="")
}
psi_1_names <- as.vector(psi_1_names)

psi_2_names <- matrix(NA, 2, L_2)
for(l in 1:L_2) {
	
	psi_2_names[1, l] <- paste("psi_2", l, "_vmp", sep="")
	psi_2_names[2, l] <- paste("psi_2", l, "_mcmc", sep="")
}
psi_2_names <- as.vector(psi_2_names)

col_names <- c(
	"N", "sim", "mu_vmp", "mu_mcmc",
	psi_1_names, psi_2_names,
	"zeta_1_vmp", "zeta_1_mcmc",
	"zeta_2_vmp", "zeta_2_mcmc"
)
n_col_acc <- length(col_names)

file_name <- paste0("hpc_files/mlfpca_acc_res/mlfpca_acc_", task_id, ".txt")
write(col_names, file_name, ncol = n_col_acc, append = FALSE)

psi_1_acc <- as.vector(rbind(psi_1_vmp_acc, psi_1_mcmc_acc))
psi_2_acc <- as.vector(rbind(psi_2_vmp_acc, psi_2_mcmc_acc))

acc_res <- c(
	N, sim, mu_vmp_acc, mu_mcmc_acc,
	psi_1_acc, psi_2_acc,
	rmse_1_vmp, rmse_1_mcmc,
	rmse_2_vmp, rmse_2_mcmc
)
write(acc_res, file_name, ncol = n_col_acc, append = TRUE)

col_names <- c("N", "sim", "iter", "vmp", "mcmc")
n_col_speed <- length(col_names)
file_name <- paste0("hpc_files/mlfpca_speed_res/mlfpca_speed_", task_id, ".txt")
write(col_names, file_name, ncol = n_col_speed, append = FALSE)
speed_res <- c(N, sim, iter, vmp_time, mcmc_time)
write(speed_res, file_name, ncol = n_col_speed, append=TRUE)

############ End of mlfpca_hpc.R ############