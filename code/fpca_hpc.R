######### R script: fpca_vmp_vs_mcmc.R ##########

# For performing a simple functional principal
# components analysis via MCMC to determine
# the FPCA basis functions. The FPCA basis functions
# are determined via semiparametric regression.

# Created: 03 AUG 2020
# Last changed: 16 AUG 2022

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
set.seed(sim)

N_index <- (task_id - 1) %/% n_sims + 1
N <- c(10, 50, 100, 250, 500)[N_index]       # number of curves
n_sample <- N                                # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))      # specific curves for the plots
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

sigma_zeta_vec <- 1/(1:L)      # sd for the scores
sigma_eps <- 1                      # sd of the residuals
sigsq_eps <- sigma_eps^2

n_g <- 1000                         # length of plotting grid

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

####################################################
#
#  VMP  SIMULATIONS
#
####################################################

start_time <- Sys.time()

vmp_res <- vmp_fpca(
	n_vmp, N, L, C, Y, sigma_zeta, mu_beta,
	Sigma_beta, A, time_g, C_g, Psi_g,
	criterion, plot_elbo=FALSE
)

eta_vec <- vmp_res$"eta_vec"
iter <- vmp_res$"iter"

# Get the posterior estimates

eta_in <- list(
	eta_vec$"p(nu|Sigma_nu)->nu", eta_vec $"p(Y|nu,zeta,sigsq_eps)->nu",
	eta_vec $"p(zeta)->zeta", eta_vec $"p(Y|nu,zeta,sigsq_eps)->zeta"
)
fpc_rotns <- fpc_orthogonalization(eta_in, N_sample, time_g, C_g, Psi_g)

end_time <- Sys.time()

vmp_time <- difftime(end_time, start_time, units="secs")

Zeta_hat_vmp <- fpc_rotns$"zeta"
gbl_hat_vmp <- fpc_rotns$"gbl_curves"

# Record accuracies:

mu_hat_vmp <- gbl_hat_vmp[, 1]
Psi_hat_vmp <- gbl_hat_vmp[, 2:(L + 1)]

mu_vmp_acc <- ise(time_g, mu_g, mu_hat_vmp)
psi_vmp_acc <- rep(NA, L)
for(l in 1:L) {
	
	psi_vmp_acc[l] <- ise(time_g, Psi_g[, l], Psi_hat_vmp[,l])
}

norm_diff <- apply(Zeta_hat_vmp - Reduce(rbind, zeta), 1, function(x) sqrt(cprod(x)))
rmse_vmp <- sqrt(mean(norm_diff))

####################################################
#
#  MCMC  SIMULATIONS
#
####################################################

start_time <- Sys.time()

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

end_time <- Sys.time()

mcmc_time <- difftime(end_time, start_time, units="secs")

gbl_hat_mcmc <- mcmc_summary$"gbl_curves"
Zeta_hat_mcmc <- mcmc_summary$"zeta"

# Record accuracies:

mu_hat_mcmc <- gbl_hat_mcmc[, 1]
Psi_hat_mcmc <- gbl_hat_mcmc[, 2:(L + 1)]

mu_mcmc_acc <- ise(time_g, mu_g, mu_hat_mcmc)
psi_mcmc_acc <- rep(NA, L)
for(l in 1:L) {
	
	psi_mcmc_acc[l] <- ise(time_g, Psi_g[, l], Psi_hat_mcmc[,l])
}

norm_diff <- apply(Zeta_hat_mcmc - Reduce(rbind, zeta), 1, function(x) sqrt(cprod(x)))
rmse_mcmc <- sqrt(mean(norm_diff))

####################################################
#
#  RESULTS
#
####################################################

psi_names <- matrix(NA, 2, L)
for(l in 1:L) {
	
	psi_names[1, l] <- paste("psi_", l, "_vmp", sep="")
	psi_names[2, l] <- paste("psi_", l, "_mcmc", sep="")
}
psi_names <- as.vector(psi_names)
col_names <- c("N", "sim", "mu_vmp", "mu_mcmc", psi_names, "zeta_vmp", "zeta_mcmc")
n_col_acc <- length(col_names)
file_name <- paste0("fpca_acc_", task_id, ".txt")
write(col_names, file_name, ncol = n_col_acc, append = FALSE)
psi_acc <- as.vector(rbind(psi_vmp_acc, psi_mcmc_acc))
acc_res <- c(N, sim, mu_vmp_acc, mu_mcmc_acc, psi_acc, rmse_vmp, rmse_mcmc)
write(acc_res, file_name, ncol = n_col_acc, append = TRUE)

col_names <- c("N", "sim", "iter", "vmp", "mcmc")
n_col_speed <- length(col_names)
file_name <- paste0("fpca_speed_", task_id, ".txt")
write(col_names, file_name, ncol = n_col_speed, append = FALSE)
speed_res <- c(N, sim, iter, vmp_time, mcmc_time)
write(speed_res, file_name, ncol = n_col_speed, append=TRUE)

############ End of fpca_vmp_vs_mcmc.R ############