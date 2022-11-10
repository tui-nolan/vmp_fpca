

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

n_g <- 1000
n <- 100
beta <- c(-0.1, 0.5)
d <- length(beta)
x <- runif(n)
X <- X_design(x)
sigma_eps <- 0.1
sigsq_eps <- sigma_eps^2

epsilon <- rnorm(n, 0, sigma_eps)
y <- as.vector(X %*% beta) + epsilon

plot(x, y, pch = 16, cex = 0.4)

n_vmp <- 500

sigsq_beta <- 1e10
A <- 1e5

eta_vec <- vector("list", length = 12)
names(eta_vec) <- c(
	"beta->p(Y|beta,sigsq_eps)", "p(Y|beta,sigsq_eps)->beta",
	"beta->p(beta)", "p(beta)->beta",
	"sigsq_eps->p(Y|beta,sigsq_eps)", "p(Y|beta,sigsq_eps)->sigsq_eps",
	"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
	"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
	"a_eps->p(a_eps)", "p(a_eps)->a_eps"
)

G <- vector("list", length = 10)
names(G) <- c(
	"sigsq_eps->p(Y|beta,sigsq_eps)", "p(Y|beta,sigsq_eps)->sigsq_eps",
	"sigsq_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->sigsq_eps",
	"a_eps->p(sigsq_eps|a_eps)", "p(sigsq_eps|a_eps)->a_eps",
	"a_eps->p(a_eps)", "p(a_eps)->a_eps"
)

eta_1 <- cprod(X, y)
eta_2 <- -0.5*cprod(duplication.matrix(d), as.vector(crossprod(X)))
eta_vec$"p(Y|beta,sigsq_eps)->beta" <- c(eta_1, eta_2)

eta_vec$"p(beta)->beta" <- gauss_prior_frag(rep(0, d), sigsq_beta*diag(d), use_vech = TRUE)

eta_1 <- -d/2
eta_2 <- -d/2
eta_vec$"p(Y|beta,sigsq_eps)->sigsq_eps" <- c(eta_1, eta_2)
G$"p(Y|beta,sigsq_eps)->sigsq_eps" <- "full"

eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- c(-3/2, -1/2)
G$"p(sigsq_eps|a_eps)->sigsq_eps" <- "full"

eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- c(-1/2, -1/2)
G$"p(sigsq_eps|a_eps)->a_eps" <- "diag"

igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))

eta_vec$"p(a_eps)->a_eps" <- igw_prior_updates[[2]]
G$"p(a_eps)->a_eps" <- igw_prior_updates[[1]]

for(iter in 1:n_vmp) {
	
	print(iter)
	
	eta_vec$"beta->p(Y|beta,sigsq_eps)" <- eta_vec$"p(beta)->beta"
	eta_vec$"beta->p(beta)" <- eta_vec$"p(Y|beta,sigsq_eps)->beta"
	
	eta_vec$"sigsq_eps->p(Y|beta,sigsq_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps"
	G$"sigsq_eps->p(Y|beta,sigsq_eps)" <- G$"p(sigsq_eps|a_eps)->sigsq_eps"
	eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(Y|beta,sigsq_eps)->sigsq_eps"
	G$"sigsq_eps->p(sigsq_eps|a_eps)" <- G$"p(Y|beta,sigsq_eps)->sigsq_eps"
	
	eta_vec$"a_eps->p(sigsq_eps|a_eps)" <- eta_vec$"p(a_eps)->a_eps"
	G$"a_eps->p(sigsq_eps|a_eps)" <- G$"p(a_eps)->a_eps"
	eta_vec$"a_eps->p(a_eps)" <- eta_vec$"p(sigsq_eps|a_eps)->a_eps"
	G$"a_eps->p(a_eps)" <- G$"p(sigsq_eps|a_eps)->a_eps"
	
	# Update p(Y|beta,sigsq) fragment:
	
	eta_in <- list(
		eta_vec$"beta->p(Y|beta,sigsq_eps)",
		eta_vec$"p(Y|beta,sigsq_eps)->beta",
		eta_vec$"sigsq_eps->p(Y|beta,sigsq_eps)",
		eta_vec$"p(Y|beta,sigsq_eps)->sigsq_eps"
	)
	
	G_in <- list(
		G$"sigsq_eps->p(Y|beta,sigsq_eps)",
		G$"p(Y|beta,sigsq_eps)->sigsq_eps"
	)
	
	gauss_lik_fragment <- gauss_lik_frag(
		eta_in, G_in, y, X, sigsq_beta
	)
	
	eta_vec$"p(Y|beta,sigsq_eps)->beta" <- gauss_lik_fragment$"eta"[[1]]
	eta_vec$"p(Y|beta,sigsq_eps)->sigsq_eps" <- gauss_lik_fragment$"eta"[[2]]
	
	G$"p(Y|beta,sigsq_eps)->sigsq_eps" <- gauss_lik_fragment$"G"[[1]]
	
	# Update p(sigsq_eps|a_eps) fragment:
	
	eta_in <- list(
		eta_vec$"sigsq_eps->p(sigsq_eps|a_eps)",
		eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps",
		eta_vec$"a_eps->p(sigsq_eps|a_eps)",
		eta_vec$"p(sigsq_eps|a_eps)->a_eps"
	)
	
	iter_igw_fragment <- iter_igw_frag(
		eta_in, G$"a_eps->p(sigsq_eps|a_eps)",
		1, G$"sigsq_eps->p(sigsq_eps|a_eps)"
	)
	
	eta_vec$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"eta"[[1]]
	eta_vec$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"eta"[[2]]
	
	G$"p(sigsq_eps|a_eps)->sigsq_eps" <- iter_igw_fragment$"G"[[1]]
	G$"p(sigsq_eps|a_eps)->a_eps" <- iter_igw_fragment$"G"[[2]]
}

# Save the original q_beta:

eta_beta <- list(eta_vec$"p(Y|beta,sigsq_eps)->beta", eta_vec$"p(beta)->beta")
q_beta <- gauss_q(eta_beta, use_vech = TRUE)
E_q_beta <- q_beta[[1]]
Cov_q_beta <- q_beta[[2]]

x_g <- seq(0, 1, length = n_g)
X_g <- X_design(x_g)
y_g <- X_g %*% E_q_beta

lines(x_g, y_g, col = "red")

sd_vec <- diag(tcrossprod(X_g %*% Cov_q_beta, X_g))

lines(x_g, y_g + qnorm(0.025)*sd_vec, lty = 2, col = "red")
lines(x_g, y_g + qnorm(0.975)*sd_vec, lty = 2, col = "red")












