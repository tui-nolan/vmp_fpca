######### R script: fpca_vmp_vs_mfvb.R ##########

# For comparing a simple functional principal
# components analysis via MFVB and VMP.

# Created: 31 JAN 2022
# Last changed: 02 FEB 2022

# Updates:
# 1. Estimating eigenfunctions only.
# 2. Establishing the multiple penalization fragment.
# 3. Including scores and orthogonal decomposition.
# 4. Including mean function and extending orthogonal decomposition.

# Load libraries:

library(MASS)
library(magic)
library(lattice)
library(pracma)

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
source("vmp_functions.r")
source("fpca_algs.r")

setwd("..")

# Establish simulation variables:

N <- 50                              # number of curves
n_sample <- 4                       # number of curves for the plots
N_sample <- sort(sample(1:N, n_sample))   # specific curves for the plots
T_vec <- round(runif(N, 50, 70))    # number of time observations for each curve
n_int_knots <- 23                   # number of interior knots
K <- n_int_knots + 2                # number of spline basis functions
L <- 2                              # number of FPCA basis functions
data_col <- "grey51"                # colour of the data in the plots

n_vmp <- 100                        # number of VMP iterations
n_mfvb <- 100                       # number of MFVB iterations
n_g <- 1000                         # length of the plotting grid
vmp_col <- "red"                    # colour of the VMP plots
mfvb_col <- "blue"                  # colour of the MFVB plots
d <- (K+2)*L                        # dimension of spline vector

sigma_zeta <- c(3, 1)               # sd for first and second scores
sigma_eps <- 1                      # sd of the residuals
sigsq_eps <- sigma_eps^2

# Set up plot-grid dimensions:

plot_dim <- c(2, 2)                 # (ncol, nrow) for curve plots
plot_bf_dim <- c(1, 2)              # (ncol, nrow) for basis function plots

plot_width <- 7
plot_height <- 14

construct_pdf <- FALSE              # save the plots in a PDF?

# Set the FPCA basis functions:

psi_1 <- function(t) return(sqrt(2)*sin(2*pi*t))
psi_2 <- function(t) return(sqrt(2)*cos(2*pi*t))

# Set up fixed parameters:

time_obs <- sapply(T_vec, runif)
time_obs <- lapply(time_obs, sort)

unique_time_obs <- sort(unique(Reduce(c, time_obs)))
int_knots <- quantile(unique_time_obs, seq(0,1,length=K)[-c(1,K)])

X <- vector("list", length=N)
Z <- vector("list", length=N)
C <- vector("list", length=N)
for(i in 1:N) {
	
	X[[i]] <- X_design(time_obs[[i]])
	Z[[i]] <- ZOSull(time_obs[[i]], range.x=c(0, 1), intKnots=int_knots)
	C[[i]] <- cbind(X[[i]], Z[[i]])
}

# Set up parameters for estimation:

zeta <- vector("list", length=N)
for(i in 1:N) {
	
	zeta[[i]] <- rep(NA, L)
	
	for(l in 1:L) {
		
		zeta[[i]][l] <- rnorm(1, 0, sigma_zeta[l])
	}
}

# Set up curve observations:

Psi_t <- vector("list", length=N)
Y <- vector("list", length=N)
for(i in 1:N) {
	
	Psi_t[[i]] <- matrix(NA, nrow=T_vec[[i]], ncol=L)
	for(l in 1:L) {
		
		char_val <- paste("psi_", l, "(time_obs[[i]])", sep="")
		Psi_t[[i]][,l] <- eval(parse(text=char_val))
	}
	
	epsilon <- rnorm(T_vec[[i]], 0, sigma_eps)
	
	Y_hat <- as.vector(Psi_t[[i]]%*%zeta[[i]])
	Y[[i]] <- Y_hat + epsilon
}

# Plot the data:

Y_sample <- Y[N_sample]
time_obs_sample <- time_obs[N_sample]
T_sample <- sapply(Y_sample, length)

Y_vec <- Reduce(c, Y_sample)
time_vec <- Reduce(c, time_obs_sample)

curve_labels <- vector("list", length = n_sample)
factor_id <- rep(NA, n_sample)
for(i in 1:n_sample) {
	
	N_i <- N_sample[i]
	factor_id[i] <- parse(text=paste("Y[", N_i, "] (t)", sep=""))
	curve_val <- eval(bquote(expression(Y[.(N_i)] (t))))
	curve_labels[[i]] <- rep(curve_val, T_vec[N_i])
}
curve_labels <- do.call(c, curve_labels)
curve_labels <- factor(curve_labels, levels=factor_id)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	fl <- factor_id
		
	strip.default(which.given,which.panel,var.name,fl,...)
}

raw_data_plots <- xyplot(
	Y_vec ~ time_vec | curve_labels, groups=curve_labels,
	data=data.frame(
		time_vec=time_vec, Y_vec=Y_vec,
		curve_labels=curve_labels
	),
	layout=plot_dim, main="",
	strip=strip.math,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1)),
	xlab="time",
	ylab="nonlinear curves",
	as.table=TRUE,
	panel=function(x, y, subscripts, groups) {
		
		iPan <- panel.number()
		i <- rep(N_sample, each=1)[iPan]
		panel.grid()
		panel.superpose(
			x[order(x)], y[order(x)], subscripts, groups,
			type="p", col=data_col, pch=16, cex=0.4
		)
	}
)

print(raw_data_plots)

wait()

# Set up plotting grid

time_g <- seq(0,1, length.out=n_g)

X_g <- X_design(time_g)
Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
C_g <- cbind(X_g, Z_g)

Psi_g <- matrix(NA, n_g, L)
for(l in 1:L) {
	
	char_val <- paste("psi_", l, "(time_g)", sep="")
	Psi_g[,l] <- eval(parse(text=char_val))
}

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
mu_beta <- rep(0, 2)
A <- 1e5
sigsq_zeta <- 1

####################################################
#
#  VMP  SIMULATIONS
#
####################################################

mu_q_recip_sigsq_eps <- 1/sigsq_eps

mu_q_zeta <- vector("list", length=N)
Sigma_q_zeta <- vector("list", length=N)
for(i in 1:N) {
	
	mu_q_zeta[[i]] <- rnorm(L, 0, sigma_zeta)
	Sigma_q_zeta[[i]] <- diag(L)
}

eta_vec <- vector("list", length=16)
names(eta_vec) <- c(
	"nu_p->p(Y|nu_p,zeta)", "p(Y|nu_p,zeta)->nu_p",
	"zeta->p(Y|nu_p,zeta)", "p(Y|nu_p,zeta)->zeta",
	"zeta->p(zeta)", "p(zeta)->zeta",
	"nu_p->p(nu_p|sigsq_p)", "p(nu_p|sigsq_p)->nu_p",
	"sigsq_p->p(nu_p|sigsq_p)", "p(nu_p|sigsq_p)->sigsq_p",
	"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
	"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
	"a_p->p(a_p)", "p(a_p)->a_p"
)

G <- vector("list", length=8)
names(G) <- c(
	"sigsq_p->p(nu_p|sigsq_p)", "p(nu_p|sigsq_p)->sigsq_p",
	"sigsq_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->sigsq_p",
	"a_p->p(sigsq_p|a_p)", "p(sigsq_p|a_p)->a_p",
	"a_p->p(a_p)", "p(a_p)->a_p"
)

eta_1_sum <- 0
eta_2_sum <- 0
for(i in 1:N) {
	
	M_q_zeta_zeta_T <- Sigma_q_zeta[[i]] + tcrossprod(mu_q_zeta[[i]])
	
	sum_val <- cprod(kronecker(t(mu_q_zeta[[i]]), C[[i]]), Y[[i]])
	eta_1_sum <- eta_1_sum + sum_val
	
	sum_val <- as.vector(kronecker(M_q_zeta_zeta_T, crossprod(C[[i]])))
	eta_2_sum <- eta_2_sum + sum_val
}
eta_1 <- 1*eta_1_sum
eta_2 <- -1/2*eta_2_sum
eta_vec$"p(Y|nu_p,zeta)->nu_p" <- c(eta_1, eta_2)

D_L <- duplication.matrix(L)
eta_1 <- Reduce(cbind, mu_q_zeta)
eta_2 <- replicate(N, -0.5*cprod(D_L, as.vector(diag(L) - 1/sigsq_zeta*diag(L))))
eta_vec$"p(Y|nu_p,zeta)->zeta" <- rbind(eta_1, eta_2)

eta_vec$"p(zeta)->zeta" <- replicate(
	N,
	gauss_prior_frag(rep(0, L), 1/sigsq_zeta*diag(L), use_vech=TRUE)
)

eta_1 <- rep(0, d)
eta_2 <- -0.5*as.vector(diag(d))
eta_vec$"p(nu_p|sigsq_p)->nu_p" <- c(eta_1, eta_2)

eta_vec$"p(nu_p|sigsq_p)->sigsq_p" <- replicate(L, c(-K/2, -K/2))
G$"p(nu_p|sigsq_p)->sigsq_p" <- rep("full", L)

eta_vec$"p(sigsq_p|a_p)->sigsq_p" <- replicate(L, c(-3/2, -1/2))
G$"p(sigsq_p|a_p)->sigsq_p" <- rep("full", L)

eta_vec$"p(sigsq_p|a_p)->a_p" <- replicate(L, c(-1/2, -1/2))
G$"p(sigsq_p|a_p)->a_p" <- rep("diag", L)

igw_prior_updates <- igw_prior_frag(list("diag", 1, 1/A^2))

eta_vec$"p(a_p)->a_p" <- replicate(L, igw_prior_updates[[2]])
G$"p(a_p)->a_p" <- rep(igw_prior_updates[[1]], L)

for(i_iter in 1:n_vmp) {
	
	cat("starting iteration", i_iter, "of", n_vmp, "\n")
	
	eta_vec$"nu_p->p(Y|nu_p,zeta)" <- eta_vec$"p(nu_p|sigsq_p)->nu_p"
	eta_vec$"nu_p->p(nu_p|sigsq_p)" <- eta_vec$"p(Y|nu_p,zeta)->nu_p"
	
	eta_vec$"zeta->p(Y|nu_p,zeta)" <- eta_vec$"p(zeta)->zeta"
	eta_vec$"zeta->p(zeta)" <- eta_vec$"p(Y|nu_p,zeta)->zeta"
	
	eta_vec$"sigsq_p->p(nu_p|sigsq_p)" <- eta_vec$"p(sigsq_p|a_p)->sigsq_p"
	G$"sigsq_p->p(nu_p|sigsq_p)" <- G$"p(sigsq_p|a_p)->sigsq_p"
	eta_vec$"sigsq_p->p(sigsq_p|a_p)" <- eta_vec$"p(nu_p|sigsq_p)->sigsq_p"
	G$"sigsq_p->p(sigsq_p|a_p)" <- G$"p(nu_p|sigsq_p)->sigsq_p"
	
	eta_vec$"a_p->p(sigsq_p|a_p)" <- eta_vec$"p(a_p)->a_p"
	G$"a_p->p(sigsq_p|a_p)" <- G$"p(a_p)->a_p"
	eta_vec$"a_p->p(a_p)" <- eta_vec$"p(sigsq_p|a_p)->a_p"
	G$"a_p->p(a_p)" <- G$"p(sigsq_p|a_p)->a_p"
	
	# Update p(Y|nu_p,zeta) fragment:
	
	eta_in <- list(
		eta_vec$"nu_p->p(Y|nu_p,zeta)",
		eta_vec$"p(Y|nu_p,zeta)->nu_p",
		eta_vec$"zeta->p(Y|nu_p,zeta)",
		eta_vec$"p(Y|nu_p,zeta)->zeta"
	)
	
	eta_in_nu_p <- list(eta_in[[1]], eta_in[[2]])
	q_nu_psi <- gauss_q(eta_in_nu_p, use_vech=FALSE)
	mu_q_nu_psi <- q_nu_psi[[1]]
	Sigma_q_nu_psi <- q_nu_psi[[2]]
	
	eta_1_sum <- 0
	eta_2_sum <- 0
	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		eta_in_zeta <- list(
			eta_in[[3]][,i],
			eta_in[[4]][,i]
		)
		q_zeta <- gauss_q(eta_in_zeta, use_vech=TRUE)
		mu_q_zeta[[i]] <- q_zeta[[1]]
		Sigma_q_zeta[[i]] <- q_zeta[[2]]
		M_q_zeta_zeta_T <- Sigma_q_zeta[[i]] + tcrossprod(mu_q_zeta[[i]])
		
		sum_val <- cprod(kronecker(t(mu_q_zeta[[i]]), C[[i]]), Y[[i]])
		eta_1_sum <- eta_1_sum + sum_val
		
		sum_val <- as.vector(kronecker(M_q_zeta_zeta_T, crossprod(C[[i]])))
		eta_2_sum <- eta_2_sum + sum_val
	}
	eta_1 <- mu_q_recip_sigsq_eps*eta_1_sum
	eta_2 <- -0.5*mu_q_recip_sigsq_eps*eta_2_sum
	eta_nu_psi <- c(eta_1, eta_2)
	
	groups <- ceiling((1:d)/(K+2))-1
	inds <- split(1:d, groups)
	
	M_q_H_psi <- vector("list", length=N)
	for(i in 1:N) {
		
		M_q_H_psi[[i]] <- matrix(NA, L, L)
		
		for(l_1 in 1:L) {
			
			inds_1 <- inds[[l_1]]
			mu_q_nu_psi_l_1 <- mu_q_nu_psi[inds_1]
			
			for(l_2 in 1:L) {
				
				inds_2 <- inds[[l_2]]
				mu_q_nu_psi_l_2 <- mu_q_nu_psi[inds_2]
				
				Cov_q_nu_psi_l_21 <- Sigma_q_nu_psi[inds_2, inds_1]
				
				tr_term <- tr(Cov_q_nu_psi_l_21%*%crossprod(C[[i]]))
				cprod_term <- cprod(mu_q_nu_psi_l_1, crossprod(C[[i]])%*%mu_q_nu_psi_l_2)
				M_q_H_psi[[i]][l_1, l_2] <- tr_term + cprod_term
			}
		}
	}
	
	M_q_V_psi <- matrix(NA, K+2, L)
	for(l in 1:L) {
		
		l_inds <- inds[[l]]
		M_q_V_psi[,l] <- mu_q_nu_psi[l_inds]
	}
	
	D_L <- duplication.matrix(L)
	
	eta_1 <- matrix(NA, L, N)
	eta_2 <- matrix(NA, 0.5*L*(L+1), N)
	for(i in 1:N) {
		
		M_q_Psi <- C[[i]]%*%M_q_V_psi
		eta_1[,i] <- mu_q_recip_sigsq_eps*cprod(M_q_Psi, Y[[i]])
		eta_2[,i] <- -0.5*mu_q_recip_sigsq_eps*cprod(D_L, as.vector(M_q_H_psi[[i]]))
	}
	eta_zeta <- rbind(eta_1, eta_2)
	
	eta_out <- list(eta_nu_psi, eta_zeta)
	names(eta_out) <- c(
		"p(Y|nu_p,zeta)->nu_p",
		"p(Y|nu_p,zeta)->zeta"
	)
	
	eta_vec$"p(Y|nu_p,zeta)->nu_p" <- eta_out[[1]]
	eta_vec$"p(Y|nu_p,zeta)->zeta" <- eta_out[[2]]
	
	# Update p(nu_p|sigsq_p):
	
	eta_in <- list(
		eta_vec$"nu_p->p(nu_p|sigsq_p)", eta_vec$"p(nu_p|sigsq_p)->nu_p",
		eta_vec$"sigsq_p->p(nu_p|sigsq_p)", eta_vec$"p(nu_p|sigsq_p)->sigsq_p"
	)
	
	G_in <- list(G$"sigsq_p->p(nu_p|sigsq_p)", G$"p(nu_p|sigsq_p)->sigsq_p")
	
	mult_gauss_pen_fragment <- mult_gauss_pen_frag(eta_in, G_in, L, sigsq_beta)
	
	eta_vec$"p(nu_p|sigsq_p)->nu_p" <- mult_gauss_pen_fragment$eta[[1]]
	eta_vec$"p(nu_p|sigsq_p)->sigsq_p" <- mult_gauss_pen_fragment$eta[[2]]
	
	G$"p(nu_p|sigsq_p)->sigsq_p" <- mult_gauss_pen_fragment$G[[1]]
	
	# Update p(sigsq_p|a_p) fragment:
	
	for(l in 1:L) {
		
		eta_in <- list(
			eta_vec$"sigsq_p->p(sigsq_p|a_p)"[,l],
			eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l],
			eta_vec$"a_p->p(sigsq_p|a_p)"[,l],
			eta_vec$"p(sigsq_p|a_p)->a_p"[,l]
		)
		
		iter_igw_fragment <- iter_igw_frag(
			eta_in, G$"a_p->p(sigsq_p|a_p)"[l],
			1, G$"sigsq_p->p(sigsq_p|a_p)"[l]
		)
		
		eta_vec$"p(sigsq_p|a_p)->sigsq_p"[,l] <- iter_igw_fragment$"eta"[[1]]
		eta_vec$"p(sigsq_p|a_p)->a_p"[,l] <- iter_igw_fragment$"eta"[[2]]
		
		G$"p(sigsq_p|a_p)->sigsq_p"[l] <- iter_igw_fragment$"G"[[1]]
		G$"p(sigsq_p|a_p)->a_p"[l] <- iter_igw_fragment$"G"[[2]]
	}
	
	# Construct orthogonal decomposition:
	
	eta_in <- list(
		eta_vec$"p(nu_p|sigsq_p)->nu_p", eta_vec$"p(Y|nu_p,zeta)->nu_p",
		eta_vec$"p(zeta)->zeta", eta_vec$"p(Y|nu_p,zeta)->zeta"
	)
	
	eta_nu_psi <- list(eta_in[[1]], eta_in[[2]])
	q_nu_psi <- gauss_q(eta_nu_psi, use_vech=FALSE)
	mu_q_nu_psi <- q_nu_psi[[1]]
	Sigma_q_nu_psi <- q_nu_psi[[2]]
	
	groups <- ceiling((1:d)/(K+2))-1
	M_q_V_psi <- Reduce(cbind, split(mu_q_nu_psi, groups))
	M_q_Psi <- C_g%*%M_q_V_psi
	
	#
	
	mu_q_zeta <- vector("list", length=N)
	Sigma_q_zeta <- vector("list", length=N)
	for(i in 1:N) {
		
		eta_zeta <- list(eta_in[[3]][,i], eta_in[[4]][,i])
		q_zeta <- gauss_q(eta_zeta, use_vech=TRUE)
		mu_q_zeta[[i]] <- q_zeta[[1]]
		Sigma_q_zeta[[i]] <- q_zeta[[2]]
	}
	
	M_q_Zeta <- Reduce(rbind, mu_q_zeta)
	
	#
	
	one_N <- rep(1, N)
	
	M_q_Psi_svd <- svd(M_q_Psi)
	U_orth <- M_q_Psi_svd$u
	D_diag <- diag(M_q_Psi_svd$d)
	V_orth <- M_q_Psi_svd$v
	
	M_q_Zeta_rotn <- M_q_Zeta%*%V_orth%*%D_diag
	eigen_M_q_Zeta_shift <- eigen(cov(M_q_Zeta_rotn))
	Q <- eigen_M_q_Zeta_shift$vectors
	Lambda <- diag(eigen_M_q_Zeta_shift$values + 1e-10)
	Lambda_inv <- diag(1/(eigen_M_q_Zeta_shift$values + 1e-10))
	S <- Q%*%sqrt(Lambda)
	S_inv <- tcrossprod(sqrt(Lambda_inv), Q)
	
	Psi_hat <- U_orth%*%S
	Zeta_hat <- tcrossprod(M_q_Zeta_rotn, S_inv)
	
	norm_const <- rep(NA, L)
	for(l in 1:L) {
		
		norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[,l])^2))
		if(norm_const[l]!=0) {
			
			Psi_hat[,l] <- Psi_hat[,l]/norm_const[l]
			Zeta_hat[,l] <- norm_const[l]*Zeta_hat[,l]
			
			cprod_sign <- sign(cprod(Psi_hat[,l], Psi_g[,l]))
			if(cprod_sign==-1) {
				
				Psi_hat[,l] <- -Psi_hat[,l]
				Zeta_hat[,l] <- -Zeta_hat[,l]
			}
		}
	}
	
	#
	
	scale_mat <- diag(norm_const)
	for(i in 1:N) {
		
		mat_transform <- S_inv%*%tcrossprod(D_diag, V_orth)
		Sigma_q_zeta[[i]] <- tcrossprod(mat_transform%*%Sigma_q_zeta[[i]], mat_transform)
		Sigma_q_zeta[[i]] <- tcrossprod(scale_mat%*%Sigma_q_zeta[[i]], scale_mat)
	}
	
	#
	
	zeta_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_mean <- Zeta_hat[i,][1:2]
		
		zeta_ellipse <- ellipse(
			Sigma_q_zeta[[i]][1:2, 1:2],
			centre=zeta_mean,
			level=0.95
		)
		
		zeta_summary[[i]] <- list(zeta_mean, zeta_ellipse)
		names(zeta_summary[[i]]) <- c("mean", "credible boundary")
	}
}

Psi_hat_vmp <- Psi_hat
mu_q_zeta_vmp <- mu_q_zeta
Sigma_q_zeta_vmp <- Sigma_q_zeta
M_q_V_psi_vmp <- M_q_V_psi
zeta_summary_vmp <- zeta_summary

Y_summary_vmp <- vector("list", length = N)
for(i in 1:N) {
	
	sd_vec <- sqrt(diag(tcrossprod(Psi_hat%*%Sigma_q_zeta[[i]], Psi_hat)))
	Y_hat <- as.vector(C_g %*% M_q_V_psi %*% mu_q_zeta[[i]])
	
	Y_summary_vmp[[i]] <- matrix(NA, nrow=n_g, ncol=3)
	Y_summary_vmp[[i]][,1] <- Y_hat + qnorm(0.025)*sd_vec
	Y_summary_vmp[[i]][,2] <- Y_hat
	Y_summary_vmp[[i]][,3] <- Y_hat + qnorm(0.975)*sd_vec
}

####################################################
#
#  MFVB  SIMULATIONS
#
####################################################

E_q_zeta <- vector("list", length = N)
Cov_q_zeta <- vector("list", length = N)
for(i in 1:N) {
	
	E_q_zeta[[i]] <- rnorm(L, 0, 1)
	Cov_q_zeta[[i]] <- diag(L)
}

E_q_inv_Sigma_psi <- vector("list", length = L)
for(l in 1:L) {
	
	E_q_inv_Sigma_psi[[l]] <- adiag(1/sigsq_beta*diag(2), diag(K))
}
E_q_inv_Sigma_psi <- Reduce(adiag, E_q_inv_Sigma_psi)

E_q_recip_a <- rep(1, L)

for(i_iter in 1:n_mfvb) {
	
	cat("starting iteration", i_iter, "of", n_mfvb, "\n")
	
	# Update q(nu_psi):
	
	Cov_sum <- 0
	E_sum <- 0
	for(i in 1:N) {
		
		E_q_zeta_zetaT <- Cov_q_zeta[[i]] + tcrossprod(E_q_zeta[[i]])
		
		sum_val <- kronecker(E_q_zeta_zetaT, crossprod(C[[i]]))
		Cov_sum <- Cov_sum + sum_val
		
		sum_val <- cprod(kronecker(t(E_q_zeta[[i]]), C[[i]]), Y[[i]])
		E_sum <- E_sum + sum_val
	}
	Cov_q_nu_psi <- solve(1/sigsq_eps*Cov_sum + E_q_inv_Sigma_psi)
	E_q_nu_psi <- 1/sigsq_eps*as.vector(Cov_q_nu_psi %*% E_sum)
	E_q_V_psi <- matrix(E_q_nu_psi, K + 2, L)
	
	inds <- matrix(1:(L*(K + 2)), K + 2, L)
	u_inds <- inds[-c(1, 2), ]
	E_q_u_psi <- vector("list", length = L)
	Cov_q_u_psi <- vector("list", length = L)
	for(l in 1:L) {
		
		E_q_u_psi[[l]] <- E_q_nu_psi[u_inds[,l]]
		Cov_q_u_psi[[l]] <- Cov_q_nu_psi[u_inds[,l], u_inds[,l]]
	}
	
	E_q_H_psi <- vector("list", length=N)
	for(i in 1:N) {
		
		E_q_H_psi[[i]] <- matrix(NA, L, L)
		
		for(l_1 in 1:L) {
			
			inds_1 <- inds[, l_1]
			E_q_nu_psi_l_1 <- E_q_nu_psi[inds_1]
			
			for(l_2 in 1:L) {
				
				inds_2 <- inds[, l_2]
				E_q_nu_psi_l_2 <- E_q_nu_psi[inds_2]
				
				Cov_q_nu_psi_l_21 <- Cov_q_nu_psi[inds_2, inds_1]
				
				tr_term <- tr(Cov_q_nu_psi_l_21%*%crossprod(C[[i]]))
				cprod_term <- cprod(E_q_nu_psi_l_1, crossprod(C[[i]])%*%E_q_nu_psi_l_2)
				E_q_H_psi[[i]][l_1, l_2] <- tr_term + cprod_term
			}
		}
	}
	
	# For i = 1, ..., N, update q(zeta[[i]]):
	
	E_q_zeta <- vector("list", length = N)
	Cov_q_zeta <- vector("list", length = N)
	for(i in 1:N) {
		
		E_q_Psi <- C[[i]] %*% E_q_V_psi
		
		Cov_q_zeta[[i]] <- solve(1/sigsq_eps*E_q_H_psi[[i]] + diag(L))
		E_q_zeta[[i]] <- 1/sigsq_eps*as.vector(Cov_q_zeta[[i]] %*% cprod(E_q_Psi, Y[[i]]))
	}
	
	# For l = 1, ..., L, update q(sigsq_psi[l]):
	
	lambda_q_sigsq_psi <- rep(NA, L)
	E_q_recip_sigsq_psi <- rep(NA, L)
	E_q_inv_Sigma_psi <- vector("list", length = L)
	for(l in 1:L) {
		
		tr_term <- tr(Cov_q_u_psi[[l]])
		cprod_term <- cprod(E_q_u_psi[[l]])
		lambda_q_sigsq_psi[l] <- tr_term + cprod_term + E_q_recip_a[l]
		
		E_q_recip_sigsq_psi[l] <- (K + 1)/lambda_q_sigsq_psi[l]
		
		E_q_inv_Sigma_psi[[l]] <- adiag(1/sigsq_beta*diag(2), E_q_recip_sigsq_psi[l]*diag(K))
	}
	E_q_inv_Sigma_psi <- Reduce(adiag, E_q_inv_Sigma_psi)
	
	# For l = 1, ..., L, updates q(a_psi[l]):
	
	lambda_q_a_psi <- rep(NA, L)
	E_q_recip_a_psi <- rep(NA, L)
	for(l in 1:L) {
		
		lambda_q_a_psi[l] <- E_q_recip_sigsq_psi[l] + 1/A^2
		
		E_q_recip_a_psi[l] <- 2/lambda_q_a_psi[l]
	}
	
	# Construct orthogonal decomposition:
	
	E_q_Psi <- C_g %*% E_q_V_psi
	E_q_Zeta <- Reduce(rbind, E_q_zeta)
	
	one_N <- rep(1, N)
	
	E_q_Psi_svd <- svd(E_q_Psi)
	U_psi <- E_q_Psi_svd$u
	D_psi <- diag(E_q_Psi_svd$d)
	V_psi <- E_q_Psi_svd$v
	
	E_q_Zeta_rotn <- E_q_Zeta %*% V_psi %*% D_psi
	eigen_E_q_Zeta_shift <- eigen(cov(E_q_Zeta_rotn))
	Q <- eigen_E_q_Zeta_shift$vectors
	Lambda <- diag(eigen_E_q_Zeta_shift$values + 1e-10)
	Lambda_inv <- diag(1/(eigen_E_q_Zeta_shift$values + 1e-10))
	S <- Q %*% sqrt(Lambda)
	S_inv <- tcrossprod(sqrt(Lambda_inv), Q)
	
	Psi_hat <- U_psi %*% S
	Zeta_hat <- tcrossprod(E_q_Zeta_rotn, S_inv)
	
	norm_const <- rep(NA, L)
	for(l in 1:L) {
		
		norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[,l])^2))
		if(norm_const[l]!=0) {
			
			Psi_hat[,l] <- Psi_hat[,l]/norm_const[l]
			Zeta_hat[,l] <- norm_const[l]*Zeta_hat[,l]
			
			cprod_sign <- sign(cprod(Psi_hat[,l], Psi_g[,l]))
			if(cprod_sign==-1) {
				
				Psi_hat[,l] <- -Psi_hat[,l]
				Zeta_hat[,l] <- -Zeta_hat[,l]
			}
		}
	}
	
	scale_mat <- diag(norm_const)
	Cov_zeta_hat <- vector("list", length = N)
	for(i in 1:N) {
		
		mat_transform <- S_inv %*% tcrossprod(D_psi, V_psi)
		Cov_zeta_hat[[i]] <- tcrossprod(mat_transform %*% Cov_q_zeta[[i]], mat_transform)
		Cov_zeta_hat[[i]] <- tcrossprod(scale_mat %*% Cov_zeta_hat[[i]], scale_mat)
	}
	
	zeta_summary <- vector("list", length=N)
	for(i in 1:N) {
		
		zeta_mean <- Zeta_hat[i,][1:2]
		
		zeta_ellipse <- ellipse(
			Cov_zeta_hat[[i]][1:2, 1:2],
			centre = zeta_mean,
			level = 0.95
		)
		
		zeta_summary[[i]] <- list(zeta_mean, zeta_ellipse)
		names(zeta_summary[[i]]) <- c("mean", "credible boundary")
	}
}

Psi_hat_mfvb <- Psi_hat
E_q_zeta_mfvb <- E_q_zeta
Cov_q_zeta_mfvb <- Cov_q_zeta
E_q_V_psi_mfvb <- E_q_V_psi
zeta_summary_mfvb <- zeta_summary

Y_summary_mfvb <- vector("list", length = N)
for(i in 1:N) {
	
	sd_vec <- sqrt(diag(tcrossprod(Psi_hat_mfvb %*% Cov_zeta_hat[[i]], Psi_hat)))
	Y_hat <- as.vector(C_g %*% E_q_V_psi %*% E_q_zeta[[i]])
	
	Y_summary_mfvb[[i]] <- matrix(NA, nrow=n_g, ncol=3)
	Y_summary_mfvb[[i]][,1] <- Y_hat + qnorm(0.025)*sd_vec
	Y_summary_mfvb[[i]][,2] <- Y_hat
	Y_summary_mfvb[[i]][,3] <- Y_hat + qnorm(0.975)*sd_vec
}

# Set up the curve plots:

fitted_data_plots <- xyplot(
	Y_vec ~ time_vec | curve_labels, groups=curve_labels,
	data=data.frame(
		time_vec=time_vec, Y_vec=Y_vec,
		curve_labels=curve_labels
	),
	layout=plot_dim, main="",
	strip=strip.math,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1)),
	xlab="time",
	ylab="nonlinear curves",
	as.table=TRUE,
	panel=function(x, y, subscripts, groups) {
		
		iPan <- panel.number()
		i <- rep(N_sample, each=1)[iPan]
		panel.grid()
		panel.superpose(
			x[order(x)], y[order(x)], subscripts, groups,
			type="p", col=data_col, pch=16, cex=0.4
		)
		panel.xyplot(
			time_g, Y_summary_vmp[[i]][,2],
			col=vmp_col, type="l", lwd=1
		)
		panel.xyplot(
			time_g, Y_summary_vmp[[i]][,1],
			col=vmp_col, type="l", lwd=1, lty = 2
		)
		panel.xyplot(
			time_g, Y_summary_vmp[[i]][,3],
			col=vmp_col, type="l", lwd=1, lty = 2
		)
		panel.xyplot(
			time_g, Y_summary_mfvb[[i]][,2],
			col=mfvb_col, type="l", lwd=1
		)
		panel.xyplot(
			time_g, Y_summary_mfvb[[i]][,1],
			col=mfvb_col, type="l", lwd=1, lty = 2
		)
		panel.xyplot(
			time_g, Y_summary_mfvb[[i]][,3],
			col=mfvb_col, type="l", lwd=1, lty = 2
		)
	}
)

print(fitted_data_plots)

wait()

# Set up the basis function plots:

time_g_psi <- rep(time_g, L)
psi_g_vec <- as.vector(Psi_g)

bf_labels <- vector("list", length=L)
bf_id <- rep(NA, L)
for(l in 1:L) {
	
	bf_id[l] <- parse(text=paste("psi[", l, "] (t)", sep=""))
	bf_val <- eval(bquote(expression(psi[.(l)] (t))))
	bf_labels[[l]] <- rep(bf_val, n_g)
}
bf_labels <- do.call(c, bf_labels)
bf_labels <- factor(bf_labels, levels=bf_id)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	fl <- bf_id
		
	strip.default(which.given,which.panel,var.name,fl,...)
}

bf_plots <- xyplot(
	psi_g_vec ~ time_g_psi | bf_labels, groups=bf_labels,
	data=data.frame(
		time_g_psi=time_g_psi, psi_g_vec=psi_g_vec,
		bf_labels=bf_labels
	),
	layout=plot_bf_dim, main="",
	strip=strip.math,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1)),
	xlab="time",
	ylab="basis functions",
	as.table=TRUE,
	panel=function(x, y, subscripts, groups) {
		
		lPan <- panel.number()
		l <- rep(1:L, each=1)[lPan]
		panel.grid()
		panel.superpose(
			x[order(x)], y[order(x)], subscripts, groups,
			type="l", col=data_col, pch=16, cex=0.4
		)
		panel.xyplot(
			time_g, Psi_hat_vmp[,l],
			col=vmp_col, type="l", lwd=1
		)
		panel.xyplot(
			time_g, Psi_hat_mfvb[,l],
			col=mfvb_col, type="l", lwd=1
		)
	}
)

print(bf_plots)

wait()

# Plot the scores:

plot_fpca_scores(N_sample, zeta_summary_vmp, zeta, plot_dim, data_col, vmp_col)

wait()

plot_fpca_scores(N_sample, zeta_summary_mfvb, zeta, plot_dim, data_col, vmp_col)

############ End of fpca_vmp_vs_mfvb.R ############