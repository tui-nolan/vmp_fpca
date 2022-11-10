######### R script: mlfpca_vmp.R ##########

# For comparing a MlFPCA via MFVB and MCMC

# Created: 09 JUL 2022
# Last changed: 09 JUL 2022

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
source("trapint.r")
source("cprod.r")
source("wait.r")
source("ise.r")
source("vmp_functions.r")
source("fourier_basis.r")
source("fpca_algs.r")

setwd("..")

set.seed(36)

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

n_mfvb <- 200                             # number of MFVB iterations
n_vmp <- 200                              # number of VMP iterations
n_g <- 1000                               # length of the plotting grid
mfvb_col <- "red"                         # colour of the MFVB plots
vmp_col <- "deepskyblue2"                 # colour of the VMP plots
d <- (K+2)*(L_1 + L_2 + 1)                # dimension of spline vector

sigma_zeta_1 <- 1/(1:L_1)                   # vector of st. dev.'s for the first level scores
sigma_zeta_2 <- 1/(1:L_2)                 # vector of st. dev.'s for the second level scores
l_zeta <- L_1 + M*L_2                     # length of score vector

sigma_eps <- 1                            # sd of the residuals
sigsq_eps <- sigma_eps^2

# Set up plot-grid dimensions:

plot_dim <- c(m_sample, n_sample)         # (ncol, nrow) for curve plots
plot_gbl_dim <- c(L + 1, 1)               # (ncol, nrow) for basis function plots
plot_score_1_dim <- c(2, 2)
plot_score_2_dim <- c(m_sample, n_sample)

plot_width <- 7
plot_height <- 14

construct_pdf <- FALSE              # save the plots in a PDF?

# Set the FPCA basis functions:

mu <- function(t) return(3*sin(pi*t) - 3/2)
Psi_1 <- fourier_basis(L_1)
Psi_2 <- fourier_basis(L_1 + L_2)[(L_1 + 1):(L_1 + L_2)]

# psi_11 <- function(t) return(sqrt(2)*sin(2*pi*t))
# psi_12 <- function(t) return(sqrt(2)*cos(2*pi*t))
# Psi_func_1 <- list(psi_11, psi_12)

# psi_21 <- function(t) return(sqrt(2)*sin(3*pi*t))
# psi_22 <- function(t) return(sqrt(2)*cos(3*pi*t))
# Psi_func_2 <- list(psi_21, psi_22)

# Generate the data:

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
X_g <- mlfpca_data$X_g
Z_g <- mlfpca_data$Z_g
C_g <- mlfpca_data$C_g
zeta_1 <- mlfpca_data$zeta_1
zeta_2 <- mlfpca_data$zeta_2
mu_g <- mlfpca_data$mu_g
Psi_1_g <- mlfpca_data$Psi_1_g
Psi_2_g <- mlfpca_data$Psi_2_g
Y <- mlfpca_data$Y

# YTY <- vector("list", length = N)
# CTY <- vector("list", length = N)
# CTC <- vector("list", length = N)
# for(i in 1:N) {
	
	# YTY[[i]] <- rep(NA, M[i])
	# CTY[[i]] <- matrix(NA, K + 2, M[i])
	# CTC[[i]] <- vector("list", length = M[i])
	# for(j in 1:M[i]) {
		
		# YTY[[i]][j] <- cprod(Y[[i]][[j]])
		# CTY[[i]][, j] <- cprod(C[[i]][[j]], Y[[i]][[j]])
		# CTC[[i]][[j]] <- crossprod(C[[i]][[j]])
	# }
# }

# Plot the data:

plot_mlfpca_data(N_sample, M_sample, time_obs, Y, plot_dim, data_col)

wait()

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
sigma_zeta <- 1
sigsq_zeta <- sigma_zeta^2
mu_beta <- rep(0, 2)
A <- 1e5

####################################################
#
#  VMP  SIMULATIONS
#
####################################################

cat("Starting the VMP algorithm. \n")

set.seed(1)

vmp_alg <- vmp_mlfpca(
	n_vmp, N, M, L_1, L_2, C, Y,
	sigsq_beta, A, criterion, plot_elbo = TRUE
)

eta_vec <- vmp_alg$"eta_vec"

# Orthogonalization:

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

wait()

####################################################
#
#  RESULTS
#
####################################################

# Plot the fits:

Y_sample <- vector("list", length = n_sample)
time_obs_sample <- vector("list", length = n_sample)
T_sample <- vector("list", length = n_sample)
Y_vec <- vector("list", length = n_sample)
time_vec <- vector("list", length = n_sample)
for(i in 1:n_sample) {
	
	N_i <- N_sample[i]
	Y_sample[[i]] <- Y[[N_i]][M_sample]
	time_obs_sample[[i]] <- time_obs[[N_i]][M_sample]
	T_sample[[i]] <- sapply(Y_sample[[i]], length)
	Y_vec[[i]] <- Reduce(c, Y_sample[[i]])
	time_vec[[i]] <- Reduce(c, time_obs_sample[[i]])
}
Y_vec <- Reduce(c, Y_vec)
time_vec <- Reduce(c, time_vec)

curve_labels_1 <- vector("list", length = n_sample)
for(i in 1:n_sample) {
	
	curve_labels_1[[i]] <- rep.int(M_sample, times =  T_sample[[i]])
}
curve_labels_1 <- Reduce(c, curve_labels_1)
curve_id_1 <- rep(NA, m_sample)
for(j in 1:m_sample) {
	
	curve_id_1[j] <- parse(text = as.character(M_sample[j]))
}
curve_labels_1 <- factor(curve_labels_1, levels=curve_id_1)

curve_labels_2 <- vector("list", length = n_sample)
curve_id_2 <- rep(NA, n_sample)
for(i in 1:n_sample) {
	
	N_i <- N_sample[i]
	curve_id_2[i] <- parse(text=paste("Y [", N_i, "] (t)", sep=""))
	curve_val <- eval(bquote(expression(Y[.(N_i)] (t))))
	curve_labels_2[[i]] <- rep(curve_val, sum(T_sample[[i]]))
}
curve_labels_2 <- do.call(c, curve_labels_2)
curve_labels_2 <- factor(curve_labels_2, levels=curve_id_2)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	if(which.given==1) {
		
		fl <- curve_id_1
		
		strip.default(which.given, which.panel, var.name, fl, ...)
	}
	
	if (which.given==2) {
		
		fl <- curve_id_2
		
		strip.default(which.given, which.panel, var.name, fl, ...)
	}
}

fitted_data_plots <- xyplot(
	Y_vec ~ time_vec | curve_labels_1*curve_labels_2, groups = curve_labels_1,
	data = data.frame(
		time_vec = time_vec, Y_vec = Y_vec,
		curve_labels_1 = curve_labels_1, curve_labels_2 = curve_labels_2
	), layout = plot_dim, main = "",
	strip=strip.math,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1)),
	xlab = "time", ylab = "functional responses", as.table = TRUE,
	panel=function(x,y,subscripts,groups) {
		
		i_pan <- panel.number()
		i <- rep(1:n_sample, each = m_sample)[i_pan]
		j <- rep(1:m_sample, n_sample)[i_pan]
		panel.grid()
		
		panel.superpose(
			x[order(x)], y[order(x)], subscripts, groups,
			type="p", col = data_col, pch = 16, cex = 0.4
		)
		
		panel.xyplot(
			time_g, vmp_res$"fits"[[i]][[j]][, 1],
			col = vmp_col, type = "l", lwd = 1, lty = 2
		)
		
		panel.xyplot(
			time_g, vmp_res$"fits"[[i]][[j]][, 2],
			col = vmp_col, type = "l", lwd = 1
		)
		
		panel.xyplot(
			time_g, vmp_res$"fits"[[i]][[j]][, 3],
			col = vmp_col, type = "l", lwd = 1, lty = 2
		)
	}
)

print(fitted_data_plots)

wait()

# Plot the global curves:

time_g_psi <- rep(time_g, L + 1)

gbl_g_vec <- c(mu_g, as.vector(Psi_1_g), as.vector(Psi_2_g))

mu_id <- parse(text = paste("mu (t)", sep=""))
mu_labels <- rep(mu_id, n_g)
mu_labels <- factor(mu_labels, levels = mu_id)

psi_1_labels <- vector("list", length = L_1)
psi_1_id <- rep(NA, L_1)
for(l in 1:L_1) {
	
	psi_1_id[l] <- parse(text = paste("psi[1", l, "] (t)", sep=""))
	psi_1_labels[[l]] <- rep(psi_1_id[l], n_g)
}
psi_1_labels <- do.call(c, psi_1_labels)
psi_1_labels <- factor(psi_1_labels, levels = psi_1_id)

psi_2_labels <- vector("list", length = L_2)
psi_2_id <- rep(NA, L_2)
for(l in 1:L_2) {
	
	psi_2_id[l] <- parse(text = paste("psi[2", l, "] (t)", sep=""))
	psi_2_labels[[l]] <- rep(psi_2_id[l], n_g)
}
psi_2_labels <- do.call(c, psi_2_labels)
psi_2_labels <- factor(psi_2_labels, levels = psi_2_id)

gbl_id <- c(mu_id, psi_1_id, psi_2_id)
gbl_labels <- c(mu_labels, psi_1_labels, psi_2_labels)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	fl <- gbl_id
		
	strip.default(which.given,which.panel,var.name,fl,...)
}

gbl_plots <- xyplot(
	gbl_g_vec ~ time_g_psi | gbl_labels, groups = gbl_labels,
	data=data.frame(
		time_g_psi = time_g_psi, gbl_g_vec = gbl_g_vec,
		gbl_labels = gbl_labels
	),
	layout = plot_gbl_dim, main="",
	strip=strip.math,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1)),
	xlab="time",
	ylab="mean function and eigenfunctions",
	as.table=TRUE,
	panel=function(x, y, subscripts, groups) {
		
		lPan <- panel.number()
		
		l <- rep(1:(L + 1), each = 1)[lPan]
		panel.grid()
		panel.superpose(
			x[order(x)], y[order(x)], subscripts, groups,
			type="l", col=data_col, pch=16, cex=0.4
		)
		
		panel.xyplot(
			time_g, vmp_res$"gbl_curves"[, l],
			col = vmp_col, type = "l", lwd = 1
		)
	}
)

print(gbl_plots)

wait()

# Compare first level scores:

norm_diff <- apply(vmp_res$zeta_1 - Reduce(rbind, zeta_1), 1, function(x) sqrt(cprod(x)))
rmse_1 <- sqrt(mean(norm_diff))

cat("The rmse for the first level scores is:", rmse_1, "\n")

wait()

# Compare second level scores:

norm_diff <- apply(
	Reduce(rbind, vmp_res$zeta_2) - Reduce(rbind, lapply(zeta_2, Reduce, f = rbind)),
	1, function(x) sqrt(cprod(x))
)
rmse_2 <- sqrt(mean(norm_diff))

cat("The rmse for the second level scores is:", rmse_2, "\n")

############ End of mlfpca_vmp.R ############