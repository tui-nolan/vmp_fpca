######### R script: fpca_us_weather_vmp.R ##########

# For performing Bayesian FPCA via VMP
# for the US weather dataset.

# Created: 14 JAN 2020
# Last changed: 03 MAR 2020

# Load libraries:

library(MASS)
library(magic)
library(lattice)
library(pracma)
library(ellipse)
library(fda)

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

# Gather the data:

n_sample <- 4                             # number of curves for the plots
N_sample <- rep(NA, n_sample)             # specific curves for the plots

us_weather_data <- read.table("us_weather_data.txt", header=TRUE)

states <- unique(us_weather_data$state)
n_states <- length(states)
sample_states <- sample(states, n_sample, replace = FALSE)
sample_states <- names(sort(sapply(sample_states, function(x) which(states == x))))
sample_rows <- vector("list", length=n_states)
for(i in 1:n_states) {
	
	state_i <- states[i]
	state_inds <- which(us_weather_data$"state" == state_i)
	n_state_sample <- length(state_inds)
	
	if(n_state_sample > 2) {
		
		sample_rows[[i]] <- sample(state_inds, 3, replace = FALSE)
	} else {
		
		sample_rows[[i]] <- state_inds
	}
	
	if(state_i %in% sample_states) {
		
		ind <- which(sample_states == state_i)
		N_sample[ind] <- sample(sample_rows[[i]], 1)
	}
}
sample_rows <- sort(Reduce(c, sample_rows))
N_sample <- sapply(N_sample, function(x) which(sample_rows == x))

dataset <- t(as.matrix(us_weather_data[sample_rows, -c(1, 367)]))

N <- dim(dataset)[2]                      # number of curves
T_vec <- rep(nrow(dataset), N)            # number of time observations for each curve
n_int_knots <- 10                         # number of interior knots
K <- n_int_knots + 2                      # number of spline basis functions
L <- 4                                    # guess of the number of FPCA basis functions
data_col <- "black"                       # colour of the data in the plots
criterion <- 1e-10                        # convergence criterion
d <- (K+2)*(L+1)                          # dimension of spline vector

n_vmp <- 250                              # number of VMP iterations
n_g <- 1000                               # length of the plotting grid
vmp_col <- "red"                          # colour of the VMP plots
vmp_lwd <- 1                              # line width for vmp plots

# Set up plot-grid dimensions:

L_present <- 2                            # number of basis functions to present

plot_dim <- c(2, 2)                       # (ncol, nrow) for curve plots
plot_gbl_dim <- c(4, 1)                   # (ncol, nrow) for basis function plots

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

# Gather a sample of the response data:

Y <- vector("list", length=N)
time_obs <- vector("list", length=N)
for(i in 1:N) {
	
	time_obs[[i]] <- sort(sample(1:365, T_vec[i]))
	Y[[i]] <- unname(dataset[,i][time_obs[[i]]])
	time_obs[[i]] <- (time_obs[[i]] - 0.5)/365
}

Y_vec <- Reduce(c, Y)

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

# Set up plotting grid

time_g <- seq(0, 1, length.out=n_g)

X_g <- X_design(time_g)
Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
C_g <- cbind(X_g, Z_g)

# Set up sampled points for FPC plots:

n_points <- 20
n_g_sample <- n_g/n_points
time_g_sample_points <- seq((n_g_sample/2), n_g, by=n_g_sample)
time_g_sample <- time_g[time_g_sample_points]

# Plot the data:

if(!print_pdf) {
	
	plot_fpca_data(N_sample, time_obs, Y, plot_dim, data_col)
	
	wait()
}

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

eta_nu <- list(eta_in[[1]], eta_in[[2]])
q_nu <- gauss_q(eta_nu, use_vech=FALSE)
mu_q_nu <- q_nu[[1]]
Sigma_q_nu <- q_nu[[2]]

mu_inds <- 1:(K+2)
psi_inds <- (1:d)[-mu_inds]
psi_groups <- ceiling(psi_inds/(K+2))-1
psi_inds <- split(psi_inds, psi_groups)

mu_q_nu_mu <- mu_q_nu[mu_inds]
Sigma_q_nu_mu <- Sigma_q_nu[mu_inds, mu_inds]
mu_q_mu <- as.vector(C_g%*%mu_q_nu_mu)

mu_q_nu_psi <- vector("list", length=L)
Sigma_q_nu_psi <- vector("list", length=L)
for(l in 1:L) {
	
	psi_l_inds <- psi_inds[[l]]
	mu_q_nu_psi[[l]] <- mu_q_nu[psi_l_inds]
	Sigma_q_nu_psi[[l]] <- Sigma_q_nu[psi_l_inds, psi_l_inds]
}

M_q_V_psi <- Reduce(cbind, mu_q_nu_psi)
M_q_Psi <- C_g%*%M_q_V_psi

mu_q_zeta <- vector("list", length=N)
Sigma_q_zeta <- vector("list", length=N)
for(i in 1:N) {
	
	eta_zeta <- list(eta_in[[3]][,i], eta_in[[4]][,i])
	q_zeta <- gauss_q(eta_zeta, use_vech=TRUE)
	mu_q_zeta[[i]] <- q_zeta[[1]]
	Sigma_q_zeta[[i]] <- q_zeta[[2]]
}

M_q_Zeta <- Reduce(rbind, mu_q_zeta)

one_N <- rep(1, N)
mu_mat <- tcrossprod(mu_q_mu, one_N)
Y_mat <- mu_mat + tcrossprod(M_q_Psi, M_q_Zeta)

M_q_Psi_svd <- svd(M_q_Psi)
U_orth <- M_q_Psi_svd$u
D_diag <- diag(M_q_Psi_svd$d)
V_orth <- M_q_Psi_svd$v

M_q_Zeta_rotn <- M_q_Zeta%*%V_orth%*%D_diag
mu_Zeta_rotn <- apply(M_q_Zeta_rotn, 2, mean)

mu_q_mu <- mu_q_mu + U_orth%*%mu_Zeta_rotn
M_q_Zeta_shift <- M_q_Zeta_rotn - tcrossprod(one_N, mu_Zeta_rotn)

eigen_M_q_Zeta_shift <- eigen(cov(M_q_Zeta_shift))
Q <- eigen_M_q_Zeta_shift$vectors
Lambda <- diag(eigen_M_q_Zeta_shift$values + 1e-10)
Lambda_inv <- diag(1/(eigen_M_q_Zeta_shift$values + 1e-10))
S <- Q%*%sqrt(Lambda)
S_inv <- tcrossprod(sqrt(Lambda_inv), Q)

Psi_hat <- U_orth%*%S
Zeta_hat <- tcrossprod(M_q_Zeta_shift, S_inv)

norm_const <- rep(NA, L)
for(l in 1:L) {
	
	norm_const[l] <- sqrt(trapint(time_g, (Psi_hat[,l])^2))
	if(norm_const[l]!=0) {
		
		Psi_hat[,l] <- Psi_hat[,l]/norm_const[l]
		Zeta_hat[,l] <- norm_const[l]*Zeta_hat[,l]
	}
}

mu_q_zeta <- split(Zeta_hat, row(Zeta_hat))

scale_mat <- diag(norm_const)
for(i in 1:N) {
	
	mat_transform <- S_inv%*%tcrossprod(D_diag, V_orth)
	Sigma_q_zeta[[i]] <- tcrossprod(mat_transform%*%Sigma_q_zeta[[i]], mat_transform)
	Sigma_q_zeta[[i]] <- tcrossprod(scale_mat%*%Sigma_q_zeta[[i]], scale_mat)
}

Y_summary <- vector("list", length=N)
for(i in 1:N) {
	
	sd_vec <- sqrt(diag(tcrossprod(Psi_hat%*%Sigma_q_zeta[[i]], Psi_hat)))
	
	Y_summary[[i]] <- matrix(NA, nrow=n_g, ncol=3)
	Y_summary[[i]][,1] <- Y_mat[,i] + qnorm(0.025)*sd_vec
	Y_summary[[i]][,2] <- Y_mat[,i]
	Y_summary[[i]][,3] <- Y_mat[,i] + qnorm(0.975)*sd_vec
}

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

gbl_summary <- cbind(mu_q_mu, Psi_hat)
gbl_summary[,2] <- -gbl_summary[,2]
function_names <- c("mu", rep("NA", L_present))
for(l in 1:L_present) {
	
	function_names[l+1] <- paste("psi_", l, sep="")
}
write(function_names, "./res/us_gbl_funcs.txt", ncol = L_present + 1, append = FALSE)
write(gbl_summary, "./res/us_gbl_funcs.txt", ncol = L_present + 1, append=TRUE)

zeta_hat_results <- t(Reduce(cbind, mu_q_zeta[N_sample]))
zeta_hat_results[,1] <- -zeta_hat_results[,1]
write(sample_states, "./res/us_zeta_post.txt", ncol = n_sample, append = FALSE)
write(zeta_hat_results, "./res/us_zeta_post.txt", ncol = n_sample, append=TRUE)

post_fit <- do.call(cbind, lapply(Y_summary, `[`, ,2))[,N_sample]
write(sample_states, "./res/us_fits.txt", ncol = n_sample, append = FALSE)
write(post_fit, "./res/us_fits.txt", ncol = n_sample, append=TRUE)

# Plot the fitted curves:

Y_sample <- Y[N_sample]
time_obs_sample <- time_obs[N_sample]
T_sample <- sapply(Y_sample, length)
length_N <- length(N_sample)

Y_vec <- Reduce(c, Y_sample)
time_vec <- Reduce(c, time_obs_sample)

curve_id <- sample_states
curve_labels <- vector("list", length=length_N)
for(i in 1:length_N) {
	
	curve_val <- curve_id[i]
	curve_labels[[i]] <- rep(curve_val, T_sample[i])
}
curve_labels <- do.call(c, curve_labels)
curve_labels <- factor(curve_labels, levels=curve_id)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	fl <- curve_id
		
	strip.default(which.given,which.panel,var.name,fl,...)
}

if(print_pdf) {
	
	pdf(
		"./res/us_fits.pdf",
		width=plot_width, height=plot_height
	)
}

fitted_data_plots <- xyplot(
	Y_vec ~ time_vec | curve_labels, groups=curve_labels,
	data=data.frame(
		time_vec=time_vec, Y_vec=Y_vec,
		curve_labels=curve_labels
	),
	layout=plot_dim, main="", type = c("p", "g"),
	strip=strip.math,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1.2)),
	xlab="time (years)",
	ylab="temperature (\u00B0C)",
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
			time_g, Y_summary[[i]][,2],
			col=vmp_col, type="l", lwd=2
		)
#		panel.xyplot(
#			time_g, Y_summary[[i]][,1],
#			col=vmp_col, type="l", lwd=1, lty=2
#		)
#		panel.xyplot(
#			time_g, Y_summary[[i]][,3],
#			col=vmp_col, type="l", lwd=1, lty=2
#		)
	}
)

print(fitted_data_plots)

if(print_pdf) {
	
	dev.off()
} else {
	
	wait()
}

# Construct plots of the global functions:

var_cont <- round(norm_const^2/sum(norm_const^2), 3)*100
var_cont <- var_cont[1:L_present]

time_g_gbl <- rep(time_g, L_present)

gbl_g_vec <- as.vector(gbl_summary[,2:(L_present+1)])
gbl_labels <- vector("list", length=L_present)
gbl_id <- rep(NA, L_present)

mu_g_gbl <- rep(mu_q_mu, L_present)
mu_labels <- vector("list", length=L_present)
mu_id <- rep(NA, L_present)

mean_mat <- matrix(rep(mu_q_mu[time_g_sample_points,], L_present), ncol=L_present)
fpc_bf_mat <- gbl_summary[time_g_sample_points, 2:(L_present+1)]
scale_mat <- diag(c(3, 1.2))
mean_plus_fpc_mat <- mean_mat + fpc_bf_mat%*%scale_mat
mean_minus_fpc_mat <- mean_mat - fpc_bf_mat%*%scale_mat

for(l in 1:L_present) {
	
	gbl_id[l] <- parse(text=paste0("psi[", l, "] (t) : ", var_cont[l], "*\'%\'"))
	gbl_labels[[l]] <- rep(gbl_id[l], n_g)
	
	mu_id[l] <- parse(text=paste0("mu (t) %+-% delta[", l, "] *\ psi[", l, "] (t)"))
	mu_labels[[l]] <- rep(mu_id[l], n_g)
}

gbl_plot_id <- c(gbl_id, mu_id)
gbl_plot_labels <- c(do.call(c, gbl_labels), do.call(c, mu_labels))
gbl_plot_labels <- factor(gbl_plot_labels, levels = gbl_plot_id)

gbl_plot_vec <- c(gbl_g_vec, mu_g_gbl)
time_plot_vec <- rep(time_g_gbl, 2)
mean_pos_shift_mat <- cbind(matrix(NA, n_points, L_present), mean_plus_fpc_mat)
mean_neg_shift_mat <- cbind(matrix(NA, n_points, L_present), mean_minus_fpc_mat)
write(
	c("time", function_names[2:(L_present + 1)]),
	"./res/us_pos_shift.txt", ncol = L_present + 1, append = FALSE
)
write(
	cbind(time_g_sample, mean_plus_fpc_mat),
	"./res/us_pos_shift.txt", ncol = L_present + 1, append=TRUE
)
write(
	c("time", function_names[2:(L_present + 1)]),
	"./res/us_neg_shift.txt", ncol = L_present + 1, append = FALSE
)
write(
	cbind(time_g_sample, mean_minus_fpc_mat),
	"./res/us_neg_shift.txt", ncol = L_present + 1, append=TRUE
)

y_lims <- rep(list(c(-2, 2), c(0, 32)), each = L_present)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	fl <- gbl_plot_id
		
	strip.default(which.given,which.panel,var.name,fl,...)
}

if(print_pdf) {
	
	pdf("./res/us_bf.pdf",width=plot_width, height=plot_height)
}

basis_plots <- xyplot(
	gbl_plot_vec ~ time_plot_vec | gbl_plot_labels, groups = gbl_plot_labels,
	data=data.frame(
		time_plot_vec = time_plot_vec, gbl_plot_vec = gbl_plot_vec,
		gbl_plot_labels = gbl_plot_labels
	),
	layout=c(L_present, 2), main="",
	strip=strip.math, col=vmp_col, type=c("l", "g"), lwd=2,
	scales = list(x = list(tick.number = 3), y = list(relation = "free", limits = y_lims)),
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1.2)),
	xlab="time (years)",
	ylab="temperature (\u00B0C)",
	as.table=TRUE,
	panel=function(x, y, subscripts, groups) {
		
		lPan <- panel.number()
		l <- rep((1:(2*L_present)), each=1)[lPan]
		panel.grid()
		panel.superpose(
			x[order(x)], y[order(x)], subscripts, groups,
			type="l", col=vmp_col, lwd=2
		)
		panel.xyplot(
			time_g_sample, mean_pos_shift_mat[,l],
			col="black", pch="+", cex=1
		)
		panel.xyplot(
			time_g_sample, mean_neg_shift_mat[,l],
			col="black", pch="-", cex=1
		)
	}
)

print(basis_plots)

if(print_pdf) {
	
	dev.off()
}

############ End of fpca_us_weather_vmp.R ############