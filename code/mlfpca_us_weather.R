######### R script: mlfpca_us_weather.R ##########

# For applying mlfpca to the US weather dataset

# Created: 30 AUG 2022
# Last changed: 05 SEP 2022

# Load libraries:

library(MASS)
library(magic)
library(lattice)
library(pracma)
library(data.table)
library(readtext)

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

files <- list.files(
	"./us_temp_data", recursive = TRUE,
	pattern = "\\.txt$", full.names = TRUE
)
us_weather_data <- lapply(files, function(x) read.table(x, header = TRUE))
n_states <- length(us_weather_data)
for(i in 1:n_states) {
	
	na_years <- apply(us_weather_data[[i]], 2, function(x) sum(!is.na(x)) < 5)
	if(any(na_years)) {
		
		na_cols <- which(na_years)
		us_weather_data[[i]] <- us_weather_data[[i]][, -na_cols]
	}
}

# Establish simulation variables:

states <- substr(files, 16, 17)
names(us_weather_data) <- states
N <- length(us_weather_data)
M <- sapply(us_weather_data, ncol) - 1
T_vec <- vector("list", length = N)
time_obs <- vector("list", length = N)
Y <- vector("list", length = N)
unique_time_obs <- vector("list", length = N)
for(i in 1:N) {
	
	time_vec <- us_weather_data[[i]][, 1]
	all_obs <- us_weather_data[[i]][, -1]
	
	T_vec[[i]] <- rep(NA, M[i])
	time_obs[[i]] <- vector("list", length = M[i])
	Y[[i]] <- vector("list", length = M[i])
	for(j in 1:M[i]) {
		
		temp_vec <- all_obs[, j]
		na_inds <- is.na(temp_vec)
		
		time_obs[[i]][[j]] <- time_vec[!na_inds]
		Y[[i]][[j]] <- temp_vec[!na_inds]
		T_vec[[i]][j] <- length(time_obs[[i]][[j]])
	}
	
	unique_time_obs[[i]] <- sort(unique(Reduce(c, time_obs[[i]])))
}
unique_time_obs <- sort(Reduce(c, unique_time_obs))

m_sample <- 5                             # number of second level curves for the plots
M_sample <- floor(seq(1, max(M), length.out = m_sample))# specific second level curves for the plots

main_states <- c("CA", "NY", "TX", "WA", "FL")
N_sample <- sort(match(main_states, states)) # specific curves for the plots
n_sample <- length(N_sample)

criterion  <- 1e-6                        # convergence criterion

n_int_knots <- 10                         # number of interior knots
K <- n_int_knots + 2                      # number of spline basis functions
L_1 <- 4                                  # number of first level eigenfunctions
L_2 <- 10:20                                  # number of second level eigenfunctions
L <- L_1 + L_2
data_col <- "grey51"                      # colour of the data in the plots
delta <- 1e-3

n_vmp <- 200                              # number of VMP iterations
n_g <- 1000                               # length of the plotting grid
vmp_col <- "red"                          # colour of the VMP plots

# Set up plot-grid dimensions:

L_1_plots <- 2                            # number of level 1 eigenfunctions for plots
L_2_plots <- 4                            # number of level 2 eigenfunctions for plots

plot_dim <- c(m_sample, n_sample)         # (ncol, nrow) for curve plots
plot_psi_1_dim <- c(L_1_plots, 1)               # (ncol, nrow) for basis function plots
plot_psi_2_dim <- c(L_2_plots, 1)               # (ncol, nrow) for basis function plots
plot_score_1_dim <- c(2, 2)
plot_score_2_dim <- c(m_sample, n_sample)

plot_width <- 7
plot_height <- 7

print_pdf <- FALSE              # save the plots in a PDF?

colour <- c("indianred", "dodgerblue", "magenta", "sandybrown", "seagreen", "springgreen4")

# Set up fixed parameters:

int_knots <- quantile(unique_time_obs, seq(0,1,length=K)[-c(1,K)])

X <- vector("list", length=N)
Z <- vector("list", length=N)
C <- vector("list", length=N)
for(i in 1:N) {
	
	X[[i]] <- vector("list", length = M[i])
	Z[[i]] <- vector("list", length = M[i])
	C[[i]] <- vector("list", length = M[i])
	
	for(j in 1:M[i]) {
		
		X[[i]][[j]] <- X_design(time_obs[[i]][[j]])
		Z[[i]][[j]] <- ZOSull(time_obs[[i]][[j]], range.x=c(0, 1), intKnots=int_knots)
		C[[i]][[j]] <- cbind(X[[i]][[j]], Z[[i]][[j]])

	}
}

# Plot the data:

if(print_pdf) {
	
	pdf(
		"/Users/tuinolan/Desktop/fulbright/images/us_ml_data.pdf",
		width=plot_width, height=plot_height
	)
}

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
	
	curve_labels_1[[i]] <- rep.int(M_sample + 1960 - 1, times =  T_sample[[i]])
}
curve_labels_1 <- Reduce(c, curve_labels_1)
curve_id_1 <- rep(NA, m_sample)
for(j in 1:m_sample) {
	
	curve_id_1[j] <- parse(text = as.character(M_sample[j] + 1960 - 1))
}
curve_labels_1 <- factor(curve_labels_1, levels=curve_id_1)

curve_labels_2 <- vector("list", length = n_sample)
curve_id_2 <- rep(NA, n_sample)
for(i in 1:n_sample) {
	
	N_i <- N_sample[i]
	curve_id_2[i] <- parse(text = states[N_i])
	curve_labels_2[[i]] <- rep(states[N_i], sum(T_sample[[i]]))
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

raw_data_plots <- xyplot(
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
		
		panel.grid() 
		panel.superpose(
			x[order(x)], y[order(x)], subscripts, groups,
			type="p", col = data_col, pch = 16, cex = 0.4
		)
	}
)

print(raw_data_plots)

if(print_pdf) {
	
	dev.off()
}

wait()

# Set up plotting grid

time_g <- seq(0,1, length.out=n_g)

X_g <- X_design(time_g)
Z_g <- ZOSull(time_g, range.x=c(0, 1), intKnots=int_knots)
C_g <- cbind(X_g, Z_g)

# Establish hyperparameters:

sigsq_beta <- 1e10
Sigma_beta <- sigsq_beta*diag(2)
sigma_zeta <- 1
sigsq_zeta <- sigma_zeta^2
mu_beta <- rep(0, 2)
A <- 1e5

# VMP algorithm:

cum_var_2 <- matrix(1, length(L_2), max(L_2))
for(l in 1:length(L_2)) {
	
	vmp_res <- vmp_mlfpca(
		n_vmp = 50, N, M, L_1, L_2[l], C, Y,
		sigsq_beta, A, criterion, plot_elbo = TRUE
	)
	
	eta_vec <- vmp_res$"eta_vec"
	
	eta_in <- list(
		eta_vec$"p(nu|Sigma_nu)->nu",
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
		eta_vec$"p(zeta)->zeta",
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
	)
	
	fpca_res <- mlfpc_rotation(
		eta_in, time_g, C_g, L_1, L_2[l],
		N_sample, M_sample, Psi_g = NULL
	)
	
	Zeta_2 <- fpca_res$zeta_2
	
	eig_vals_2 <- apply(Reduce(rbind, Zeta_2), 2, var)
	exp_var <- eig_vals_2/sum(eig_vals_2)
	cum_var_vec <- cumsum(exp_var)
	cum_var_2[l, 1:L_2[l]] <- cum_var_vec
}

L_2 <- match(length(L_2), round(apply(cum_var_2, 2, sum), 4))

L_1 <- 2:10

cum_var_1 <- matrix(1, length(L_1), max(L_1))
for(l in 1:length(L_1)) {
	
	vmp_res <- vmp_mlfpca(
		n_vmp = 50, N, M, L_1[l], L_2, C, Y,
		sigsq_beta, A, criterion, plot_elbo = TRUE
	)
	
	eta_vec <- vmp_res$"eta_vec"
	
	eta_in <- list(
		eta_vec$"p(nu|Sigma_nu)->nu",
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
		eta_vec$"p(zeta)->zeta",
		eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
	)
	
	fpca_res <- mlfpc_rotation(
		eta_in, time_g, C_g, L_1[l], L_2,
		N_sample, M_sample, Psi_g = NULL
	)
	
	Zeta_1 <- fpca_res$zeta_1
	
	eig_vals_1 <- apply(Zeta_1, 2, var)
	exp_var <- eig_vals_1/sum(eig_vals_1)
	cum_var_vec <- cumsum(exp_var)
	cum_var_1[l, 1:L_1[l]] <- cum_var_vec
}

L_1 <- match(length(L_1), round(apply(cum_var_1, 2, sum), 4))

# Run final VMP algorithm for US temperature data:

vmp_alg <- vmp_mlfpca(
	n_vmp, N, M, L_1, L_2, C, Y,
	sigsq_beta, A, criterion, plot_elbo = TRUE
)

eta_vec <- vmp_alg$"eta_vec"

# Get the posterior estimates

eta_in <- list(
	eta_vec$"p(nu|Sigma_nu)->nu",
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->nu",
	eta_vec$"p(zeta)->zeta",
	eta_vec$"p(Y|nu,zeta,sigsq_eps)->zeta"
)

vmp_res <- mlfpc_rotation(
	eta_in, time_g, C_g, L_1, L_2,
	N_sample, M_sample, Psi_g = NULL
)

Psi_1_hat <- vmp_res$"gbl_curves"[, 1:L_1 + 1]
Psi_2_hat <- vmp_res$"gbl_curves"[, 1:L_2 + L_1 + 1]
Zeta_1_hat <- vmp_res$"zeta_1"
Zeta_2_hat <- vmp_res$"zeta_2"
Cov_zeta_1_hat <- vmp_res$"Cov_zeta_1"
Cov_zeta_2_hat <- vmp_res$"Cov_zeta_2"

gamma_1_hat <- apply(Zeta_1_hat, 2, var)
var_cont_1 <- gamma_1_hat/sum(gamma_1_hat)*100

gamma_2_hat <- apply(Reduce(rbind, Zeta_2_hat), 2, var)
var_cont_2 <- gamma_2_hat/sum(gamma_2_hat)*100

# Plot the first level FPC eigenfunctions:

sign_Psi_1 <- sign(Psi_1_hat[, 1][1])
Psi_1_hat[, 1] <- sign_Psi_1*Psi_1_hat[, 1]
Zeta_1_hat[, 1] <- sign_Psi_1*Zeta_1_hat[, 1]

sign_Psi_2 <- sign(Psi_1_hat[, 2][1])
Psi_1_hat[, 2] <- sign_Psi_2*Psi_1_hat[, 2]
Zeta_1_hat[, 2] <- sign_Psi_2*Zeta_1_hat[, 2]

time_g_psi_1 <- rep(time_g, L_1_plots)
psi_1_vec <- as.vector(Psi_1_hat[, 1:L_1_plots])
psi_1_labels <- vector("list", length = L_1_plots)
psi_1_id <- rep(NA, L_1_plots)
for(l in 1:L_1_plots) {
	
	psi_1_id[l] <- parse(text=paste0("psi[", l, "]^(1)~(t) : ", round(var_cont_1[l]), "*\'%\'"))
	psi_1_labels[[l]] <- rep(psi_1_id[l], n_g)
}
psi_1_labels <- do.call(c, psi_1_labels)
psi_1_labels <- factor(psi_1_labels, levels = psi_1_id)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	fl <- psi_1_id
		
	strip.default(which.given,which.panel,var.name,fl,...)
}

if(print_pdf) {
	
	pdf(
		"res/us_ml_bf_1.pdf",
		width = 3, height = 4
	)
}

psi_1_plots <- xyplot(
	psi_1_vec ~ time_g_psi_1 | psi_1_labels, groups = psi_1_labels,
	data=data.frame(
		time_g_psi_1 = time_g_psi_1, psi_1_vec = psi_1_vec,
		psi_1_labels = psi_1_labels
	),
	layout = plot_psi_1_dim, main="",
	strip=strip.math, col = vmp_col, type="l", lwd = 2,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1.4)),
	xlab="normalized time",
	ylab="level one FPC eigenfunctions",
	as.table=TRUE
)

print(psi_1_plots)

if(print_pdf) {
	
	dev.off()
}

wait()

# Plot the second level FPC eigenfunctions:

sign_Psi_1 <- sign(Psi_2_hat[, 1][1])
diag_mat <- diag(c(sign_Psi_1, rep(1, L_2 - 1)))
Psi_2_hat[, 1] <- sign_Psi_1*Psi_2_hat[, 1]
Zeta_2_hat <- lapply(Zeta_2_hat, function(X) X %*% diag_mat)

sign_Psi_2 <- sign(Psi_2_hat[, 2][1])
diag_mat <- diag(c(1, sign_Psi_2, rep(1, L_2 - 2)))
Psi_2_hat[, 2] <- sign_Psi_2*Psi_2_hat[, 2]
Zeta_2_hat <- lapply(Zeta_2_hat, function(X) X %*% diag_mat)

time_g_psi_2 <- rep(time_g, L_2_plots)
psi_2_vec <- as.vector(Psi_2_hat[, 1:L_2_plots])
psi_2_labels <- vector("list", length = L_2_plots)
psi_2_id <- rep(NA, L_2_plots)
for(l in 1:L_2_plots) {
	
	psi_2_id[l] <- parse(text=paste0("psi[", l, "]^(2)~(t) : ", round(var_cont_2[l]), "*\'%\'"))
	psi_2_labels[[l]] <- rep(psi_2_id[l], n_g)
}
psi_2_labels <- do.call(c, psi_2_labels)
psi_2_labels <- factor(psi_2_labels, levels = psi_2_id)

strip.math <- function(
	which.given, which.panel, var.name, factor.levels, ...
) {
	
	fl <- psi_2_id
		
	strip.default(which.given,which.panel,var.name,fl,...)
}

if(print_pdf) {
	
	pdf(
		"res/us_ml_bf_2.pdf",
		width=6, height=4
	)
}

psi_2_plots <- xyplot(
	psi_2_vec ~ time_g_psi_2 | psi_2_labels, groups = psi_2_labels,
	data=data.frame(
		time_g_psi_2 = time_g_psi_2, psi_2_vec = psi_2_vec,
		psi_2_labels = psi_2_labels
	),
	layout = plot_psi_2_dim, main="",
	strip=strip.math, col = vmp_col, type="l", lwd = 2,
	par.strip.text=list(cex=0.8),
	par.settings = list(layout.heights = list(strip = 1.4)),
	xlab="normalized time",
	ylab="level two FPC eigenfunctions",
	as.table=TRUE
)

print(psi_2_plots)

if(print_pdf) {
	
	dev.off()
}

wait()

# Plot level 1 scores:

if(print_pdf) {
	
	pdf(
		"res/us_scores_1.pdf",
		width=4, height=4
	)
}

x_range <- range(Zeta_1_hat[N_sample, 1])
x_lim_1 <- x_range[1] - 0.1*diff(x_range)
x_lim_2 <- x_range[2] + 0.1*diff(x_range)
x_lim <- c(x_lim_1, x_lim_2)

y_range <- range(Zeta_1_hat[N_sample, 2])
y_lim_1 <- y_range[1] - 0.1*diff(y_range)
y_lim_2 <- y_range[2] + 0.1*diff(y_range)
y_lim <- c(y_lim_1, y_lim_2)

plot(
	Zeta_1_hat[N_sample, 1:2], pch = 16, cex = 0.4, col = vmp_col,
	xlab = "first eigenfunction's score",
	ylab = "second eigenfunction's score",
	xlim = c(-7, 10), ylim = y_lim
)
abline(h = 0, lwd = 1, col = "gray")
abline(v = 0, lwd = 1, col = "gray")
for(i in 1:n_sample) {
	
	N_i <- N_sample[i]
	cred_bound <- ellipse(Cov_zeta_1_hat[[N_i]][1:2, 1:2], centre = Zeta_1_hat[N_i, 1:2])
	lines(cred_bound, col = vmp_col)
}
text(Zeta_1_hat[N_sample, ], labels = states[N_sample], pos = 1)

if(print_pdf) {
	
	dev.off()
}

wait()

# Plot level 2 scores:

if(print_pdf) {
	
	pdf(
		"/res/us_scores_2.pdf",
		width=plot_width, height=plot_height
	)
}

Zeta_2_sample <- lapply(Zeta_2_hat[N_sample], function(X) X[M_sample, 1:2])

x_range <- range(Reduce(rbind, Zeta_2_sample)[, 1])
x_lim_1 <- x_range[1] - 0.1*diff(x_range)
x_lim_2 <- x_range[2] + 0.1*diff(x_range)
x_lim <- c(x_lim_1, x_lim_2)

y_range <- range(Reduce(rbind, Zeta_2_sample)[, 2])
y_lim_1 <- y_range[1] - 0.1*diff(y_range)
y_lim_2 <- y_range[2] + 0.1*diff(y_range)
y_lim <- c(y_lim_1, y_lim_2)

plot(
	Zeta_2_sample[[1]],
	xlab = "first eigenfunction's score",
	ylab = "second eigenfunction's score",
	xlim = x_lim, ylim = y_lim, type = "l"
)
points(Zeta_2_sample[[1]], pch = 16, cex = 0.8, col = colour[1:m_sample])
for(i in 2:n_sample) {
	
	lines(Zeta_2_sample[[i]])
	points(Zeta_2_sample[[i]], pch = 16, cex = 0.8, col = colour[1:m_sample])
}
abline(h = 0, lwd = 1, col = "gray")
abline(v = 0, lwd = 1, col = "gray")
text(Zeta_1_hat[sample_1, ], labels = states[sample_1], pos = 1)










zeta_21_hat <- lapply(Zeta_2_hat[N_sample], function(X) X[M_sample, 1])
years <- (1960:1994)[M_sample]

y_range <- range(unlist(zeta_21_hat))
y_lim_1 <- y_range[1] - 0.1*diff(y_range)
y_lim_2 <- y_range[2] + 0.1*diff(y_range)
y_lim <- c(y_lim_1, y_lim_2)

plot(
	years, zeta_21_hat[[1]], type = "l", col = colour[1],
	xlab = "year", ylab = "score", ylim = y_lim
)
for(j in 2:n_sample) {
	
	lines(years, zeta_21_hat[[j]], col = colour[j])
}

legend(
	"topleft", legend = states[N_sample],
	col = colour[1:n_sample], lty = 1
)

if(print_pdf) {
	
	dev.off()
}

############ End of mlfpca_us_weather.R ############