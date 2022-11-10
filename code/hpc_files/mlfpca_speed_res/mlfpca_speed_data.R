######### R script: fpca_speed_data.R ##########

# For transforming the FPCA accuracy data from the HPC
# into the appropriate txt format

# Created: 16 AUG 2022
# Last changed: 16 AUG 2022

library(lattice)

source("../../functions/fpca_algs.r")
source("../../functions/construct_txt.r")

file_name <- "../../res/mlfpca_speed.txt"
construct_txt(file_name = file_name, header = TRUE)

plot_width <- 9
plot_height <- 4

plot_dim <- c(5, 1)              # c(ncol, nrow)

saved_location <- "../../res"

hist_plot_fpca_sims(
	file_name, plot_dim, plot_width, plot_height,
	fpca_mod = "standard", save_pdf = FALSE, saved_location = saved_location
)

speed_res <- fpca_speed_summary(file_name)

conv_res <- fpca_converged(file_name, 500)

############ End of fpca_acc_data.R ############