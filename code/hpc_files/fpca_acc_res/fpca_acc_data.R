######### R script: fpca_acc_data.R ##########

# For transforming the FPCA accuracy data from the HPC
# into the appropriate txt format

# Created: 16 AUG 2022
# Last changed: 16 AUG 2022

library(lattice)

source("../../functions/fpca_algs.r")
source("../../functions/construct_txt.r")

file_name <- "../../res/fpca_acc.txt"
construct_txt(file_name = file_name, header = TRUE)

plot_width <- 9
plot_height <- 4

plot_dim <- c(5, 1)              # c(ncol, nrow)

saved_location <- "../../res"

box_plot_fpca_sims(
	file_name, plot_dim, plot_width, plot_height,
	save_pdf = TRUE, log_acc = TRUE, saved_location = saved_location
)

score_res <- fpca_score_summary(file_name)

############ End of fpca_acc_data.R ############