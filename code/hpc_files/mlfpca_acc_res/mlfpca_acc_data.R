######### R script: mlfpca_acc_data.R ##########

# For transforming the FPCA accuracy data from the HPC
# into the appropriate txt format

# Created: 23 AUG 2022
# Last changed: 23 AUG 2022

library(lattice)

source("../../functions/fpca_algs.r")
source("../../functions/construct_txt.r")

file_name <- "../../res/mlfpca_acc_new.txt"
construct_txt(file_name = file_name, header = TRUE)

plot_width <- 8
plot_height <- 5.6

plot_dim <- c(4, 2)              # c(ncol, nrow)

saved_location <- "../../res"

box_plot_mlfpca_sims(
	file_name, plot_dim, plot_width, plot_height,
	save_pdf = TRUE, log_acc = TRUE, saved_location = saved_location
)

score_res <- mlfpca_score_summary(file_name)

############ End of mlfpca_acc_data.R ############