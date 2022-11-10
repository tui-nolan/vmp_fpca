######### R script: us_ml_temp_data_old.R ##########

# For gathering the US temperature data from the rnoaa package

# Created: 18 APR 2020
# Last changed: 02 SEP 2020

library(rnoaa)
library(dplyr)

load("US_temp_data.RData")

time_vec <- (1:365 - 0.5)/365

all_stations <- ghcnd_stations()
n_states <- length(state.abb)
n_year <- 1994 - 1960 + 1
for(i in 1:n_states) {
	
	state <- state.abb[i]
	
	cat("collecting data for", state, "\n")
	
	inds <- sample(which(all_stations$"state" == state))
	if(length(inds) == 0) {
		
		next
	}
	
	avail_data <- FALSE
	j <- 0
	while(!avail_data) {
		
		j <- j + 1
		if(j > length(inds)) {
			
			break
		}
		
		id_val <- all_stations$"id"[inds[j]]
		id_rows <- (weather_dat$"id" == id_val)
		
		if(!any(id_rows)) {
			
			next
		}
		
		id_data <- weather_dat[id_rows, c(2, 4)]
		
		dates <- as.character(pull(id_data, date))
		t_max <- pull(id_data, tmax)/10
		prop_na <- sum(is.na(t_max))/length(t_max)
		if(prop_na > 0.1) {
			
			inds[j] <- NA
			next
		} else {
			
			avail_data <- TRUE
		}
	}
	
	if(!avail_data) {
		
		next
	}
	
	# inds <- inds[!is.na(inds)]
	
	# if(!avail_data) {
		
		# ind <- sample(inds, 1)
		
		# id_val <- all_stations$"id"[ind]
		# id_rows <- (weather_dat$"id" == id_val)
		
		# id_data <- weather_dat[id_rows, c(2, 4)]
		
		# dates <- as.character(pull(id_data, date))
		# t_max <- pull(id_data, tmax)/10
	# }
	
	data_summary <- matrix(NA, 365, n_year)
	for(j in 1:n_year) {
		
		year <- as.character(1960 + j - 1)
		
		year_dates <- as.character(
			seq(
				as.Date(paste0(year, "-01-01")),
				as.Date(paste0(year, "-12-31")),
				by = "+1 day"
			)
		)
		
	    year_rows <- grepl(year, dates, fixed = TRUE)
		is_leap_year <- ((as.numeric(year) %% 4) == 0)
		if(is_leap_year) {
			
			feb_29 <- paste0(year, "-02-29")
			
			feb_29_ind <- which(year_dates == feb_29)
			year_dates <- year_dates[-feb_29_ind]
			
			feb_29_row <- which(dates == feb_29)
			year_rows[feb_29_row] <- FALSE
		}
		
		obs_dates <- dates[year_rows]
		
		obs_temp <- t_max[year_rows]
		# mean_temp <- mean(obs_temp, na.rm = TRUE)
		
		# below_average <- which(obs_temp < mean_temp)
		# above_average <- which(obs_temp > mean_temp)
		
		# below_diffs <- mean_temp - obs_temp[below_average]
		# obs_temp[below_average] <- mean_temp + below_diffs
		
		# above_diffs <- obs_temp[above_average] - mean_temp
		# obs_temp[above_average] <- mean_temp - above_diffs
		
		time_inds <- match(obs_dates, year_dates)
		
		data_summary[, j][time_inds] <- obs_temp
	}
	state_data <- cbind(time_vec, data_summary)
	
	col_names <- c("time", as.character(1960:1994))
	n_col <- length(col_names)
	file_name <- paste0("us_temp_data_old/", state, ".txt")
	write(col_names, file_name, ncol = n_col, append = FALSE)
	write(t(state_data), file_name, ncol = n_col, append = TRUE)
}

############ End of us_ml_temp_data_old.R ############