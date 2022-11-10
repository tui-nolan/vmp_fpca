######### R script: us_ml_temp_data.R ##########

# For gathering the US temperature data from the rnoaa package

# Created: 18 APR 2020
# Last changed: 06 SEP 2020

library(rnoaa)
library(dplyr)

load("US_temp_data.RData")

time_vec <- (1:365 - 0.5)/365
max_prop_year_mis <- 0.1

all_stations <- ghcnd_stations()
n_states <- length(state.abb)
n_year <- 1994 - 1960 + 1
for(i in 1:n_states) {
	
	state <- state.abb[i]
	
	cat("collecting data for", state, "\n")
	
	inds <- which(all_stations$"state" == state)
	inds <- inds[pull(all_stations[inds, ], first_year) == 1960]
	if(length(inds) == 0) {
		
		cat("there is no available data for", state, "\n")
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
		unique_years <- unique(as.numeric(substr(dates, 1, 4)))
		if(!all(1960:1994 %in% unique_years)) {
			
			next
		}
		
		t_max <- pull(id_data, tmax)/10
		
		avail_year <- rep(FALSE, n_year)
		for(k in 1:n_year) {
			
			year <- unique_years[k]
			year_inds <- which(as.numeric(substr(dates, 1, 4)) == year)
			t_year_data <- t_max[year_inds]
			n_obs <- length(t_year_data[!is.na(t_year_data)])
			n_mis <- max(365 - n_obs, 0)
			prop_mis <- min(n_mis/365, 1)
			
			if(prop_mis > max_prop_year_mis) {
				
				break
			} else {
				
				avail_year[k] <- TRUE
			}
		}
		avail_data <- all(avail_year)
	}
	
	if(avail_data) {
		
		cat("there is available data for", state, "\n")
	} else {
		
		cat("there is no available data for", state, "\n")
		next
	}
	
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
		
		time_inds <- match(obs_dates, year_dates)
		
		data_summary[, j][time_inds] <- obs_temp
	}
	state_data <- cbind(time_vec, data_summary)
	
	col_names <- c("time", as.character(1960:1994))
	n_col <- length(col_names)
	file_name <- paste0("us_ml_temp_data/", state, ".txt")
	write(col_names, file_name, ncol = n_col, append = FALSE)
	write(t(state_data), file_name, ncol = n_col, append = TRUE)
}