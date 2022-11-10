library(rnoaa)
library(tidyverse)

library(furrr)

# Get a list of all station IDs with data between 1960 and 1994
all_stations <- ghcnd_stations()

use_stations <-  
	all_stations %>% 
	filter(first_year < 1960, last_year > 1994) %>% 
	select(id:name) %>% 
	distinct()

plan(multisession, workers = 6)

# Pull the desired weather data for all of these stations
weather_dat <- 
	future_map(use_stations$id, ~meteo_pull_monitors(.x, date_min = "1960-01-01", date_max = "1994-12-31", var = c("PRCP", "TMAX", "TMIN"))) %>% 
	bind_rows
	

# Save the resulting data
save(use_stations, weather_dat, file = here::here("data", "US_temp_data.RData"))
