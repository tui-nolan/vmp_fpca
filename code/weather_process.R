library(tidyverse)
library(maps)
library(furrr)

load(here::here("data", "US_temp_data.RData"))

## join temperature and station information

weather_dat = 
  left_join(use_stations, weather_dat) %>% 
  select(id, latitude, longitude, state, name, date, tmax)

## recode observations coded as missing; 
## recode data entered below plausible min temps or above highest recorded temps;
## remove data from a station with a single observation
## covert to degrees celcius
## drop missing values

weather_dat = 
  weather_dat %>% 
  mutate(
    tmax = na_if(tmax, 9990),
    tmax = na_if(tmax, 990),
    tmax = replace(tmax, tmax > 540, NA),
    tmax = replace(tmax, tmax < -800, NA),
    tmax = tmax / 10) %>% 
  drop_na(date, tmax) %>% 
  filter(id != "USC00049582")

# options(future.globals.maxSize= 3 * 1024^3)

plan(multisession, workers = 6)

ann_summ = function(df) {
  
  df %>% 
    separate(date, into = c("y", "m", "d")) %>% 
    filter(!(m == "02" & d == "29")) %>% 
    group_by(m, d) %>% 
    summarize(
      n_days = sum(!is.na(tmax)),
      avg_temp = mean(tmax, na.rm = TRUE), 
      .groups = "drop") %>% 
    mutate(
      day = row_number(),
      min_n_days = min(n_days))
  
}

nest_weather_df = 
  weather_dat %>% 
  nest(df = date:tmax)

aggregate_temp = 
  nest_weather_df %>% 
  # slice(1:100) %>% 
  mutate(
    df = purrr::map(df, ann_summ)
  ) %>% 
  unnest(df) %>% 
  filter(min_n_days >= 30) %>% 
  select(id, latitude, longitude, name, day, avg_temp)

beepr::beep(sound = 5)

agg_temp_wide = 
  aggregate_temp %>% 
  select(id, avg_temp, day) %>% 
  pivot_wider(
    names_from = day, 
    values_from = avg_temp
  )

## export data

write_csv(aggregate_temp, here::here("data", "aggregate_temp.csv"))
write_csv(agg_temp_wide, here::here("data", "agg_temp_wide.csv"))


## export two exploratory graphs
ggp_temp_curves = 
  aggregate_temp %>% 
  ggplot(aes(x = day, y = avg_temp, group = id)) + 
  geom_line(alpha = .2) +
  theme_minimal()


station_locations = 
  aggregate_temp %>% 
  select(id, latitude, longitude) %>% 
  distinct()

world_map = map_data("world")

ggp_station_map = 
  ggplot(world_map, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "lightgray", colour = "white") +
  geom_point(data = station_locations, aes(x = longitude, y = latitude, group = NULL))

ggsave(here::here("images", "ggp_curves.png"), ggp_temp_curves, width = 8, height = 6)
ggsave(here::here("images", "ggp_map.png"), ggp_station_map, width = 8, height = 6)
