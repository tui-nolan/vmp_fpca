library(tidyverse)
library(furrr)

load(here::here("data", "US_temp_data.RData"))

plan(multisession, workers = 6)

ann_summ = function(df) {
  
  df %>% 
    group_by(m, d) %>% 
    summarize(avg_temp = mean(tmax, na.rm = TRUE), .groups = "drop") %>% 
    filter(!(m == 2 & d == 29)) %>% 
    mutate(day = row_number())
  
}

nest_weather_df = 
  weather_dat %>% 
  nest(df = -id)
  
aggregate_temp = 
  nest_weather_df %>% 
#  slice(1:100) %>% 
  mutate(
    df = future_map(df, separate, col = "date", into = c("y", "m", "d")),
    df = future_map(df, ann_summ)
  ) 
  
aggregate_temp = 
  aggregate_temp %>% 
  unnest(df)

write_csv(aggregate_temp, here::here("data", "aggregate_temp.csv"))

agg_temp_wide = 
  aggregate_temp %>% 
  select(id, avg_temp, day) %>% 
  pivot_wider(
    names_from = day, 
    values_from = avg_temp
  )

write_csv(agg_temp_wide, here::here("data", "agg_temp_wide.csv"))


aggregate_temp %>% 
  ggplot(aes(x = day, y = avg_temp, group = id)) + 
  geom_line(alpha = .2)
  