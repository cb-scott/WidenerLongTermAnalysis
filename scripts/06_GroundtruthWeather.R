### Read in data from our weather station
library(readxl)
library(tidyverse)
#list the weather files 
weather.widener = list.files("data/Widener Sensor Data/", full.names = T)
inwide = read_excel(weather.widener[1], skip = 5, na = "--")
all_weather = data.frame()
for(file in weather.widener){
  tmp = read_excel(file, skip = 5, na = "--")
  all_weather = bind_rows(all_weather, tmp)
}

colnames(all_weather) = janitor::make_clean_names(colnames(all_weather))
all_weather = all_weather %>% select(date_time, temp_c, hum_percent,avg_wind_speed_m_s, soil_moisture_1_cb, high_wind_speed_m_s,  high_rain_rate_mm)
all_weather = all_weather %>% mutate(date = date(date_time)) %>% group_by(date) %>%
  summarise(across(temp_c:soil_moisture_1_cb, mean, na.rm = T),
            total_precip_mm = mean(high_rain_rate_mm, na.rm = T), 
            max_wind_speed_m_s = max(high_wind_speed_m_s, na.rm =T))
colnames(all_weather)=paste("WIDE_", colnames(all_weather), sep = "")


weather.station = read_xlsx("data/KT4CQ4Z5_1.xlsx", skip = 12)
colnames(weather.station) = paste("DURH_", janitor::make_clean_names(colnames(weather.station)), sep = "") 
DURH2WIDE= weather.station %>%
  mutate(DURH_date = as.Date(DURH_date, format = "%Y-%m-%d")) %>%
  left_join(all_weather, by = c("DURH_date" = "WIDE_date")) %>%
  mutate(across(!DURH_date, as.numeric))

get_corr <- function(vec1, vec2){
  tmp = data.frame(cbind(vec1, vec2)) %>% drop_na()
  x =cor(tmp)
  return(x)
}
library(ggplot2)
library(ggpubr)
library(patchwork)
p1 <- DURH2WIDE %>% ggplot(aes(x = DURH_total_precipitation_mm, y = WIDE_total_precip_mm)) +
  geom_point(color = "#3a86d4", size = 2.5, alpha = 0.7) +  # Plot data points
  geom_smooth(method = "lm", color = "red", fill = "lightgray", alpha = 0.4) + # Trend line + 95% CI
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 2) + # Auto-calculated correlation label
  theme_minimal() + 
  theme(text = element_text(size = 8))

p2 <- DURH2WIDE %>% ggplot(aes(x = DURH_average_air_temperature_c, y = WIDE_temp_c)) +
  geom_point(color = "#3a86d4", size = 2.5, alpha = 0.7) +  # Plot data points
  geom_smooth(method = "lm", color = "red", fill = "lightgray", alpha = 0.4) + # Trend line + 95% CI
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 2) + # Auto-calculated correlation label
  theme_minimal() + 
  theme(text = element_text(size = 8))

p3 <- DURH2WIDE %>% ggplot(aes(x = DURH_average_relative_humidity_percent, y = WIDE_hum_percent)) +
  geom_point(color = "#3a86d4", size = 2.5, alpha = 0.7) +  # Plot data points
  geom_smooth(method = "lm", color = "red", fill = "lightgray", alpha = 0.4) + # Trend line + 95% CI
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 2) + # Auto-calculated correlation label
  theme_minimal() + 
  theme(text = element_text(size = 8))


p4 <- DURH2WIDE %>% ggplot(aes(x = DURH_average_wind_speed_ms, y = WIDE_avg_wind_speed_m_s)) +
  geom_point(color = "#3a86d4", size = 2.5, alpha = 0.7) +  # Plot data points
  geom_smooth(method = "lm", color = "red", fill = "lightgray", alpha = 0.4) + # Trend line + 95% CI
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 2) + # Auto-calculated correlation label
  theme_minimal() + 
  theme(text = element_text(size = 8))

 p2 + p3 + p4 
ggsave("results/DURH2WIDE_Comparison.pdf", width = 6.5, height = 4)
