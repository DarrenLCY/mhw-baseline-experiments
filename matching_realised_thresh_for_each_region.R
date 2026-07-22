# Script to produce fig. to match the realised thresholds to the severity threshold of each region

#Load packages----
pkgs <- c("tidyverse", "lubridate", "car", "visreg", "cowplot", "survminer", "survival",
          "emmeans", "lme4")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)
library(minpack.lm)
library(grid)
library(gridExtra)

#Load packages----
pkgs <- c("rTPC", "nls.multstart", "broom", "tidyverse", "cowplot", "ggrepel")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

rm(list = ls())
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/")

#Load csvs ----
growth_base <- read.csv("growth_variables.csv")
# growth_base <- read.csv("growth_variables_corrected.csv")
# meso_survival <- read.csv("data/snail_RVs/meso_collated_survival.csv")
# meso_feed <- read.csv("data/snail_RVs/meso_percap_feeding_rate.csv")

pre_T11_logger_data <- read.csv("Loggers data/D11 2025-07-26 11_24_45 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T12.5_logger_data <- read.csv("Loggers data/D125 2025-07-26 17_15_02 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T14_logger_data <- read.csv("Loggers data/D14 2025-07-24 21_03_35 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T15.5_logger_data <- read.csv("Loggers data/D155 2025-07-24 16_46_34 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T17_logger_data <- read.csv("Loggers data/D17 2025-07-25 14_58_41 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T18_logger_data <- read.csv("Loggers data/D18 2025-07-25 18_11_00 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T19_logger_data <- read.csv("Loggers data/D19 2025-07-24 11_17_09 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T20_logger_data <- read.csv("Loggers data/D20 2025-07-26 21_01_04 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T20.5_logger_data <- read.csv("Loggers data/D205 2025-07-28 15_27_32 PDT (Data PDT).csv", fileEncoding = "latin1")
T21_final_logger_data <- read.csv("Loggers data/D21 2025-07-25 11_23_35 PDT (Data PDT).csv", fileEncoding = "latin1")
pre_T21.5_logger_data <- read.csv("Loggers data/D215 2025-07-26 15_10_23 PDT (Data PDT).csv", fileEncoding = "latin1")
T22_final_logger_data <- read.csv("Loggers data/D22 2025-07-24 15_02_17 PDT (Data PDT).csv", fileEncoding = "latin1")

# these 2 loggers data were downloaded from broken loggers (which stopped recording after some time), after which they were replaced
T21_first_logger_data <- read.csv("Loggers data/D21 2025-06-25 18_35_22 PDT (Data PDT).csv", fileEncoding = "latin1")
T22_first_logger_data <- read.csv("Loggers data/D22 2025-06-28 10_30_52 PDT (Data PDT).csv", fileEncoding = "latin1")

# merge for T21 and T22
pre_T21_logger_data = rbind(T21_first_logger_data, T21_final_logger_data)
pre_T22_logger_data = rbind(T22_first_logger_data, T22_final_logger_data)

# ------------- clean loggers data
# NOTES:
# On 6th June 2025 (Friday) at night, snails were moved from sea tables (constant temperature of ~ 12.3 degC) to 
# incubators, set at constant temperatures of 14 deg for high temperature treatments (top 5) and 13 deg otherwise

# 7th June: Temperatures were changed at a max rate of 2 deg/day for 2 days to reach the climatological temperatures starting values
# If you want to check temperatures of loggers and how constant they are (i.e., at 13 and 14 degC), 
# check from morning of 7th to around 2pm when I started changing the temperatures 

# acclimation at seasonal temperatures from 9 to 13th June 2025

# wk1 feeding rate - 20 June 2025
# safe to discard all temperatures prior to 13th June 2025, which is when the experiments started
# remove temperature spikes which occurred when displacing and rearranging tanks randomly: wk1, 3 and 5
# wk 1 displacement: all on the same day, on 20 June 2025 
# wk 3 displacement: 
# Thursday 3 Jul: T19, 22, 15.5, 14, 17,18
# Friday 4 Jul: T20,21,21.5,20.5,11,12.5
# wk 5 displacement: 
# Thursday 17 Jul: T19, 22, 15.5, 14, 17,18
# Friday 18 Jul: T20,21,21.5,20.5,11,12.5 

T11_logger_data = pre_T11_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T11_logger_data$Temperature....C[1083:1133] = NA
T11_logger_data$Temperature....C[2968:3020] = NA # remove 3rd Jul - likely displaced T11 to make space for others
T11_logger_data$Temperature....C[3109:3168] = NA # remove 4 Jul - displacement
T11_logger_data$Temperature....C[4966:5025] = NA # remove 17 Jul - likely displaced T11 to make space for others
T11_logger_data$Temperature....C[4559:4580] = NA # remove random fluctuations likely due to water changes
T11_logger_data$Temperature....C[5570:5599] = NA # remove random fluctuations likely due to water changes

T12.5_logger_data = pre_T12.5_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T12.5_logger_data$Temperature....C[1129:1165] = NA
T12.5_logger_data$Temperature....C[3139:3195] = NA # remove 4 Jul - displacement
T12.5_logger_data$Temperature....C[4995:5040] = NA # remove 17 Jul

T14_logger_data = pre_T14_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T14_logger_data$Temperature....C[1095:1144] = NA
T14_logger_data$Temperature....C[2984:3035] = NA # remove 3 Jul - displacement
T14_logger_data$Temperature....C[5132:5160] = NA # remove 18 Jul - displacement

T15.5_logger_data = pre_T15.5_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T15.5_logger_data$Temperature....C[1136:1175] = NA
T15.5_logger_data$Temperature....C[2954:2976] = NA # remove 3 Jul - displacement
T15.5_logger_data$Temperature....C[4999:5065] = NA # remove 17 Jul - displacement
T15.5_logger_data$Temperature....C[5110:5160] = NA # remove 18 Jul - displacement
T15.5_logger_data$Temperature....C[2553:2600] = NA # remove random fluctuations likely due to water changes

T17_logger_data = pre_T17_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T17_logger_data$Temperature....C[1133:1165] = NA
T17_logger_data$Temperature....C[2968:3000] = NA

T18_logger_data = pre_T18_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T18_logger_data$Temperature....C[1103:1130] = NA
T18_logger_data$Temperature....C[4993:5025] = NA # remove 17 Jul - displacement

T19_logger_data = pre_T19_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T19_logger_data$Temperature....C[1104:1160] = NA

T20_logger_data = pre_T20_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T20_logger_data$Temperature....C[1099:1130] = NA
T20_logger_data$Temperature....C[3124:3155] = NA
T20_logger_data$Temperature....C[4971:5034] = NA # remove 17 Jul - displacement

T20.5_logger_data = pre_T20.5_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T20.5_logger_data$Temperature....C[1102:1127] = NA
T20.5_logger_data$Temperature....C[3146:3180] = NA
T20.5_logger_data$Temperature....C[5018:5055] = NA # remove 17 Jul - displacement

T21_logger_data = pre_T21_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T21_logger_data$Temperature....C[214:218] = NA # correction cos of change in logger
T21_logger_data$Temperature....C[1500:1560] = NA

T21.5_logger_data = pre_T21.5_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T21.5_logger_data$Temperature....C[1085:1096] = NA

T22_logger_data = pre_T22_logger_data %>% 
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
  filter(Date.Time..PDT. > as.Date("2025-06-13"))
T22_logger_data$Temperature....C[2126:2148] = NA # correction cos of change in logger
T22_logger_data$Temperature....C[2988:3000] = NA

# this part is for plotting temps prior to expt and to remove the first 40 values for plotting temps prior to expt
pre_T11_logger_data = pre_T11_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T12.5_logger_data = pre_T12.5_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T14_logger_data = pre_T14_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T15.5_logger_data = pre_T15.5_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T17_logger_data = pre_T17_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T18_logger_data = pre_T18_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T19_logger_data = pre_T19_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T20_logger_data = pre_T20_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T20.5_logger_data = pre_T20.5_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T21_logger_data = pre_T21_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T21.5_logger_data = pre_T21.5_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
pre_T22_logger_data = pre_T22_logger_data %>%
  mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>%
  filter(Date.Time..PDT. <= as.Date("2025-06-13"))
T11_logger_data = rbind(pre_T11_logger_data, T11_logger_data)
T12.5_logger_data = rbind(pre_T12.5_logger_data, T12.5_logger_data)
T14_logger_data = rbind(pre_T14_logger_data, T14_logger_data)
T15.5_logger_data = rbind(pre_T15.5_logger_data, T15.5_logger_data)
T17_logger_data = rbind(pre_T17_logger_data, T17_logger_data)
T18_logger_data = rbind(pre_T18_logger_data, T18_logger_data)
T19_logger_data = rbind(pre_T19_logger_data, T19_logger_data)
T20_logger_data = rbind(pre_T20_logger_data, T20_logger_data)
T20.5_logger_data = rbind(pre_T20.5_logger_data, T20.5_logger_data)
T21_logger_data = rbind(pre_T21_logger_data, T21_logger_data)
T21.5_logger_data = rbind(pre_T21.5_logger_data, T21.5_logger_data)
T22_logger_data = rbind(pre_T22_logger_data, T22_logger_data)
T11_logger_data$Temperature....C[1:40] = NA
T12.5_logger_data$Temperature....C[1:40] = NA
T14_logger_data$Temperature....C[1:40] = NA
T15.5_logger_data$Temperature....C[1:40] = NA
T17_logger_data$Temperature....C[1:40] = NA
T18_logger_data$Temperature....C[1:40] = NA
T19_logger_data$Temperature....C[1:40] = NA
T20_logger_data$Temperature....C[1:40] = NA
T20.5_logger_data$Temperature....C[1:40] = NA
T21_logger_data$Temperature....C[1:40] = NA
T21.5_logger_data$Temperature....C[1:40] = NA
T22_logger_data$Temperature....C[1:40] = NA

ggplot(T11_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T12.5_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T14_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T15.5_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T17_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T18_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T19_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T20_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T20.5_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T21_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T21.5_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T22_logger_data, aes(x = Date.Time..PDT., y = Temperature....C)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 

# after cleaning, we can compute the daily average temps
T11_logger_data_daily <- T11_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "11")
T12.5_logger_data_daily <- T12.5_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "12.5")
T14_logger_data_daily <- T14_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "14")
T15.5_logger_data_daily <- T15.5_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "15.5")
T17_logger_data_daily <- T17_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "17")
T18_logger_data_daily <- T18_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "18")
T19_logger_data_daily <- T19_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "19")
T20_logger_data_daily <- T20_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "20")
T20.5_logger_data_daily <- T20.5_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "20.5")
# for T21, let's add the missing dates and best estimate of missed temperatures
# Create full sequence of dates
T21_all_dates <- data.frame(date = seq.Date(from = as.Date("2025-06-13"), # case when we add temps prior to expt for plotting... change to 6 June
                                            to = as.Date("2025-07-25"),
                                            by = "day"))
T21_set_temps = read.csv("Programmed loggers temps/T21_only.csv")
T21_logger_data_daily <- T21_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  right_join(T21_all_dates, by = "date") %>%
  arrange(date)  %>%
  mutate(Treat = "21")
# T21_set_temps$diff = T21_set_temps$sev21_avg - T21_logger_data_daily$daily_avg_temp[8:50] # case when we add temps prior to expt for plotting... COMMENT OUT normally
T21_set_temps$diff = T21_set_temps$sev21_avg - T21_logger_data_daily$daily_avg_temp
T21_set_temps$diff[3:7] = -0.5 # we have 2 different systematic errors cos of different incubators (changed on 20 June)
T21_set_temps$diff[8:12] = -0.1
T21_logger_data_daily$daily_avg_temp[3:12] = T21_set_temps$sev21_avg[3:12] - T21_set_temps$diff[3:12]
# T21_logger_data_daily$daily_avg_temp[10:19] = T21_set_temps$sev21_avg[10:19] - T21_set_temps$diff[10:19] # case when we add temps prior to expt for plotting... COMMENT OUT normally

T21.5_logger_data_daily <- T21.5_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "21.5")
T22_logger_data_daily <- T22_logger_data %>%
  mutate(date = as.Date(format(Date.Time..PDT., "%Y-%m-%d"))) %>% # Extract date
  group_by(date) %>% # Group by date
  summarise(daily_avg_temp = mean(Temperature....C, na.rm = TRUE)) %>%
  mutate(Treat = "22")

ggplot(T11_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T12.5_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T14_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T15.5_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T17_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T18_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T19_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T20_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T20.5_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T21_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T21.5_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 
ggplot(T22_logger_data_daily, aes(x = date, y = daily_avg_temp)) + geom_line() +
  labs(title = "Temperature Over Time", x = "Date", y = "Temperature (°C)") 

# combine all daily temps into one df
all_daily_temps = rbind(T11_logger_data_daily, T12.5_logger_data_daily, T14_logger_data_daily, T15.5_logger_data_daily,
                        T17_logger_data_daily, T18_logger_data_daily, T19_logger_data_daily, T20_logger_data_daily,
                        T20.5_logger_data_daily, T21_logger_data_daily, T21.5_logger_data_daily, T22_logger_data_daily)

# ---------------- match the daily temps to the Hobday thresholds to determine the 'realised' severity threshold
# we will match them for 3 different periods (weeks 0-2, 2-4 and 4-6), due to the variable realised temperature

#let's see how it looks like first... Run experimental_temps.R until you have the final plot (located in C:\Users\dlcyli\OneDrive\Development of thesis\Nucella experiments)
# then combine the final plot with the realised temps
# change the dates from 2020 to 2025...
ts_sev11 <- ts_sev11 %>% mutate(t = update(t, year = 2025))
ts_sev12.5 <- ts_sev12.5 %>% mutate(t = update(t, year = 2025))
ts_sev14 <- ts_sev14 %>% mutate(t = update(t, year = 2025))
ts_sev15.5 <- ts_sev15.5 %>% mutate(t = update(t, year = 2025))
ts_sev17 <- ts_sev17 %>% mutate(t = update(t, year = 2025))
ts_sev18 <- ts_sev18 %>% mutate(t = update(t, year = 2025))
ts_sev19 <- ts_sev19 %>% mutate(t = update(t, year = 2025))
ts_sev20 <- ts_sev20 %>% mutate(t = update(t, year = 2025))
ts_sev20.5 <- ts_sev20.5 %>% mutate(t = update(t, year = 2025))
ts_sev21 <- ts_sev21 %>% mutate(t = update(t, year = 2025))
ts_sev21.5 <- ts_sev21.5 %>% mutate(t = update(t, year = 2025))
ts_sev22 <- ts_sev22 %>% mutate(t = update(t, year = 2025))
check_point_dates = c("2025-06-09", "2025-06-13", "2025-06-26", "2025-07-10", "2025-07-26")
p = ggplot() +
  geom_line(data = ts_sev15.5, aes(x = t, y = sev15.5_avg), colour = "red", size = 1.5) +
  geom_line(data = T15.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_rect(aes(xmin = as.Date("2025-06-13"), xmax = as.Date("2025-07-26"), ymin = - Inf,
                ymax = Inf, fill = "red"), alpha = 0.2)  +
  geom_vline(xintercept = as.Date(check_point_dates)) +
  labs(y = "Temperature (°C)", x = "Date") +
  guides(fill = 'none')  +
  ggtitle("Averaged") +
  theme_classic(base_size = 25) +
  theme(plot.title = element_blank(), 
        text = element_text(size=25)) +
  coord_cartesian(ylim = c(8, 25))

sog_sev11 <- sog_sev11 %>% mutate(t = update(t, year = 2025))
sog_sev12.5 <- sog_sev12.5 %>% mutate(t = update(t, year = 2025))
sog_sev14 <- sog_sev14 %>% mutate(t = update(t, year = 2025))
sog_sev15.5 <- sog_sev15.5 %>% mutate(t = update(t, year = 2025))
sog_sev17 <- sog_sev17 %>% mutate(t = update(t, year = 2025))
sog_sev18 <- sog_sev18 %>% mutate(t = update(t, year = 2025))
sog_sev19 <- sog_sev19 %>% mutate(t = update(t, year = 2025))
sog_sev20 <- sog_sev20 %>% mutate(t = update(t, year = 2025))
sog_sev20.5 <- sog_sev20.5 %>% mutate(t = update(t, year = 2025))
sog_sev21 <- sog_sev21 %>% mutate(t = update(t, year = 2025))
sog_sev21.5 <- sog_sev21.5 %>% mutate(t = update(t, year = 2025))
sog_sev22 <- sog_sev22 %>% mutate(t = update(t, year = 2025))

cc_sev11 <- cc_sev11 %>% mutate(t = update(t, year = 2025))
cc_sev12.5 <- cc_sev12.5 %>% mutate(t = update(t, year = 2025))
cc_sev14 <- cc_sev14 %>% mutate(t = update(t, year = 2025))
cc_sev15.5 <- cc_sev15.5 %>% mutate(t = update(t, year = 2025))
cc_sev17 <- cc_sev17 %>% mutate(t = update(t, year = 2025))
cc_sev18 <- cc_sev18 %>% mutate(t = update(t, year = 2025))
cc_sev19 <- cc_sev19 %>% mutate(t = update(t, year = 2025))
cc_sev20 <- cc_sev20 %>% mutate(t = update(t, year = 2025))
cc_sev20.5 <- cc_sev20.5 %>% mutate(t = update(t, year = 2025))
cc_sev21 <- cc_sev21 %>% mutate(t = update(t, year = 2025))
cc_sev21.5 <- cc_sev21.5 %>% mutate(t = update(t, year = 2025))
cc_sev22 <- cc_sev22 %>% mutate(t = update(t, year = 2025))

#Estimate shell weight (ShW) & tissue weight (TiW) based on the following submerged regressions for each population where x is SW (submerged weight):
#Pruth	y = 1.5889x + 0.1392
#Kwakshua	y = 1.5958x + 0.0646
#Cedar	y = 1.61x + 0.0266
#Heron	y = 1.6104x - .1292
#Calculate TiW (tissue weight) based on Shell weight and Total weight
#Remove any rows with NAs in the L, T, SG, TiW & ShW columns (these died, but I kept them in my original datasheet)
#Separate the Treat into Temp & pH columns

pruth_reg <- function(x){
  ShW <- 2.1027890*x - 0.2122236
  return(ShW)
}
kwak_reg <- function(x){
  ShW <- 1.90570714*x + 0.01766196
  return(ShW)
}
cedar_reg <- function(x){
  ShW <- 1.71697847*x - 0.02090532
  return(ShW)
}
heron_reg <- function(x){
  ShW <- 1.74215419*x - 0.05069835
  return(ShW)
}

growth_clean <- growth_base %>% 
  rename(L = Length, SG = Linear_shell_growth, TW = Total_weight, SW = Submerged_weight) %>% 
  mutate(ShW = ifelse(SP == "Cedar", cedar_reg(SW), 
                      ifelse(SP == "Heron", heron_reg(SW),
                             ifelse(SP == "Kwak", kwak_reg(SW),
                                    ifelse(SP == "Pruth", pruth_reg(SW), NA)))),
         SR = as.factor(ifelse(SP == "Cedar" | SP == "Heron", "Strait of Georgia", "Central Coast")),
         SP = as.factor(SP),
         Stage = as.factor(Stage), 
         TiW = TW-ShW)

growth_rates = data.frame(ID = vector(), Treat = vector(), SP = vector(), SR = vector(), 
                          wk2_L_growth = vector(), wk4_L_growth = vector(), wk6_L_growth = vector(), 
                          wk2_SG_growth = vector(), wk4_SG_growth = vector(), wk6_SG_growth = vector(), 
                          wk2_ShW_growth = vector(), wk4_ShW_growth = vector(), wk6_ShW_growth = vector(), 
                          wk2_TiW_growth = vector(), wk4_TiW_growth = vector(), wk6_TiW_growth = vector())
for(treat in unique(growth_clean$Treat)) {
  treat_temp_df = growth_clean[which(growth_clean$Treat == treat),]
  wk0_temp_df = treat_temp_df[which(treat_temp_df$Stage == "wk0"),]
  wk2_temp_df = treat_temp_df[which(treat_temp_df$Stage == "wk2"),]
  wk4_temp_df = treat_temp_df[which(treat_temp_df$Stage == "wk4"),]
  wk6_temp_df = treat_temp_df[which(treat_temp_df$Stage == "wk6"),]
  growth_rates = rbind(growth_rates, data.frame(ID = wk0_temp_df$ID, Treat = wk0_temp_df$Treat, 
                                                SP = wk0_temp_df$SP, SR = wk0_temp_df$SR,
                                                wk2_L_growth = wk2_temp_df$L - wk0_temp_df$L, 
                                                wk4_L_growth = wk4_temp_df$L - wk2_temp_df$L, 
                                                wk6_L_growth = wk6_temp_df$L - wk4_temp_df$L, 
                                                wk2_SG_growth = wk2_temp_df$SG - wk0_temp_df$SG, 
                                                wk4_SG_growth = wk4_temp_df$SG - wk2_temp_df$SG, 
                                                wk6_SG_growth = wk6_temp_df$SG - wk4_temp_df$SG, 
                                                wk2_ShW_growth = wk2_temp_df$ShW - wk0_temp_df$ShW, 
                                                wk4_ShW_growth = wk4_temp_df$ShW - wk2_temp_df$ShW, 
                                                wk6_ShW_growth = wk6_temp_df$ShW - wk4_temp_df$ShW, 
                                                wk2_TiW_growth = wk2_temp_df$TiW - wk0_temp_df$TiW, 
                                                wk4_TiW_growth = wk4_temp_df$TiW - wk2_temp_df$TiW, 
                                                wk6_TiW_growth = wk6_temp_df$TiW - wk4_temp_df$TiW))
}

growth_rates_avg <- growth_rates %>% 
  unite(unique_ID, c(SP, Treat), sep = "_", remove = FALSE) %>% 
  group_by(SP, Treat) %>% 
  summarize(meanL_growth_wk2 = mean(wk2_L_growth, na.rm = TRUE), sdL_wk2 = sd(wk2_L_growth, na.rm = TRUE),
            meanL_growth_wk4 = mean(wk4_L_growth, na.rm = TRUE), sdL_wk4 = sd(wk4_L_growth, na.rm = TRUE), 
            meanL_growth_wk6 = mean(wk6_L_growth, na.rm = TRUE), sdL_wk6 = sd(wk6_L_growth, na.rm = TRUE),
            meanSG_growth_wk2 = mean(wk2_SG_growth, na.rm = TRUE), sdSG_wk2 = sd(wk2_SG_growth, na.rm = TRUE),
            meanSG_growth_wk4 = mean(wk4_SG_growth, na.rm = TRUE), sdSG_wk4 = sd(wk4_SG_growth, na.rm = TRUE), 
            meanSG_growth_wk6 = mean(wk6_SG_growth, na.rm = TRUE), sdSG_wk6 = sd(wk6_SG_growth, na.rm = TRUE),
            meanShW_growth_wk2 = mean(wk2_ShW_growth, na.rm = TRUE), sdShW_wk2 = sd(wk2_ShW_growth, na.rm = TRUE),
            meanShW_growth_wk4 = mean(wk4_ShW_growth, na.rm = TRUE), sdShW_wk4 = sd(wk4_ShW_growth, na.rm = TRUE), 
            meanShW_growth_wk6 = mean(wk6_ShW_growth, na.rm = TRUE), sdShW_wk6 = sd(wk6_ShW_growth, na.rm = TRUE),
            meanTiW_growth_wk2 = mean(wk2_TiW_growth, na.rm = TRUE), sdTiW_wk2 = sd(wk2_TiW_growth, na.rm = TRUE),
            meanTiW_growth_wk4 = mean(wk4_TiW_growth, na.rm = TRUE), sdTiW_wk4 = sd(wk4_TiW_growth, na.rm = TRUE), 
            meanTiW_growth_wk6 = mean(wk6_TiW_growth, na.rm = TRUE), sdTiW_wk6 = sd(wk6_TiW_growth, na.rm = TRUE),n = n()) %>% 
  ungroup() %>%
  mutate(SR = as.factor(ifelse(SP == "Cedar" | SP == "Heron", "Strait of Georgia", "Central Coast")))

all_daily_temps$Treat = as.numeric(all_daily_temps$Treat)

L_rates_avg_with_temps = data.frame(period = vector(), temp = vector(), L_growth = vector(), L_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
                                    starting_date = vector(), ending_date = vector(), difference_days = vector())
all_sevs = seq(-3,10,by=0.1)
for (i in 1:3) {
  weekly_temp_period_i = data.frame(group = vector(), start_date = vector(), end_date = vector(), Treat = vector(), 
                                    avg_temp = vector(), sog_sev_thresh = vector(), cc_sev_thresh = vector(), diff_days = vector())
  
  for (j in unique(all_daily_temps$Treat)) {
    all_daily_temps_j = all_daily_temps[which(all_daily_temps$Treat == j),]
    if (i == 1) {
      growth_clean_filtered = growth_clean[which((growth_clean$Stage == "wk0" |
                                                    growth_clean$Stage == "wk2") &
                                                   growth_clean$Treat == j),]
      
    } else if (i == 2) {
      growth_clean_filtered = growth_clean[which((growth_clean$Stage == "wk2" |
                                                    growth_clean$Stage == "wk4") &
                                                   growth_clean$Treat == j),]
    } else {
      growth_clean_filtered = growth_clean[which((growth_clean$Stage == "wk4" |
                                                    growth_clean$Stage == "wk6") &
                                                   growth_clean$Treat == j),]
    }
    
    current_diff = 100
    subset_target_temp = sog_sev11 %>% filter(t >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") & t <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
    subset_realised_temp = all_daily_temps_j %>% filter(date >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") & date <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
    for (sev in all_sevs) {
      sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
      if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
        sog_chosen_sev = sev
        current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
      }
    }
    
    current_diff = 100
    subset_target_temp = cc_sev11 %>% filter(t >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") & t <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
    for (sev in all_sevs) {
      sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
      if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
        cc_chosen_sev = sev
        current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
      }
    }
    
    d_days = as.numeric(as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y")-
                          as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y"))
    weekly_temp_period_i = rbind(weekly_temp_period_i, data.frame(group = i, start_date = growth_clean_filtered$Date[1], 
                                                                  end_date = growth_clean_filtered$Date[length(growth_clean_filtered$Date)], 
                                                                  Treat = j, 
                                                                  avg_temp = mean(all_daily_temps_j$daily_avg_temp[
                                                                    which(all_daily_temps_j$date >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") &
                                                                            all_daily_temps_j$date <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
                                                                  ]), sog_sev_thresh = sog_chosen_sev, cc_sev_thresh = cc_chosen_sev, diff_days = d_days))
    
  }
  
  all_sevs2 = c(weekly_temp_period_i$sog_sev_thresh, weekly_temp_period_i$sog_sev_thresh, weekly_temp_period_i$cc_sev_thresh, weekly_temp_period_i$cc_sev_thresh)
  L_rates_avg_with_temps = rbind(L_rates_avg_with_temps, data.frame(period = i, temp = weekly_temp_period_i$avg_temp, 
                                                                    L_growth = growth_rates_avg[[(3+(i-1)*2)]], L_growth_std_dev = growth_rates_avg[[(4+(i-1)*2)]], 
                                                                    Treat = weekly_temp_period_i$Treat,
                                                                    SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                    starting_date = weekly_temp_period_i$start_date, ending_date = weekly_temp_period_i$end_date, 
                                                                    difference_days = weekly_temp_period_i$diff_days))
}

L_rates_avg_with_temps$L_growth_standardised = L_rates_avg_with_temps$L_growth/L_rates_avg_with_temps$difference_days*14
L_rates_avg_with_temps$L_growth_sd_standardised = L_rates_avg_with_temps$L_growth_std_dev/L_rates_avg_with_temps$difference_days*14

# whole experiment --------------------- 
growth_clean <- growth_base %>% 
  rename(L = Length, SG = Linear_shell_growth, TW = Total_weight, SW = Submerged_weight) %>% 
  mutate(ShW = ifelse(SP == "Cedar", cedar_reg(SW), 
                      ifelse(SP == "Heron", heron_reg(SW),
                             ifelse(SP == "Kwak", kwak_reg(SW),
                                    ifelse(SP == "Pruth", pruth_reg(SW), NA)))),
         SR = as.factor(ifelse(SP == "Cedar" | SP == "Heron", "Strait of Georgia", "Central Coast")),
         SP = as.factor(SP),
         Stage = as.factor(Stage), 
         TiW = TW-ShW)

growth_rates = data.frame(ID = vector(), Treat = vector(), SP = vector(), SR = vector(), 
                          cum_L_growth = vector(), cum_SG_growth = vector(), 
                          cum_ShW_growth = vector(), cum_TiW_growth = vector())
for(treat in unique(growth_clean$Treat)) {
  treat_temp_df = growth_clean[which(growth_clean$Treat == treat),]
  wk0_temp_df = treat_temp_df[which(treat_temp_df$Stage == "wk0"),]
  wk6_temp_df = treat_temp_df[which(treat_temp_df$Stage == "wk6"),]
  growth_rates = rbind(growth_rates, data.frame(ID = wk0_temp_df$ID, Treat = wk0_temp_df$Treat, 
                                                SP = wk0_temp_df$SP, SR = wk0_temp_df$SR,
                                                cum_L_growth = wk6_temp_df$L - wk0_temp_df$L,
                                                cum_SG_growth = wk6_temp_df$SG - wk0_temp_df$SG, 
                                                cum_ShW_growth = wk6_temp_df$ShW - wk0_temp_df$ShW,
                                                cum_TiW_growth = wk6_temp_df$TiW - wk0_temp_df$TiW))
}

growth_rates_avg <- growth_rates %>% 
  unite(unique_ID, c(SP, Treat), sep = "_", remove = FALSE) %>% 
  group_by(SP, Treat) %>% 
  summarize(meanL_growth_cum = mean(cum_L_growth, na.rm = TRUE), sdL_cum = sd(cum_L_growth, na.rm = TRUE),
            meanSG_growth_cum = mean(cum_SG_growth, na.rm = TRUE), sdSG_cum = sd(cum_SG_growth, na.rm = TRUE),
            meanShW_growth_cum = mean(cum_ShW_growth, na.rm = TRUE), sdShW_cum = sd(cum_ShW_growth, na.rm = TRUE),
            meanTiW_growth_cum = mean(cum_TiW_growth, na.rm = TRUE), sdTiW_cum = sd(cum_TiW_growth, na.rm = TRUE),
            n = n()) %>% 
  ungroup() %>%
  mutate(SR = as.factor(ifelse(SP == "Cedar" | SP == "Heron", "Strait of Georgia", "Central Coast")))

all_daily_temps$Treat = as.numeric(all_daily_temps$Treat)

L_rates_avg_with_temps_whole_period = data.frame(temp = vector(), L_growth = vector(), L_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
                                                 starting_date = vector(), ending_date = vector(), difference_days = vector())
all_sevs = seq(-3,10,by=0.1)

sev_thresh_by_treat = data.frame(start_date = vector(), end_date = vector(), Treat = vector(), 
                                 avg_temp = vector(), sog_sev_thresh = vector(), cc_sev_thresh = vector(), diff_days = vector())

for (j in unique(all_daily_temps$Treat)) {
  all_daily_temps_j = all_daily_temps[which(all_daily_temps$Treat == j),]
  growth_clean_filtered = growth_clean[which((growth_clean$Stage == "wk0" |
                                                growth_clean$Stage == "wk6") &
                                               growth_clean$Treat == j),]
  
  current_diff = 100
  subset_target_temp = sog_sev11 %>% filter(t >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") & t <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
  subset_realised_temp = all_daily_temps_j %>% filter(date >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") & date <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      sog_chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  
  current_diff = 100
  subset_target_temp = cc_sev11 %>% filter(t >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") & t <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      cc_chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  
  d_days = as.numeric(as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y")-
                        as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y"))
  sev_thresh_by_treat = rbind(sev_thresh_by_treat, data.frame(start_date = growth_clean_filtered$Date[1], 
                                                              end_date = growth_clean_filtered$Date[length(growth_clean_filtered$Date)], 
                                                              Treat = j, 
                                                              avg_temp = mean(all_daily_temps_j$daily_avg_temp[
                                                                which(all_daily_temps_j$date >= as.Date(growth_clean_filtered$Date[1], format = "%d-%b-%y") &
                                                                        all_daily_temps_j$date <= as.Date(growth_clean_filtered$Date[length(growth_clean_filtered$Date)], format = "%d-%b-%y"))
                                                              ]), sog_sev_thresh = sog_chosen_sev, cc_sev_thresh = cc_chosen_sev, diff_days = d_days))
  
}

all_sevs2 = c(sev_thresh_by_treat$sog_sev_thresh, sev_thresh_by_treat$sog_sev_thresh, sev_thresh_by_treat$cc_sev_thresh, sev_thresh_by_treat$cc_sev_thresh)
L_rates_avg_with_temps_whole_period = rbind(L_rates_avg_with_temps_whole_period, data.frame(temp = sev_thresh_by_treat$avg_temp, 
                                                                                            L_growth = growth_rates_avg$meanL_growth_cum, L_growth_std_dev = growth_rates_avg$sdL_cum, 
                                                                                            Treat = sev_thresh_by_treat$Treat,
                                                                                            SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                                            starting_date = sev_thresh_by_treat$start_date, ending_date = sev_thresh_by_treat$end_date, 
                                                                                            difference_days = sev_thresh_by_treat$diff_days))

L_rates_avg_with_temps_whole_period$L_growth_standardised = L_rates_avg_with_temps_whole_period$L_growth/L_rates_avg_with_temps_whole_period$difference_days*42
L_rates_avg_with_temps_whole_period$L_growth_sd_standardised = L_rates_avg_with_temps_whole_period$L_growth_std_dev/L_rates_avg_with_temps_whole_period$difference_days*42


whole_period_alpha_warm = L_rates_avg_with_temps_whole_period$sev_thresh[which(L_rates_avg_with_temps_whole_period$Treat == "15.5" &
                                                                                 L_rates_avg_with_temps_whole_period$SR == "Strait of Georgia")][1]
whole_period_alpha_cold = L_rates_avg_with_temps_whole_period$sev_thresh[which(L_rates_avg_with_temps_whole_period$Treat == "15.5" &
                                                                                 L_rates_avg_with_temps_whole_period$SR == "Central Coast")][1]
whole_period_thresh_warm = sog_sev11
colnames(whole_period_thresh_warm)[4] = "sev_of_interest"
whole_period_thresh_warm$sev_of_interest = whole_period_thresh_warm$seas + whole_period_thresh_warm$thresh_less_seas*whole_period_alpha_warm
whole_period_thresh_warm <- whole_period_thresh_warm %>%
  filter(t >= as.Date("2025-06-13"),
         t <= as.Date("2025-07-26"))

whole_period_thresh_cold = cc_sev11
colnames(whole_period_thresh_cold)[4] = "sev_of_interest"
whole_period_thresh_cold$sev_of_interest = whole_period_thresh_cold$seas + whole_period_thresh_cold$thresh_less_seas*whole_period_alpha_cold
whole_period_thresh_cold <- whole_period_thresh_cold %>%
  filter(t >= as.Date("2025-06-13"),
         t <= as.Date("2025-07-26"))

# period 1
period1_alpha_warm = L_rates_avg_with_temps$sev_thresh[which(L_rates_avg_with_temps$period == 1 &
                                                               L_rates_avg_with_temps$Treat == "15.5" &
                                                               L_rates_avg_with_temps$SR == "Strait of Georgia")][1]
period1_alpha_cold = L_rates_avg_with_temps$sev_thresh[which(L_rates_avg_with_temps$period == 1 &
                                                               L_rates_avg_with_temps$Treat == "15.5" &
                                                               L_rates_avg_with_temps$SR == "Central Coast")][1]
period1_thresh_warm = sog_sev11
colnames(period1_thresh_warm)[4] = "sev_of_interest"
period1_thresh_warm$sev_of_interest = period1_thresh_warm$seas + period1_thresh_warm$thresh_less_seas*period1_alpha_warm
period1_thresh_warm <- period1_thresh_warm %>%
  filter(t >= as.Date("2025-06-13"),
         t <= as.Date("2025-06-26"))

period1_thresh_cold = cc_sev11
colnames(period1_thresh_cold)[4] = "sev_of_interest"
period1_thresh_cold$sev_of_interest = period1_thresh_cold$seas + period1_thresh_cold$thresh_less_seas*period1_alpha_cold
period1_thresh_cold <- period1_thresh_cold %>%
  filter(t >= as.Date("2025-06-13"),
         t <= as.Date("2025-06-26"))

# period 2
period2_alpha_warm = L_rates_avg_with_temps$sev_thresh[which(L_rates_avg_with_temps$period == 2 &
                                                               L_rates_avg_with_temps$Treat == "15.5" &
                                                               L_rates_avg_with_temps$SR == "Strait of Georgia")][1]
period2_alpha_cold = L_rates_avg_with_temps$sev_thresh[which(L_rates_avg_with_temps$period == 2 &
                                                               L_rates_avg_with_temps$Treat == "15.5" &
                                                               L_rates_avg_with_temps$SR == "Central Coast")][1]
period2_thresh_warm = sog_sev11
colnames(period2_thresh_warm)[4] = "sev_of_interest"
period2_thresh_warm$sev_of_interest = period2_thresh_warm$seas + period2_thresh_warm$thresh_less_seas*period2_alpha_warm
period2_thresh_warm <- period2_thresh_warm %>%
  filter(t >= as.Date("2025-06-26"),
         t <= as.Date("2025-07-10"))

period2_thresh_cold = cc_sev11
colnames(period2_thresh_cold)[4] = "sev_of_interest"
period2_thresh_cold$sev_of_interest = period2_thresh_cold$seas + period2_thresh_cold$thresh_less_seas*period2_alpha_cold
period2_thresh_cold <- period2_thresh_cold %>%
  filter(t >= as.Date("2025-06-26"),
         t <= as.Date("2025-07-10"))

# period 3
period3_alpha_warm = L_rates_avg_with_temps$sev_thresh[which(L_rates_avg_with_temps$period == 3 &
                                                               L_rates_avg_with_temps$Treat == "15.5" &
                                                               L_rates_avg_with_temps$SR == "Strait of Georgia")][1]
period3_alpha_cold = L_rates_avg_with_temps$sev_thresh[which(L_rates_avg_with_temps$period == 3 &
                                                               L_rates_avg_with_temps$Treat == "15.5" &
                                                               L_rates_avg_with_temps$SR == "Central Coast")][1]
period3_thresh_warm = sog_sev11
colnames(period3_thresh_warm)[4] = "sev_of_interest"
period3_thresh_warm$sev_of_interest = period3_thresh_warm$seas + period3_thresh_warm$thresh_less_seas*period3_alpha_warm
period3_thresh_warm <- period3_thresh_warm %>%
  filter(t >= as.Date("2025-07-10"),
         t <= as.Date("2025-07-26"))

period3_thresh_cold = cc_sev11
colnames(period3_thresh_cold)[4] = "sev_of_interest"
period3_thresh_cold$sev_of_interest = period3_thresh_cold$seas + period3_thresh_cold$thresh_less_seas*period3_alpha_cold
period3_thresh_cold <- period3_thresh_cold %>%
  filter(t >= as.Date("2025-07-10"),
         t <= as.Date("2025-07-26"))

# plot for matching severity threshs throughout the whole experiment 
p_whole = ggplot() +
  geom_line(data = whole_period_thresh_warm, aes(x = t, y = sev_of_interest), colour = "red", size = 3) +
  geom_line(data = whole_period_thresh_cold, aes(x = t, y = sev_of_interest), colour = "blue", size = 3) +
  geom_line(data = T15.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "black", size = 3) +
  
  geom_rect(aes(xmin = as.Date("2025-06-13"), xmax = as.Date("2025-07-26"),
                ymin = -Inf, ymax = Inf, fill = "red"),
            alpha = 0.2)  +
  
  geom_vline(xintercept = as.Date(check_point_dates), 
             linetype = "dashed", size = 3) +
  
  # ✅ custom x-axis ticks
  scale_x_date(breaks = as.Date(check_point_dates),
               labels = format(as.Date(check_point_dates), "%b %d")) +
  
  labs(y = "Temperature (°C)", x = "Date") +
  guides(fill = 'none')  +
  ggtitle("Matching of severity levels (whole experiment)") +
  theme_classic(base_size = 35) +
  theme(axis.title.x = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 40)) +
  coord_cartesian(ylim = c(13.5, 17))


# plot for matching severity threshs during subsets of the experiment 
p = ggplot() +
  geom_line(data = period1_thresh_warm, aes(x = t, y = sev_of_interest), colour = "red", size = 3) +
  geom_line(data = period1_thresh_cold, aes(x = t, y = sev_of_interest), colour = "blue", size = 3) +
  geom_line(data = period2_thresh_warm, aes(x = t, y = sev_of_interest), colour = "red", size = 3) +
  geom_line(data = period2_thresh_cold, aes(x = t, y = sev_of_interest), colour = "blue", size = 3) +
  geom_line(data = period3_thresh_warm, aes(x = t, y = sev_of_interest), colour = "red", size = 3) +
  geom_line(data = period3_thresh_cold, aes(x = t, y = sev_of_interest), colour = "blue", size = 3) +
  geom_line(data = T15.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "black", size = 3) +
  geom_rect(aes(xmin = as.Date("2025-06-13"), xmax = as.Date("2025-07-26"), ymin = - Inf,
                ymax = Inf, fill = "red"), alpha = 0.2)  +
  # ✅ custom x-axis ticks
  scale_x_date(breaks = as.Date(check_point_dates),
               labels = format(as.Date(check_point_dates), "%b %d"),) +
  geom_vline(xintercept = as.Date(check_point_dates),
             linetype = "dashed", size = 3) +
  labs(y = "Temperature (°C)", x = "Date") +
  guides(fill = 'none')  +
  ggtitle("Matching of severity levels (subsets of experiment)") +
  theme_classic(base_size = 35) +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size=40)) +
  coord_cartesian(ylim = c(13.5, 17))

whole_period_alpha_warm
whole_period_alpha_cold
period1_alpha_warm
period1_alpha_cold
period2_alpha_warm
period2_alpha_cold
period3_alpha_warm
period3_alpha_cold

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub")
ggsave(filename="match_realised_thresh_fig.png", height=8, width=15, plot=p)
ggsave(filename="match_realised_thresh_fig_v2.png", height=15, width=15, 
       plot=grid.arrange(p_whole, p,
                         ncol = 1))
ggsave(filename="match_realised_thresh_fig_v3.png", height=15, width=15, 
       plot=grid.arrange(p_whole, p,
                         ncol = 1))
ggsave(filename="match_realised_thresh_fig_v4.png", height=15, width=15, 
       plot=grid.arrange(p_whole, p,
                         ncol = 1))
ggsave(filename="match_realised_thresh_fig_v5.png", height=15, width=21, 
       plot=grid.arrange(p_whole, p,
                         ncol = 1))
