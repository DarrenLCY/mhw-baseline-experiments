# Script to analyze CTmax and TSM changes for snail growth response variables
# first part of the code is identical to growth_variables_v5.R script
# Actual analyses start on: # Central Coast CTmax changes ----------

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

# COMMENT OUT... this is just for plotting temps prior to expt and to remove the first 40 values for plotting temps prior to expt
# pre_T11_logger_data = pre_T11_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T12.5_logger_data = pre_T12.5_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T14_logger_data = pre_T14_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T15.5_logger_data = pre_T15.5_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T17_logger_data = pre_T17_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T18_logger_data = pre_T18_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T19_logger_data = pre_T19_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T20_logger_data = pre_T20_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T20.5_logger_data = pre_T20.5_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T21_logger_data = pre_T21_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T21.5_logger_data = pre_T21.5_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# pre_T22_logger_data = pre_T22_logger_data %>% 
#   mutate(Date.Time..PDT. = as.POSIXct(Date.Time..PDT., format = "%m/%d/%Y %H:%M:%S")) %>% 
#   filter(Date.Time..PDT. <= as.Date("2025-06-13"))
# T11_logger_data = rbind(pre_T11_logger_data, T11_logger_data)
# T12.5_logger_data = rbind(pre_T12.5_logger_data, T12.5_logger_data)
# T14_logger_data = rbind(pre_T14_logger_data, T14_logger_data)
# T15.5_logger_data = rbind(pre_T15.5_logger_data, T15.5_logger_data)
# T17_logger_data = rbind(pre_T17_logger_data, T17_logger_data)
# T18_logger_data = rbind(pre_T18_logger_data, T18_logger_data)
# T19_logger_data = rbind(pre_T19_logger_data, T19_logger_data)
# T20_logger_data = rbind(pre_T20_logger_data, T20_logger_data)
# T20.5_logger_data = rbind(pre_T20.5_logger_data, T20.5_logger_data)
# T21_logger_data = rbind(pre_T21_logger_data, T21_logger_data)
# T21.5_logger_data = rbind(pre_T21.5_logger_data, T21.5_logger_data)
# T22_logger_data = rbind(pre_T22_logger_data, T22_logger_data)
# T11_logger_data$Temperature....C[1:40] = NA
# T12.5_logger_data$Temperature....C[1:40] = NA
# T14_logger_data$Temperature....C[1:40] = NA
# T15.5_logger_data$Temperature....C[1:40] = NA
# T17_logger_data$Temperature....C[1:40] = NA
# T18_logger_data$Temperature....C[1:40] = NA
# T19_logger_data$Temperature....C[1:40] = NA
# T20_logger_data$Temperature....C[1:40] = NA
# T20.5_logger_data$Temperature....C[1:40] = NA
# T21_logger_data$Temperature....C[1:40] = NA
# T21.5_logger_data$Temperature....C[1:40] = NA
# T22_logger_data$Temperature....C[1:40] = NA

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
ggplot() +
  geom_line(data = ts_sev11, aes(x = t, y = sev11_avg), colour = "red", size = 1.5) +
  geom_line(data = T11_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev12.5, aes(x = t, y = sev12.5_avg), colour = "red", size = 1.5) +
  geom_line(data = T12.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev14, aes(x = t, y = sev14_avg), colour = "red", size = 1.5) +
  geom_line(data = T14_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev15.5, aes(x = t, y = sev15.5_avg), colour = "red", size = 1.5) +
  geom_line(data = T15.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev17, aes(x = t, y = sev17_avg), colour = "red", size = 1.5) +
  geom_line(data = T17_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev18, aes(x = t, y = sev18_avg), colour = "red", size = 1.5) +
  geom_line(data = T18_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev19, aes(x = t, y = sev19_avg), colour = "red", size = 1.5) +
  geom_line(data = T19_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev20, aes(x = t, y = sev20_avg), colour = "red", size = 1.5) +
  geom_line(data = T20_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev20.5, aes(x = t, y = sev20.5_avg), colour = "red", size = 1.5) +
  geom_line(data = T20.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev21, aes(x = t, y = sev21_avg), colour = "red", size = 1.5) +
  geom_line(data = T21_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev21.5, aes(x = t, y = sev21.5_avg), colour = "red", size = 1.5) +
  geom_line(data = T21.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_line(data = ts_sev22, aes(x = t, y = sev22_avg), colour = "red", size = 1.5) +
  geom_line(data = T22_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
  geom_rect(aes(xmin = as.Date("2025-06-13"), xmax = as.Date("2025-07-26"), ymin = - Inf,
                ymax = Inf, fill = "red"), alpha = 0.2)  +
  geom_vline(xintercept = as.Date(check_point_dates)) +
  labs(y = "Temperatures", x = "Time") +
  guides(fill = 'none')  +
  ggtitle("Averaged") +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=15)) +
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

SG_rates_avg_with_temps = data.frame(period = vector(), temp = vector(), SG_growth = vector(), SG_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
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
  SG_rates_avg_with_temps = rbind(SG_rates_avg_with_temps, data.frame(period = i, temp = weekly_temp_period_i$avg_temp, 
                                                                      SG_growth = growth_rates_avg[[(9+(i-1)*2)]], SG_growth_std_dev = growth_rates_avg[[(10+(i-1)*2)]], 
                                                                      Treat = weekly_temp_period_i$Treat,
                                                                      SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                      starting_date = weekly_temp_period_i$start_date, ending_date = weekly_temp_period_i$end_date, 
                                                                      difference_days = weekly_temp_period_i$diff_days))
}

SG_rates_avg_with_temps$SG_growth_standardised = SG_rates_avg_with_temps$SG_growth/SG_rates_avg_with_temps$difference_days*14
SG_rates_avg_with_temps$SG_growth_sd_standardised = SG_rates_avg_with_temps$SG_growth_std_dev/SG_rates_avg_with_temps$difference_days*14

ShW_rates_avg_with_temps = data.frame(period = vector(), temp = vector(), ShW_growth = vector(), ShW_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
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
  ShW_rates_avg_with_temps = rbind(ShW_rates_avg_with_temps, data.frame(period = i, temp = weekly_temp_period_i$avg_temp, 
                                                                        ShW_growth = growth_rates_avg[[(15+(i-1)*2)]], ShW_growth_std_dev = growth_rates_avg[[(16+(i-1)*2)]], 
                                                                        Treat = weekly_temp_period_i$Treat,
                                                                        SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                        starting_date = weekly_temp_period_i$start_date, ending_date = weekly_temp_period_i$end_date, 
                                                                        difference_days = weekly_temp_period_i$diff_days))
}

ShW_rates_avg_with_temps$ShW_growth_standardised = ShW_rates_avg_with_temps$ShW_growth/ShW_rates_avg_with_temps$difference_days*14
ShW_rates_avg_with_temps$ShW_growth_sd_standardised = ShW_rates_avg_with_temps$ShW_growth_std_dev/ShW_rates_avg_with_temps$difference_days*14

TiW_rates_avg_with_temps = data.frame(period = vector(), temp = vector(), TiW_growth = vector(), TiW_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
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
  TiW_rates_avg_with_temps = rbind(TiW_rates_avg_with_temps, data.frame(period = i, temp = weekly_temp_period_i$avg_temp, 
                                                                        TiW_growth = growth_rates_avg[[(21+(i-1)*2)]], TiW_growth_std_dev = growth_rates_avg[[(22+(i-1)*2)]], 
                                                                        Treat = weekly_temp_period_i$Treat,
                                                                        SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                        starting_date = weekly_temp_period_i$start_date, ending_date = weekly_temp_period_i$end_date, 
                                                                        difference_days = weekly_temp_period_i$diff_days))
}

TiW_rates_avg_with_temps$TiW_growth_standardised = TiW_rates_avg_with_temps$TiW_growth/TiW_rates_avg_with_temps$difference_days*14
TiW_rates_avg_with_temps$TiW_growth_sd_standardised = TiW_rates_avg_with_temps$TiW_growth_std_dev/TiW_rates_avg_with_temps$difference_days*14

# RV vs absolute temps --------------- 

# Shell length ----------------- 

# Create an empty list to store plots
plots_SL <- list()

topt_sev_results = data.frame(param = vector(), op_sev = vector(), conf_lower = vector(), 
                              conf_upper = vector(), SR = vector(), model = vector(), period = vector(),
                              RV = vector())

for (i in 1:length(unique(L_rates_avg_with_temps$period))) {
  
  # CC Shell length 
  
  CC_rv <- L_rates_avg_with_temps %>% 
    filter(SR == "Central Coast" & period == unique(L_rates_avg_with_temps$period)[i]) %>% 
    select(temp, L_growth_standardised, SR) %>% 
    rename("rate" = "L_growth_standardised")
  CC_rv$temp = as.numeric(as.character(CC_rv$temp))
  
  d_fits <- CC_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), CC_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Central Coast') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #CC_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  cc_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), CC_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Shell length growth (mm)',
         title = 'Central Coast') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(CC_rv, aes(Treat, L_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #CC_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_CC <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Central Coast")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
  #   geom_point(aes(temp, rate), CC_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Central Coast')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = CC_rv,
                       start = beta_start,
                       lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #CC_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_CC <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Central Coast")
  
  #CC_rv: quadratic: Estimate parameters & CI intervals
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_CC_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Central Coast",
           RV = "SL",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[1],
                                                        op_sev = ci_params_select_CC_fr$estimate[1], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[1],
                                                        SR = ci_params_select_CC_fr$SR[1],
                                                        model = ci_params_select_CC_fr$model[1],
                                                        period = i, RV = "SL"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                        op_sev = ci_params_select_CC_fr$estimate[2], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                        SR = ci_params_select_CC_fr$SR[2],
                                                        model = ci_params_select_CC_fr$model[2],
                                                        period = i, RV = "SL"))
  
  
  # SoG Shell length
  
  SoG_rv <- L_rates_avg_with_temps %>% 
    filter(SR == "Strait of Georgia" & period == unique(L_rates_avg_with_temps$period)[i]) %>% 
    select(temp, L_growth_standardised, SR) %>% 
    rename("rate" = "L_growth_standardised")
  SoG_rv$temp = as.numeric(as.character(SoG_rv$temp))
  
  d_fits <- SoG_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), SoG_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Strait of Georgia') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #SoG_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  SoG_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), SoG_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Shell length growth (mm)',
         title = 'Strait of Georgia') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(SoG_rv, aes(Treat, L_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #SoG_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_SoG <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
  #   geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Strait of Georgia')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = SoG_rv,
                       start = beta_start,
                       lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #SoG_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_SoG <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot bootstrapped CIs
  p = ggplot() +
    geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
    theme_bw(base_size = 17) +
    labs(x = 'Temperature (ºC)',
         y = 'Shell length growth (mm)',
         title = paste0('TPC during period ', unique(L_rates_avg_with_temps$period)[i]))
  
  # plot --------
  p = ggplot() +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_CC, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_SoG, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
    scale_colour_manual(values = c("skyblue", "coral")) +
    scale_fill_manual(values = c("skyblue3", "coral3")) +
    labs(x = 'Temperature (ºC)',
         y = if (i == 1) "Shell length growth (mm)" else NULL,
         col = "Source Region",
         fill = "Source Region",
         title = paste0('TPC during period ', i)) + 
    theme_cowplot(16) + 
    scale_x_continuous(breaks = c(10,12,14,16,18,20,22)) +
    expand_limits(x = c(10,22.5)) +
    scale_y_continuous(breaks = c(0,1,2,3)) +
    expand_limits(y = c(-1,3.5)) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.title.x = element_blank(),
          axis.text  = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  # Store the plot in the list
  plots_SL[[i]] <- p
  
  #SoG_rv: quadratic: Estimate parameters & CI intervals 
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_SoG_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Strait of Georgia",
           RV = "SL",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                        op_sev = ci_params_select_SoG_fr$estimate[1], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                        SR = ci_params_select_SoG_fr$SR[1],
                                                        model = ci_params_select_SoG_fr$model[1],
                                                        period = i, RV = "SL"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                        op_sev = ci_params_select_SoG_fr$estimate[2], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                        SR = ci_params_select_SoG_fr$SR[2],
                                                        model = ci_params_select_SoG_fr$model[2],
                                                        period = i, RV = "SL"))
}

# Shell lip ----------
plots_SG <- list()
for (i in 1:length(unique(SG_rates_avg_with_temps$period))) {
  
  # CC Shell lip 
  
  CC_rv <- SG_rates_avg_with_temps %>% 
    filter(SR == "Central Coast" & period == unique(SG_rates_avg_with_temps$period)[i]) %>% 
    select(temp, SG_growth_standardised, SR) %>% 
    rename("rate" = "SG_growth_standardised")
  CC_rv$temp = as.numeric(as.character(CC_rv$temp))
  
  d_fits <- CC_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), CC_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Central Coast') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #CC_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  cc_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), CC_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Shell lip growth (mm)',
         title = 'Central Coast') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(CC_rv, aes(Treat, SG_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #CC_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_CC <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Central Coast")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
  #   geom_point(aes(temp, rate), CC_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Central Coast')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = CC_rv,
                       start = beta_start,
                       lower = c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       upper = c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #CC_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_CC <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Central Coast")
  
  #CC_rv: quadratic: Estimate parameters & CI intervals
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_CC_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Central Coast",
           RV = "SG",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[1],
                                                        op_sev = ci_params_select_CC_fr$estimate[1], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[1],
                                                        SR = ci_params_select_CC_fr$SR[1],
                                                        model = ci_params_select_CC_fr$model[1],
                                                        period = i, RV = "SG"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                        op_sev = ci_params_select_CC_fr$estimate[2], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                        SR = ci_params_select_CC_fr$SR[2],
                                                        model = ci_params_select_CC_fr$model[2],
                                                        period = i, RV = "SG"))
  
  
  # SoG Shell lip
  
  SoG_rv <- SG_rates_avg_with_temps %>% 
    filter(SR == "Strait of Georgia" & period == unique(SG_rates_avg_with_temps$period)[i]) %>% 
    select(temp, SG_growth_standardised, SR) %>% 
    rename("rate" = "SG_growth_standardised")
  SoG_rv$temp = as.numeric(as.character(SoG_rv$temp))
  
  d_fits <- SoG_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), SoG_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Strait of Georgia') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #SoG_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  SoG_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), SoG_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Shell lip growth (mm)',
         title = 'Strait of Georgia') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(SoG_rv, aes(Treat, SG_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #SoG_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_SoG <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
  #   geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Strait of Georgia')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = SoG_rv,
                       start = beta_start,
                       lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #SoG_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_SoG <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot bootstrapped CIs
  p = ggplot() +
    geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
    theme_bw(base_size = 17) +
    labs(x = 'Temperature (ºC)',
         y = 'Shell lip growth (mm)',
         title = paste0('TPC during period ', unique(SG_rates_avg_with_temps$period)[i]))
  
  # plot --------
  p = ggplot() +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_CC, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_SoG, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
    scale_colour_manual(values = c("skyblue", "coral")) +
    scale_fill_manual(values = c("skyblue3", "coral3")) +
    labs(x = 'Temperature (ºC)',
         y = if (i == 1) "Shell lip growth (mm)" else NULL,
         col = "Source Region",
         fill = "Source Region",
         title = paste0('TPC during period ', i)) + 
    theme_cowplot(16) + 
    scale_x_continuous(breaks = c(10,12,14,16,18,20,22)) +
    expand_limits(x = c(10,22.5)) +
    scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14)) +
    expand_limits(y = c(-4,14)) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text  = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  # Store the plot in the list
  plots_SG[[i]] <- p
  
  #SoG_rv: quadratic: Estimate parameters & CI intervals 
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_SoG_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Strait of Georgia",
           RV = "SG",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                        op_sev = ci_params_select_SoG_fr$estimate[1], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                        SR = ci_params_select_SoG_fr$SR[1],
                                                        model = ci_params_select_SoG_fr$model[1],
                                                        period = i, RV = "SG"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                        op_sev = ci_params_select_SoG_fr$estimate[2], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                        SR = ci_params_select_SoG_fr$SR[2],
                                                        model = ci_params_select_SoG_fr$model[2],
                                                        period = i, RV = "SG"))
}

# Shell weight ----------------- 
plots_ShW <- list()
for (i in 1:length(unique(ShW_rates_avg_with_temps$period))) {
  
  # CC Shell weight 
  
  CC_rv <- ShW_rates_avg_with_temps %>% 
    filter(SR == "Central Coast" & period == unique(ShW_rates_avg_with_temps$period)[i]) %>% 
    select(temp, ShW_growth_standardised, SR) %>% 
    rename("rate" = "ShW_growth_standardised")
  CC_rv$temp = as.numeric(as.character(CC_rv$temp))
  
  d_fits <- CC_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), CC_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Central Coast') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #CC_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  cc_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), CC_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Shell weight growth (g)',
         title = 'Central Coast') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(CC_rv, aes(Treat, ShW_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #CC_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_CC <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Central Coast")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
  #   geom_point(aes(temp, rate), CC_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Central Coast')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = CC_rv,
                       start = beta_start,
                       lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #CC_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_CC <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Central Coast")
  
  #CC_rv: quadratic: Estimate parameters & CI intervals
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_CC_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Central Coast",
           RV = "ShW",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[1],
                                                        op_sev = ci_params_select_CC_fr$estimate[1], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[1],
                                                        SR = ci_params_select_CC_fr$SR[1],
                                                        model = ci_params_select_CC_fr$model[1],
                                                        period = i, RV = "ShW"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                        op_sev = ci_params_select_CC_fr$estimate[2], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                        SR = ci_params_select_CC_fr$SR[2],
                                                        model = ci_params_select_CC_fr$model[2],
                                                        period = i, RV = "ShW"))
  
  
  # SoG Shell weight
  
  SoG_rv <- ShW_rates_avg_with_temps %>% 
    filter(SR == "Strait of Georgia" & period == unique(ShW_rates_avg_with_temps$period)[i]) %>% 
    select(temp, ShW_growth_standardised, SR) %>% 
    rename("rate" = "ShW_growth_standardised")
  SoG_rv$temp = as.numeric(as.character(SoG_rv$temp))
  
  d_fits <- SoG_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), SoG_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Strait of Georgia') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #SoG_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  SoG_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), SoG_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Shell weight growth (g)',
         title = 'Strait of Georgia') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(SoG_rv, aes(Treat, ShW_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #SoG_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_SoG <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
  #   geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Strait of Georgia')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = SoG_rv,
                       start = beta_start,
                       lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #SoG_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_SoG <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot bootstrapped CIs
  p = ggplot() +
    geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
    theme_bw(base_size = 17) +
    labs(x = 'Temperature (ºC)',
         y = 'Shell weight growth (g)',
         title = paste0('TPC during period ', unique(ShW_rates_avg_with_temps$period)[i]))
  
  # plot --------
  p = ggplot() +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_CC, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_SoG, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
    scale_colour_manual(values = c("skyblue", "coral")) +
    scale_fill_manual(values = c("skyblue3", "coral3")) +
    labs(x = 'Temperature (ºC)',
         y = if (i == 1) "Shell weight growth (g)" else NULL,
         col = "Source Region",
         fill = "Source Region",
         title = paste0('TPC during period ', i)) + 
    theme_cowplot(16) + 
    scale_x_continuous(breaks = c(10,12,14,16,18,20,22)) +
    expand_limits(x = c(10,22.5)) +
    scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5)) +
    expand_limits(y = c(-0.15,0.5)) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.title.x = element_blank(),
          axis.text  = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  # Store the plot in the list
  plots_ShW[[i]] <- p
  
  #SoG_rv: quadratic: Estimate parameters & CI intervals 
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_SoG_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Strait of Georgia",
           RV = "ShW",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                        op_sev = ci_params_select_SoG_fr$estimate[1], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                        SR = ci_params_select_SoG_fr$SR[1],
                                                        model = ci_params_select_SoG_fr$model[1],
                                                        period = i, RV = "ShW"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                        op_sev = ci_params_select_SoG_fr$estimate[2], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                        SR = ci_params_select_SoG_fr$SR[2],
                                                        model = ci_params_select_SoG_fr$model[2],
                                                        period = i, RV = "ShW"))
}

# Tissue weight only ----------------- 
plots_TiW <- list()
for (i in 1:length(unique(TiW_rates_avg_with_temps$period))) {
  
  # CC Tissue weight 
  
  CC_rv <- TiW_rates_avg_with_temps %>% 
    filter(SR == "Central Coast" & period == unique(TiW_rates_avg_with_temps$period)[i]) %>% 
    select(temp, TiW_growth_standardised, SR) %>% 
    rename("rate" = "TiW_growth_standardised")
  CC_rv$temp = as.numeric(as.character(CC_rv$temp))
  
  d_fits <- CC_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), CC_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Central Coast') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #CC_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  cc_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), CC_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Tissue weight growth (g)',
         title = 'Central Coast') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(CC_rv, aes(Treat, TiW_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #CC_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(CC_rv$temp, CC_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(CC_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_CC <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Central Coast")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
  #   geom_point(aes(temp, rate), CC_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Central Coast')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = CC_rv,
                       start = beta_start,
                       lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(CC_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = CC_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(CC_rv$temp, CC_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(CC_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #CC_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(CC_rv$temp), max(CC_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_CC <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Central Coast")
  
  #CC_rv: quadratic: Estimate parameters & CI intervals
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_CC_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Central Coast",
           RV = "TiW",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[1],
                                                        op_sev = ci_params_select_CC_fr$estimate[1], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[1],
                                                        SR = ci_params_select_CC_fr$SR[1],
                                                        model = ci_params_select_CC_fr$model[1],
                                                        period = i, RV = "TiW"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                        op_sev = ci_params_select_CC_fr$estimate[2], 
                                                        conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                        SR = ci_params_select_CC_fr$SR[2],
                                                        model = ci_params_select_CC_fr$model[2],
                                                        period = i, RV = "TiW"))
  
  
  # SoG Tissue weight
  
  SoG_rv <- TiW_rates_avg_with_temps %>% 
    filter(SR == "Strait of Georgia" & period == unique(TiW_rates_avg_with_temps$period)[i]) %>% 
    select(temp, TiW_growth_standardised, SR) %>% 
    rename("rate" = "TiW_growth_standardised")
  SoG_rv$temp = as.numeric(as.character(SoG_rv$temp))
  
  d_fits <- SoG_rv %>%
    nest(data = c(temp, rate)) %>%
    mutate(
      beta = map(data, ~{
        params <- c("a","b","c","d", "e")  # parameters in your formula
        starts <- get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
        # Force valid bounds
        lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
        uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
        
        # Subset to exactly the parameters you estimate
        starts <- starts[params]
        lowers <- lowers[params]
        uppers <- uppers[params]
        
        # Defensive checks
        stopifnot(length(starts) == length(params),
                  length(lowers) == length(params),
                  length(uppers) == length(params))
        
        nls_multstart(
          rate ~ beta_2012(temp, a, b, c, d, e),
          data = .x,
          iter = rep(4, length(params)),
          start_lower = starts - 10,
          start_upper = starts + 10,
          lower = lowers,
          upper = uppers,
          supp_errors = "Y",
          convergence_count = FALSE
        )
      }),
      
      gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                          data = .x,
                                          iter = c(4,4,4),
                                          start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                          start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                          lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                          supp_errors = 'Y',
                                          convergence_count = FALSE)),
      
      quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                           data = .x,
                                           iter = c(4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, temp < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(temp, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(temp, rate), SoG_rv) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Strait of Georgia') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #SoG_rv: Now start the AIC process
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  # best model is set to quadratic everywhere
  best_model = "quadratic"
  
  # get colour code
  col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]
  
  # plot
  SoG_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), SoG_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Tissue weight growth (g)',
         title = 'Strait of Georgia') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(SoG_rv, aes(Treat, TiW_growth_standardised)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Rate')
  
  models_df <- data.frame( model = c("quadratic", "gaussian", "beta"), 
                           model_name = c("quadratic_2008", "gaussian_1987", "beta_2012") )
  
  chosen_model_name <- models_df$model_name[models_df$model == best_model]
  
  # Define parameter sets and iter lengths for each model
  model_specs <- list( quadratic = list(params = c("a", "b", "c"), 
                                        iter = c(4, 4, 4)), 
                       beta = list(params = c("a", "b", "c", "d", "e"), 
                                   iter = c(4, 4, 4, 4, 4)), 
                       gaussian = list(params = c("rmax", "topt", "a"), 
                                       iter = c(4, 4, 4)) )
  spec <- model_specs[[best_model]]
  
  # Build the formula dynamically 
  fit_formula <- as.formula( paste0("rate ~ ", chosen_model_name, 
                                    "(temp = temp, ", 
                                    paste(spec$params, collapse = ", "), ")") )
  
  #SoG_rv: Fit data
  if (best_model == "quadratic" | best_model == "gaussian") {
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( fit = map(data, ~nls_multstart( 
        formula = fit_formula, data = .x, iter = spec$iter, 
        start_lower = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) - 10, 
        start_upper = get_start_vals(.x$temp, .x$rate, model_name = chosen_model_name) + 10, 
        lower = get_lower_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        upper = get_upper_lims(.x$temp, .x$rate, model_name = chosen_model_name), 
        supp_errors = "Y", convergence_count = FALSE )), 
        # create new temperature data 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        # predict over that data 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  } else {
    beta_params <- c("a", "b", "c", "d", "e") 
    beta_start = get_start_vals(SoG_rv$temp, SoG_rv$rate, model_name = "beta_2012")
    # Force valid bounds
    lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
    uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
    starts <- beta_start[beta_params]
    lowers <- lowers[beta_params]
    uppers <- uppers[beta_params]
    
    d_fit <- nest(SoG_rv, data = c(temp, rate)) %>% 
      mutate( 
        fit = map(data, ~nls_multstart(
          formula = fit_formula, 
          data = .x, 
          iter = spec$iter, 
          start_lower = beta_start - 10, 
          start_upper = beta_start + 10, 
          lower = lowers, # or your own limits 
          upper = uppers, 
          supp_errors = "Y", convergence_count = FALSE )), 
        new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))), 
        preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
  }
  
  # unnest predictions
  d_preds_SoG <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
  #   geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC) (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Strait of Georgia')
  
  if (best_model == "quadratic") {
    fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = temp, a, b, c), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "beta") {
    beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
    fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = temp, a, b, c, d, e), 
                       data = SoG_rv,
                       start = beta_start,
                       lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                       uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                       weights = rep(1, times = nrow(SoG_rv)))
  } else if (best_model == "gaussian") {
    fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = temp, rmax, topt, a), 
                       data = SoG_rv,
                       start = coef(d_fit[["fit"]][[1]]),
                       lower = get_lower_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       upper = get_upper_lims(SoG_rv$temp, SoG_rv$rate, model_name = chosen_model_name),
                       weights = rep(1, times = nrow(SoG_rv)))
  }
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  # Get the function object 
  chosen_fun <- get(chosen_model_name)
  
  # Extract parameter names from bootstrap results 
  param_names <- colnames(boot1$t)
  
  #SoG_rv: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    mutate(iter = 1:n()) %>% 
    group_by(iter) %>% 
    do({ 
      # Build temp sequence 
      temp_seq <- seq(min(SoG_rv$temp), max(SoG_rv$temp), length.out = 100) 
      # Extract parameter values for this iteration 
      params <- as.list(.[param_names][1, ]) 
      # Build argument list: temp plus parameters 
      args <- c(list(temp = temp_seq), params) 
      # Call chosen_fun with correct arguments 
      data.frame(temp = temp_seq, 
                 pred = do.call(chosen_fun, args)) 
    }) %>% 
    ungroup()
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_SoG <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot bootstrapped CIs
  p = ggplot() +
    geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
    theme_bw(base_size = 17) +
    labs(x = 'Temperature (ºC)',
         y = 'Tissue weight growth (g)',
         title = paste0('TPC during period ', unique(TiW_rates_avg_with_temps$period)[i]))
  
  # plot --------
  p = ggplot() +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_CC, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_SoG, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
    scale_colour_manual(values = c("skyblue", "coral")) +
    scale_fill_manual(values = c("skyblue3", "coral3")) +
    labs(x = 'Temperature (ºC)',
         y = if (i == 1) "Tissue weight growth (g)" else NULL,
         col = "Source Region",
         fill = "Source Region",
         title = paste0('TPC during period ', i)) + 
    theme_cowplot(16) + 
    scale_x_continuous(breaks = c(10,12,14,16,18,20,22)) +
    expand_limits(x = c(10,22.5)) +
    scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5)) +
    expand_limits(y = c(-0.2,0.5)) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.text  = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  # Store the plot in the list
  plots_TiW[[i]] <- p
  
  #SoG_rv: quadratic: Estimate parameters & CI intervals 
  extra_params <- calc_params(fit_nlsLM) %>%
    pivot_longer(everything(), names_to =  'param', values_to = 'estimate')
  
  ci_extra_params <- Boot(fit_nlsLM, f = function(x){unlist(calc_params(x))}, labels = names(calc_params(fit_nlsLM)), R = 200, method = 'case') %>%
    confint(., method = 'bca') %>%
    as.data.frame() %>%
    rename(conf_lower = 1, conf_upper = 2) %>%
    rownames_to_column(., var = 'param') %>%
    mutate(method = 'case bootstrap')
  
  ci_extra_params <- left_join(ci_extra_params, extra_params)
  
  ci_params_select_SoG_fr <- ci_extra_params %>%
    filter(param == "ctmax" | param == "topt") %>%
    mutate(SR = "Strait of Georgia",
           RV = "TiW",
           model = best_model)
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                        op_sev = ci_params_select_SoG_fr$estimate[1], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                        SR = ci_params_select_SoG_fr$SR[1],
                                                        model = ci_params_select_SoG_fr$model[1],
                                                        period = i, RV = "TiW"))
  
  topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                        op_sev = ci_params_select_SoG_fr$estimate[2], 
                                                        conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                        conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                        SR = ci_params_select_SoG_fr$SR[2],
                                                        model = ci_params_select_SoG_fr$model[2],
                                                        period = i, RV = "TiW"))
}

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/Figures")
ggsave(filename="TPCs_temp.png", height=20, width=15, 
       plot=grid.arrange(plots_SL[[1]],plots_SL[[2]],plots_SL[[3]],
                         plots_SG[[1]],plots_SG[[2]],plots_SG[[3]],
                         plots_ShW[[1]],plots_ShW[[2]],plots_ShW[[3]],
                         plots_TiW[[1]],plots_TiW[[2]],plots_TiW[[3]],
                         ncol = 3))

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/")
# write.csv(topt_sev_results, "Topt_and_CTmax.csv")
topt_sev_results = read.csv("Topt_and_CTmax.csv")

# change a few things here cos it was modified to add cumulative 
colnames(topt_sev_results)[3] = "temp"
topt_sev_results = topt_sev_results[1:48,]

# Central Coast CTmax changes ----------

get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

# Shell length
L_rates_avg_with_temps_CC = L_rates_avg_with_temps[which(L_rates_avg_with_temps$SR == "Central Coast"),]

# get modes for the different periods
get_mode(L_rates_avg_with_temps_CC$starting_date[which(L_rates_avg_with_temps_CC$period == 1)]) # 13-Jun-25
get_mode(L_rates_avg_with_temps_CC$starting_date[which(L_rates_avg_with_temps_CC$period == 2)]) # 29-Jun-25
get_mode(L_rates_avg_with_temps_CC$starting_date[which(L_rates_avg_with_temps_CC$period == 3)]) # 12-Jul-25

get_mode(L_rates_avg_with_temps_CC$ending_date[which(L_rates_avg_with_temps_CC$period == 1)]) # 29-Jun-25
get_mode(L_rates_avg_with_temps_CC$ending_date[which(L_rates_avg_with_temps_CC$period == 2)]) # 12-Jul-25
get_mode(L_rates_avg_with_temps_CC$ending_date[which(L_rates_avg_with_temps_CC$period == 3)]) # 26-Jul-25

ctmax_results_CC_L = topt_sev_results[which(topt_sev_results$RV == "SL" &
                                              topt_sev_results$param == "ctmax" &
                                              topt_sev_results$SR == "Central Coast"),]
ctmax_results_CC_L$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_CC_L$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_CC_L$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_CC_L$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_CC_L, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

# add seasonally varying climatology to the plot
library(ncdf4)
library(heatwaveR)
library(patchwork) # To display 2 charts together
library(hrbrthemes)
OISST_nc = nc_open("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/oisst_of_snails_beaty.nc")
oisst_lat = ncvar_get(OISST_nc, "matched_lat")
oisst_lon = ncvar_get(OISST_nc, "matched_lon")
sst = ncvar_get(OISST_nc, "sst")
# time is wrongly saved. Let's get it back
start_date <- as.Date("1983-01-01")
end_date <- as.Date("2024-06-19")
oisst_dates <- seq(from = start_date, to = end_date, by = "days")
all_coords = data.frame(oisst_lon, oisst_lat)
colnames(all_coords)[1] = "lon"
colnames(all_coords)[2] = "lat"
OISST_of_interest = data.frame(t = oisst_dates,
                               temp = sst[,1])

ts = ts2clm(OISST_of_interest, climatologyPeriod = c("1983-01-01", "2012-12-31"))
start_date <- as.Date("2018-06-29")
end_date <- as.Date("2018-07-26")
ts2 = ts %>%
  filter(t >= start_date & t <= end_date)

# change the dates from 2018 to 2025 in ts2 (doesn't make any difference cos it's just the clim which is the same every year anyway)
ts2$t <- seq(
  from = as.Date("2025-06-29"),
  to   = as.Date("2025-07-26"),
  by   = "day"
)

ctmax_plots_CC_L = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+11), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_CC_L$temp, x=ctmax_results_CC_L$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date") + ggtitle("Central Coast")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_CC_L$date,
  temp  = ctmax_results_CC_L$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell lip
SG_rates_avg_with_temps_CC = SG_rates_avg_with_temps[which(SG_rates_avg_with_temps$SR == "Central Coast"),]

# get modes for the different periods
get_mode(SG_rates_avg_with_temps_CC$ending_date[which(SG_rates_avg_with_temps_CC$period == 1)]) # 29-Jun-25
get_mode(SG_rates_avg_with_temps_CC$ending_date[which(SG_rates_avg_with_temps_CC$period == 2)]) # 12-Jul-25
get_mode(SG_rates_avg_with_temps_CC$ending_date[which(SG_rates_avg_with_temps_CC$period == 3)]) # 26-Jul-25

ctmax_results_CC_SG = topt_sev_results[which(topt_sev_results$RV == "SG" &
                                              topt_sev_results$param == "ctmax" &
                                              topt_sev_results$SR == "Central Coast"),]
ctmax_results_CC_SG$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_CC_SG$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_CC_SG$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_CC_SG$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_CC_SG, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ctmax_plots_CC_SG = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+11), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_CC_SG$temp, x=ctmax_results_CC_SG$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_CC_SG$date,
  temp  = ctmax_results_CC_SG$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell weight
ShW_rates_avg_with_temps_CC = ShW_rates_avg_with_temps[which(ShW_rates_avg_with_temps$SR == "Central Coast"),]

# get modes for the different periods
get_mode(ShW_rates_avg_with_temps_CC$ending_date[which(ShW_rates_avg_with_temps_CC$period == 1)]) # 29-Jun-25
get_mode(ShW_rates_avg_with_temps_CC$ending_date[which(ShW_rates_avg_with_temps_CC$period == 2)]) # 12-Jul-25
get_mode(ShW_rates_avg_with_temps_CC$ending_date[which(ShW_rates_avg_with_temps_CC$period == 3)]) # 26-Jul-25

ctmax_results_CC_ShW = topt_sev_results[which(topt_sev_results$RV == "ShW" &
                                               topt_sev_results$param == "ctmax" &
                                               topt_sev_results$SR == "Central Coast"),]
ctmax_results_CC_ShW$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_CC_ShW$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_CC_ShW$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_CC_ShW$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_CC_ShW, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ctmax_plots_CC_ShW = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+11), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_CC_ShW$temp, x=ctmax_results_CC_ShW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_CC_ShW$date,
  temp  = ctmax_results_CC_ShW$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Tissue weight
TiW_rates_avg_with_temps_CC = TiW_rates_avg_with_temps[which(TiW_rates_avg_with_temps$SR == "Central Coast"),]

# get modes for the different periods
get_mode(TiW_rates_avg_with_temps_CC$ending_date[which(TiW_rates_avg_with_temps_CC$period == 1)]) # 29-Jun-25
get_mode(TiW_rates_avg_with_temps_CC$ending_date[which(TiW_rates_avg_with_temps_CC$period == 2)]) # 12-Jul-25
get_mode(TiW_rates_avg_with_temps_CC$ending_date[which(TiW_rates_avg_with_temps_CC$period == 3)]) # 26-Jul-25

ctmax_results_CC_TiW = topt_sev_results[which(topt_sev_results$RV == "TiW" &
                                                topt_sev_results$param == "ctmax" &
                                                topt_sev_results$SR == "Central Coast"),]
ctmax_results_CC_TiW$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_CC_TiW$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_CC_TiW$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_CC_TiW$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_CC_TiW, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ctmax_plots_CC_TiW = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+11), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_CC_TiW$temp, x=ctmax_results_CC_TiW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_CC_TiW$date,
  temp  = ctmax_results_CC_TiW$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Strait of Georgia CTmax changes ----------

# Shell length
L_rates_avg_with_temps_SoG = L_rates_avg_with_temps[which(L_rates_avg_with_temps$SR == "Strait of Georgia"),]

# get modes for the different periods
get_mode(L_rates_avg_with_temps_SoG$ending_date[which(L_rates_avg_with_temps_SoG$period == 1)]) # 29-Jun-25
get_mode(L_rates_avg_with_temps_SoG$ending_date[which(L_rates_avg_with_temps_SoG$period == 2)]) # 12-Jul-25
get_mode(L_rates_avg_with_temps_SoG$ending_date[which(L_rates_avg_with_temps_SoG$period == 3)]) # 26-Jul-25

ctmax_results_SoG_L = topt_sev_results[which(topt_sev_results$RV == "SL" &
                                               topt_sev_results$param == "ctmax" &
                                               topt_sev_results$SR == "Strait of Georgia"),]
ctmax_results_SoG_L$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_SoG_L$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_SoG_L$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_SoG_L$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_SoG_L, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

# add seasonally varying climatology to the plot
library(ncdf4)
library(heatwaveR)
library(patchwork) # To display 2 charts together
library(hrbrthemes)
OISST_nc = nc_open("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/oisst_of_snails_beaty.nc")
oisst_lat = ncvar_get(OISST_nc, "matched_lat")
oisst_lon = ncvar_get(OISST_nc, "matched_lon")
sst = ncvar_get(OISST_nc, "sst")
# time is wrongly saved. Let's get it back
start_date <- as.Date("1983-01-01")
end_date <- as.Date("2024-06-19")
oisst_dates <- seq(from = start_date, to = end_date, by = "days")
all_coords = data.frame(oisst_lon, oisst_lat)
colnames(all_coords)[1] = "lon"
colnames(all_coords)[2] = "lat"
OISST_of_interest = data.frame(t = oisst_dates,
                               temp = sst[,2])

ts = ts2clm(OISST_of_interest, climatologyPeriod = c("1983-01-01", "2012-12-31"))
start_date <- as.Date("2018-06-29")
end_date <- as.Date("2018-07-26")
ts3 = ts %>%
  filter(t >= start_date & t <= end_date)

# change the dates from 2018 to 2025 in ts3 (doesn't make any difference cos it's just the clim which is the same every year anyway)
ts3$t <- seq(
  from = as.Date("2025-06-29"),
  to   = as.Date("2025-07-26"),
  by   = "day"
)

ctmax_plots_SoG_L = ctmax_plots + 
  geom_line(data = ts3, aes(x = t, y = seas+7), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_SoG_L$temp, x=ctmax_results_SoG_L$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-7, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "coral3",
                                  size = 35)) +
  labs(x = "Date") + ggtitle("Strait of Georgia")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_SoG_L$date,
  temp  = ctmax_results_SoG_L$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell lip
SG_rates_avg_with_temps_SoG = SG_rates_avg_with_temps[which(SG_rates_avg_with_temps$SR == "Strait of Georgia"),]

# get modes for the different periods
get_mode(SG_rates_avg_with_temps_SoG$ending_date[which(SG_rates_avg_with_temps_SoG$period == 1)]) # 29-Jun-25
get_mode(SG_rates_avg_with_temps_SoG$ending_date[which(SG_rates_avg_with_temps_SoG$period == 2)]) # 12-Jul-25
get_mode(SG_rates_avg_with_temps_SoG$ending_date[which(SG_rates_avg_with_temps_SoG$period == 3)]) # 26-Jul-25

ctmax_results_SoG_SG = topt_sev_results[which(topt_sev_results$RV == "SG" &
                                                topt_sev_results$param == "ctmax" &
                                                topt_sev_results$SR == "Strait of Georgia"),]
ctmax_results_SoG_SG$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_SoG_SG$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_SoG_SG$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_SoG_SG$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_SoG_SG, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ctmax_plots_SoG_SG = ctmax_plots + 
  geom_line(data = ts3, aes(x = t, y = seas+7), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_SoG_SG$temp, x=ctmax_results_SoG_SG$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-7, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "coral3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_SoG_SG$date,
  temp  = ctmax_results_SoG_SG$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell weight
ShW_rates_avg_with_temps_SoG = ShW_rates_avg_with_temps[which(ShW_rates_avg_with_temps$SR == "Strait of Georgia"),]

# get modes for the different periods
get_mode(ShW_rates_avg_with_temps_SoG$ending_date[which(ShW_rates_avg_with_temps_SoG$period == 1)]) # 29-Jun-25
get_mode(ShW_rates_avg_with_temps_SoG$ending_date[which(ShW_rates_avg_with_temps_SoG$period == 2)]) # 12-Jul-25
get_mode(ShW_rates_avg_with_temps_SoG$ending_date[which(ShW_rates_avg_with_temps_SoG$period == 3)]) # 26-Jul-25

ctmax_results_SoG_ShW = topt_sev_results[which(topt_sev_results$RV == "ShW" &
                                                 topt_sev_results$param == "ctmax" &
                                                 topt_sev_results$SR == "Strait of Georgia"),]
ctmax_results_SoG_ShW$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_SoG_ShW$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_SoG_ShW$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_SoG_ShW$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_SoG_ShW, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ctmax_plots_SoG_ShW = ctmax_plots + 
  geom_line(data = ts3, aes(x = t, y = seas+11), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_SoG_ShW$temp, x=ctmax_results_SoG_ShW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "coral3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_SoG_ShW$date,
  temp  = ctmax_results_SoG_ShW$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Tissue weight
TiW_rates_avg_with_temps_SoG = TiW_rates_avg_with_temps[which(TiW_rates_avg_with_temps$SR == "Strait of Georgia"),]

# get modes for the different periods
get_mode(TiW_rates_avg_with_temps_SoG$ending_date[which(TiW_rates_avg_with_temps_SoG$period == 1)]) # 29-Jun-25
get_mode(TiW_rates_avg_with_temps_SoG$ending_date[which(TiW_rates_avg_with_temps_SoG$period == 2)]) # 12-Jul-25
get_mode(TiW_rates_avg_with_temps_SoG$ending_date[which(TiW_rates_avg_with_temps_SoG$period == 3)]) # 26-Jul-25

ctmax_results_SoG_TiW = topt_sev_results[which(topt_sev_results$RV == "TiW" &
                                                 topt_sev_results$param == "ctmax" &
                                                 topt_sev_results$SR == "Strait of Georgia"),]
ctmax_results_SoG_TiW$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_SoG_TiW$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_SoG_TiW$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_SoG_TiW$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_SoG_TiW, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ctmax_plots_SoG_TiW = ctmax_plots + 
  geom_line(data = ts3, aes(x = t, y = seas+7), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_SoG_TiW$temp, x=ctmax_results_SoG_TiW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-7, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        plot.title = element_text(hjust = 0.5, color = "coral3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_SoG_TiW$date,
  temp  = ctmax_results_SoG_TiW$temp
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(temp ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")
ggsave(filename="CTmax_vs_seasonal_temp.png", height=20, width=15,
       plot=grid.arrange(ctmax_plots_SoG_L, ctmax_plots_CC_L,
                         ctmax_plots_SoG_SG, ctmax_plots_CC_SG,
                         ctmax_plots_SoG_ShW, ctmax_plots_CC_ShW,
                         ctmax_plots_SoG_TiW, ctmax_plots_CC_TiW, nrow=4))

# Central Coast tsm changes ----------

# load lighthouse data
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")
CC_lighthouse_biweekly_95th_temp = readRDS("CC_lighthouse_biweekly_95th_temp_growth_var.Rds")
SoG_lighthouse_biweekly_95th_temp = readRDS("SoG_lighthouse_biweekly_95th_temp_growth_var.Rds")

# Shell length
L_rates_avg_with_temps_CC = L_rates_avg_with_temps[which(L_rates_avg_with_temps$SR == "Central Coast"),]

# get modes for the different periods
get_mode(L_rates_avg_with_temps_CC$ending_date[which(L_rates_avg_with_temps_CC$period == 1)]) # 29-Jun-25
get_mode(L_rates_avg_with_temps_CC$ending_date[which(L_rates_avg_with_temps_CC$period == 2)]) # 12-Jul-25
get_mode(L_rates_avg_with_temps_CC$ending_date[which(L_rates_avg_with_temps_CC$period == 3)]) # 26-Jul-25

ctmax_results_CC_L = topt_sev_results[which(topt_sev_results$RV == "SL" &
                                              topt_sev_results$param == "ctmax" &
                                              topt_sev_results$SR == "Central Coast"),]
ctmax_results_CC_L$date = as.Date("29-Jun-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_CC_L$date[1] = as.Date("29-Jun-25", format = "%d-%b-%y")
ctmax_results_CC_L$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_CC_L$date[3] = as.Date("26-Jul-25", format = "%d-%b-%y")

CC_TSM_df_L <- ctmax_results_CC_L %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(CC_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(CC_TSM_df_L, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_CC_L = tsm_plots + 
  geom_line(data = ts2, aes(x = t, y = seas-4), color="red", size = 2) +
  geom_smooth(aes(y=CC_TSM_df_L$TSM, x=CC_TSM_df_L$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+4, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date") + ggtitle("Central Coast")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = CC_TSM_df_L$date,
  TSM  = CC_TSM_df_L$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell lip
CC_TSM_df_SG <- ctmax_results_CC_SG %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(CC_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(CC_TSM_df_SG, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_CC_SG = tsm_plots + 
  geom_line(data = ts2, aes(x = t, y = seas-4), color="red", size = 2) +
  geom_smooth(aes(y=CC_TSM_df_SG$TSM, x=CC_TSM_df_SG$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+4, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = CC_TSM_df_SG$date,
  TSM  = CC_TSM_df_SG$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell weight
CC_TSM_df_ShW <- ctmax_results_CC_ShW %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(CC_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(CC_TSM_df_ShW, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_CC_ShW = tsm_plots + 
  geom_line(data = ts2, aes(x = t, y = seas-4), color="red", size = 2) +
  geom_smooth(aes(y=CC_TSM_df_ShW$TSM, x=CC_TSM_df_ShW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+4, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = CC_TSM_df_ShW$date,
  TSM  = CC_TSM_df_ShW$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Tissue weight
CC_TSM_df_TiW <- ctmax_results_CC_TiW %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(CC_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(CC_TSM_df_TiW, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_CC_TiW = tsm_plots + 
  geom_line(data = ts2, aes(x = t, y = seas-4), color="red", size = 2) +
  geom_smooth(aes(y=CC_TSM_df_TiW$TSM, x=CC_TSM_df_TiW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+4, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = CC_TSM_df_TiW$date,
  TSM  = CC_TSM_df_TiW$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Strait of Georgia tsm changes ----------

# Shell length
SoG_TSM_df_L <- ctmax_results_SoG_L %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(SoG_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(SoG_TSM_df_L, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_SoG_L = tsm_plots + 
  geom_line(data = ts3, aes(x = t, y = seas-12), color="red", size = 2) +
  geom_smooth(aes(y=SoG_TSM_df_L$TSM, x=SoG_TSM_df_L$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+12, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "coral3",
                                  size = 35)) +
  labs(x = "Date") + ggtitle("Strait of Georgia")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = SoG_TSM_df_L$date,
  TSM  = SoG_TSM_df_L$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell lip
SoG_TSM_df_SG <- ctmax_results_SoG_SG %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(SoG_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(SoG_TSM_df_SG, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_SoG_SG = tsm_plots + 
  geom_line(data = ts3, aes(x = t, y = seas-11), color="red", size = 2) +
  geom_smooth(aes(y=SoG_TSM_df_SG$TSM, x=SoG_TSM_df_SG$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = SoG_TSM_df_SG$date,
  TSM  = SoG_TSM_df_SG$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Shell weight
SoG_TSM_df_ShW <- ctmax_results_SoG_ShW %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(SoG_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(SoG_TSM_df_ShW, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_SoG_ShW = tsm_plots + 
  geom_line(data = ts3, aes(x = t, y = seas-4), color="red", size = 2) +
  geom_smooth(aes(y=SoG_TSM_df_ShW$TSM, x=SoG_TSM_df_ShW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+4, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = SoG_TSM_df_ShW$date,
  TSM  = SoG_TSM_df_ShW$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

# Tissue weight
SoG_TSM_df_TiW <- ctmax_results_SoG_TiW %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(SoG_lighthouse_biweekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

tsm_plots = ggplot(SoG_TSM_df_TiW, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

tsm_plots_SoG_TiW = tsm_plots + 
  geom_line(data = ts3, aes(x = t, y = seas-10), color="red", size = 2) +
  geom_smooth(aes(y=SoG_TSM_df_TiW$TSM, x=SoG_TSM_df_TiW$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(TSM ~ "(°C)"), sec.axis = sec_axis( trans=~.+10, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = SoG_TSM_df_TiW$date,
  TSM  = SoG_TSM_df_TiW$TSM
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM ~ date, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["date"]]
pval  <- sm$coefficients["date", "Pr(>|t|)"]
ci    <- confint(m, "date", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))
beta1
pval

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")
ggsave(filename="TSM_vs_seasonal_temp.png", height=20, width=15,
       plot=grid.arrange(tsm_plots_SoG_L, tsm_plots_CC_L,
                         tsm_plots_SoG_SG, tsm_plots_CC_SG,
                         tsm_plots_SoG_ShW, tsm_plots_CC_ShW,
                         tsm_plots_SoG_TiW, tsm_plots_CC_TiW, nrow=4))