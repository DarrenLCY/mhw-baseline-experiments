# =============================================================================
# Survival analysis of Nucella lamellosa mortality under MHW severity treatments
# Implements: Kaplan-Meier curves + Cox Proportional Hazards (factorial)
# Reference: Covernton & Harley (DOI: 10.3354/meps13355)
#
# Data format required by survival package (one row per snail):
#   - SP       : source population (Cedar, Heron, Kwak, Pruth)
#   - Treat    : MHW severity / temperature treatment
#   - time     : days from experiment start to death (or end of experiment if survived)
#   - status   : 1 = died (event occurred), 0 = censored (survived to end of experiment)
#
# NOTE: R's survival package uses status = 1 for event (death), 0 for censored (alive).
#       This is the OPPOSITE of the JMP convention described by Chris Harley
#       (where 0 = died, 1 = survived). The biology is the same; only the coding differs.
# =============================================================================

# Load packages ---------------------------------------------------------------
pkgs <- c("tidyverse", "lubridate", "survival", "survminer", "cowplot",
          "car", "broom", "ggrepel", "coxme", "viridis")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

rm(list = ls())
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/")

# =============================================================================
# 1. LOAD AND PROCESS LOGGER DATA (unchanged from original script)
# =============================================================================

growth_base <- read.csv("growth_variables_v2.csv")

# [All logger loading and cleaning code here - unchanged from original script]
# ...
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
# 
# # ---------------- match the daily temps to the Hobday thresholds to determine the 'realised' severity threshold
# # we will match them for 3 different periods (weeks 0-2, 2-4 and 4-6), due to the variable realised temperature
# 
# #let's see how it looks like first... Run experimental_temps.R until you have the final plot (located: C:\Users\dlcyli\OneDrive\Development of thesis\Nucella experiments)
# # then combine the final plot with the realised temps
# # change the dates from 2020 to 2025...
# ts_sev11 <- ts_sev11 %>% mutate(t = update(t, year = 2025))
# ts_sev12.5 <- ts_sev12.5 %>% mutate(t = update(t, year = 2025))
# ts_sev14 <- ts_sev14 %>% mutate(t = update(t, year = 2025))
# ts_sev15.5 <- ts_sev15.5 %>% mutate(t = update(t, year = 2025))
# ts_sev17 <- ts_sev17 %>% mutate(t = update(t, year = 2025))
# ts_sev18 <- ts_sev18 %>% mutate(t = update(t, year = 2025))
# ts_sev19 <- ts_sev19 %>% mutate(t = update(t, year = 2025))
# ts_sev20 <- ts_sev20 %>% mutate(t = update(t, year = 2025))
# ts_sev20.5 <- ts_sev20.5 %>% mutate(t = update(t, year = 2025))
# ts_sev21 <- ts_sev21 %>% mutate(t = update(t, year = 2025))
# ts_sev21.5 <- ts_sev21.5 %>% mutate(t = update(t, year = 2025))
# ts_sev22 <- ts_sev22 %>% mutate(t = update(t, year = 2025))
# check_point_dates = c("2025-06-09", "2025-06-13", "2025-06-26", "2025-07-10", "2025-07-26")
# ggplot() +
#   geom_line(data = ts_sev11, aes(x = t, y = sev11_avg), colour = "red", size = 1.5) +
#   geom_line(data = T11_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev12.5, aes(x = t, y = sev12.5_avg), colour = "red", size = 1.5) +
#   geom_line(data = T12.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev14, aes(x = t, y = sev14_avg), colour = "red", size = 1.5) +
#   geom_line(data = T14_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev15.5, aes(x = t, y = sev15.5_avg), colour = "red", size = 1.5) +
#   geom_line(data = T15.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev17, aes(x = t, y = sev17_avg), colour = "red", size = 1.5) +
#   geom_line(data = T17_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev18, aes(x = t, y = sev18_avg), colour = "red", size = 1.5) +
#   geom_line(data = T18_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev19, aes(x = t, y = sev19_avg), colour = "red", size = 1.5) +
#   geom_line(data = T19_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev20, aes(x = t, y = sev20_avg), colour = "red", size = 1.5) +
#   geom_line(data = T20_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev20.5, aes(x = t, y = sev20.5_avg), colour = "red", size = 1.5) +
#   geom_line(data = T20.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev21, aes(x = t, y = sev21_avg), colour = "red", size = 1.5) +
#   geom_line(data = T21_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev21.5, aes(x = t, y = sev21.5_avg), colour = "red", size = 1.5) +
#   geom_line(data = T21.5_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_line(data = ts_sev22, aes(x = t, y = sev22_avg), colour = "red", size = 1.5) +
#   geom_line(data = T22_logger_data_daily, aes(x = date, y = daily_avg_temp), colour = "blue", size = 1.5) +
#   geom_rect(aes(xmin = as.Date("2025-06-13"), xmax = as.Date("2025-07-26"), ymin = - Inf,
#                 ymax = Inf, fill = "red"), alpha = 0.2)  +
#   geom_vline(xintercept = as.Date(check_point_dates)) +
#   labs(y = "Temperatures", x = "Time") +
#   guides(fill = 'none')  +
#   ggtitle("Averaged") +
#   theme(plot.title = element_text(hjust = 0.5), text = element_text(size=15)) +
#   coord_cartesian(ylim = c(8, 25))
# 
# sog_sev11 <- sog_sev11 %>% mutate(t = update(t, year = 2025))
# sog_sev12.5 <- sog_sev12.5 %>% mutate(t = update(t, year = 2025))
# sog_sev14 <- sog_sev14 %>% mutate(t = update(t, year = 2025))
# sog_sev15.5 <- sog_sev15.5 %>% mutate(t = update(t, year = 2025))
# sog_sev17 <- sog_sev17 %>% mutate(t = update(t, year = 2025))
# sog_sev18 <- sog_sev18 %>% mutate(t = update(t, year = 2025))
# sog_sev19 <- sog_sev19 %>% mutate(t = update(t, year = 2025))
# sog_sev20 <- sog_sev20 %>% mutate(t = update(t, year = 2025))
# sog_sev20.5 <- sog_sev20.5 %>% mutate(t = update(t, year = 2025))
# sog_sev21 <- sog_sev21 %>% mutate(t = update(t, year = 2025))
# sog_sev21.5 <- sog_sev21.5 %>% mutate(t = update(t, year = 2025))
# sog_sev22 <- sog_sev22 %>% mutate(t = update(t, year = 2025))
# 
# cc_sev11 <- cc_sev11 %>% mutate(t = update(t, year = 2025))
# cc_sev12.5 <- cc_sev12.5 %>% mutate(t = update(t, year = 2025))
# cc_sev14 <- cc_sev14 %>% mutate(t = update(t, year = 2025))
# cc_sev15.5 <- cc_sev15.5 %>% mutate(t = update(t, year = 2025))
# cc_sev17 <- cc_sev17 %>% mutate(t = update(t, year = 2025))
# cc_sev18 <- cc_sev18 %>% mutate(t = update(t, year = 2025))
# cc_sev19 <- cc_sev19 %>% mutate(t = update(t, year = 2025))
# cc_sev20 <- cc_sev20 %>% mutate(t = update(t, year = 2025))
# cc_sev20.5 <- cc_sev20.5 %>% mutate(t = update(t, year = 2025))
# cc_sev21 <- cc_sev21 %>% mutate(t = update(t, year = 2025))
# cc_sev21.5 <- cc_sev21.5 %>% mutate(t = update(t, year = 2025))
# cc_sev22 <- cc_sev22 %>% mutate(t = update(t, year = 2025))

#Estimate shell weight (ShW) & tissue weight (TiW) based on the following submerged regressions for each population where x is SW (submerged weight):
#Pruth	y = 1.5889x + 0.1392
#Kwakshua	y = 1.5958x + 0.0646
#Cedar	y = 1.61x + 0.0266
#Heron	y = 1.6104x - .1292
#Calculate TiW (tissue weight) based on Shell weight and Total weight
#Remove any rows with NAs in the L, T, SG, TiW & ShW columns (these died, but I kept them in my original datasheet)
#Separate the Treat into Temp & pH columns
# ...

# Experiment dates ------------------------------------------------------------
expt_start <- as.Date("2025-06-13")
expt_end   <- as.Date("2025-07-25")   # last full day of experiment; adjust if needed
expt_duration <- as.numeric(expt_end - expt_start)  # total days = 42

# =============================================================================
# 2. BUILD ONE-ROW-PER-SNAIL SURVIVAL DATASET
# =============================================================================
# The growth_base dataframe has multiple rows per snail (weekly measurements).
# We need exactly one row per snail, recording:
#   - their population (SP) and temperature treatment (Treat)
#   - the day they died, OR the last day of the experiment if they survived
#   - a status flag: 1 = died, 0 = censored (survived)
#
# KEY FIX from original code: the original code aggregated deaths into cumulative
# counts per population, losing individual snail resolution. Survival analysis
# requires one row per individual.

# Step 2a: identify unique snails and their fate ------------------------------
# Assumes growth_base has columns: ID (or similar unique ID), SP, Treat,
# Died ("1" = died this week, "" = alive, "N/A" = already dead), Date.dead


survival_data <- growth_base %>%

  # Convert date columns
  mutate(
    Date.dead = as.Date(Date.dead, format = "%d-%b-%y"),
    Treat     = as.numeric(as.character(Treat))   # ensure numeric for Cox model
  ) %>%

  # Summarise to one row per snail: take the last recorded observation
  group_by(ID, SP, Treat) %>%      
  summarise(
    # If the snail ever has Died == "1", it died; use its Date.dead
    died      = any(Died == "1", na.rm = TRUE),
    Date.dead = first(Date.dead[Died == "1"]),  # date of death (NA if survived)
    .groups   = "drop"
  ) %>%

  # Compute survival time in days from experiment start
  mutate(
    # If alive, time = total experiment duration (censored observation)
    # If dead,  time = days from start to death
    time = if_else(
      died,
      as.numeric(Date.dead - expt_start),
      expt_duration
    ),

    # status: 1 = event (died), 0 = censored (survived to end of experiment)
    # NOTE: this is the R convention; in JMP the coding is reversed (0 = died, 1 = censored)
    status = if_else(died, 1L, 0L),

    # Make SP and Treat factors with meaningful reference levels
    SP    = factor(SP, levels = c("Heron", "Cedar", "Kwak", "Pruth")),
    Treat = as.numeric(Treat)
  ) %>%

  # Remove any snails with impossible times (< 0 or > experiment duration)
  filter(time >= 0, time <= expt_duration)

# Quick sanity check
cat("Total snails:", nrow(survival_data), "\n")
cat("Deaths by population:\n")
print(table(survival_data$SP, survival_data$status))
cat("\nCensoring rate:", round(mean(survival_data$status == 0) * 100, 1), "%\n")

# =============================================================================
# 3. KAPLAN-MEIER SURVIVAL CURVES
# =============================================================================

# --- 3a. KM by population only -----------------------------------------------
km_fit_SP <- survfit(Surv(time, status) ~ SP, data = survival_data)

p_km_SP <- ggsurvplot(
  km_fit_SP,
  data            = survival_data,
  conf.int        = TRUE,
  pval.method     = TRUE,
  risk.table      = TRUE,          # number at risk table below
  risk.table.col  = "strata",
  risk.table.height = 0.25,
  palette         = c("blue", "red", "orange", "green"),
  legend.labs     = levels(survival_data$SP),
  legend.title    = "Source population",
  xlab            = "Days from experiment start",
  ylab            = "Survival probability",
  ggtheme         = theme_bw(base_size = 24),
  surv.median.line = "hv",          # show median survival time
  ylim            = c(0.6, 1)
)
print(p_km_SP)

output_dir <- "C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/"

# Save KM plots
ggsave(filename = paste0(output_dir, "KM_by_population.png"),
       plot     = p_km_SP$plot, height = 8, width = 10)


# --- 3b. KM by MHW severity (temperature treatment) -------------------------
# Treat as factor for KM plot; continuous for Cox model (see below)
survival_data_factor <- survival_data %>%
  mutate(Treat_f = factor(Treat))

km_fit_treat <- survfit(Surv(time, status) ~ Treat_f, data = survival_data_factor)

p_km_treat <- ggsurvplot(
  km_fit_treat,
  data         = survival_data_factor,
  color        = "Treat_f",
  palette      = viridis(12, option = "turbo"),
  conf.int     = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.30,
  xlab              = "Days from experiment start",
  ylab              = "Survival probability",
  ggtheme           = theme_bw(base_size = 24),
  legend.title      = "Treatment (°C)"
)

print(p_km_treat)

# Overall log-rank test: does survival differ among treatments? ------
logrank_Treat <- survdiff(Surv(time, status) ~ Treat_f, data = survival_data_factor)
cat("\n=== Log-rank test: survival by treatment ===\n")
print(logrank_Treat)

ggsave(filename = paste0(output_dir, "KM_by_treatment.png"),
       plot     = p_km_treat$plot, height = 8, width = 10)

# separate the plots ------------ 
library(purrr)
library(patchwork)

# ensure Treat is a factor
survival_data_factor <- survival_data %>%
  mutate(Treat_f = factor(Treat))

# create ggplot objects
plots <- survival_data_factor %>%
  split(.$Treat_f) %>%
  map(~ {
    
    fit <- survfit(Surv(time, status) ~ 1, data = .x)
    
    p <- ggsurvplot(
      fit,
      data = .x,
      conf.int = TRUE,
      pval = FALSE,
      risk.table = FALSE,
      legend = "none",
      xlab = "Days from experiment start",
      ylab = "Survival probability",
      title = NULL,
      ggtheme = theme_bw(base_size = 20)
    )
    
    # extract ggplot object and add plotmath title
    p$plot +
      ggtitle(
        bquote(T[j] == .(as.character(unique(.x$Treat_f))))
      ) +
      theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title = element_blank(),
        axis.text  = element_text(size = 20)
      )
  })

# combine into 3 columns x 4 rows
library(cowplot)

combined_plot <- wrap_plots(plots, ncol = 3, nrow = 4) &
  theme(
    plot.margin = margin(5, 5, 20, 20)  # top, right, bottom, left
  )

km_by_treat = ggdraw() +
  draw_plot(combined_plot) +
  draw_label("Days from experiment start", x = 0.5, y = 0.02, vjust = 0, size = 20) +
  draw_label("Survival probability", x = 0.02, y = 0.5, angle = 90, vjust = 1, size = 20)

output_dir <- "C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/"

# Save KM plots
ggsave(filename = paste0(output_dir, "KM_by_treat.png"),
       plot     = km_by_treat, height = 13, width = 12)

# --- 3c. KM faceted by population, coloured by severity ----------------------
# This is the most informative plot, comparable to Fig. 11 in Covernton & Harley
levels(survival_data_factor$SP)[levels(survival_data_factor$SP) == "Kwak"] <- "Kwakshua"
km_fit_combined <- survfit(
  Surv(time, status) ~ SP + Treat_f,
  data = survival_data_factor
)

p_km_combined <- ggsurvplot_facet(
  km_fit_combined,
  data         = survival_data_factor,
  facet.by     = "SP",
  color        = "Treat_f",
  palette      = viridis(12, option = "turbo"),
  conf.int     = FALSE,
  pval         = FALSE,
  xlab         = "Days from experiment start",
  ylab         = "Survival probability",
  legend.title = "Treatment (°C)",
  ggtheme      = theme_bw(base_size = 24),
  short.panel.labs = TRUE
)

p_km_combined

# Save KM plots
ggsave(filename = paste0(output_dir, "KM_by_treat_SP.png"),
       plot     = p_km_combined, height = 13, width = 12)

# =============================================================================
# 4. LOG-RANK TESTS
# =============================================================================

# --- 4a. Overall log-rank test: does survival differ among populations? ------
logrank_SP <- survdiff(Surv(time, status) ~ SP, data = survival_data)
cat("\n=== Log-rank test: survival by population ===\n")
print(logrank_SP)

# --- 4b. Pairwise log-rank tests with Bonferroni correction ------------------
pairwise_SP <- pairwise_survdiff(
  Surv(time, status) ~ SP,
  data   = survival_data,
  p.adjust.method = "bonferroni"
)
cat("\n=== Pairwise log-rank tests (Bonferroni correction) ===\n")
print(pairwise_SP)

# =============================================================================
# 5. COX PROPORTIONAL HAZARDS MODEL — FACTORIAL ANALYSIS
# =============================================================================
# This is the analysis Chris described as the "factorial analysis".
# The Cox PH model estimates hazard ratios (instantaneous risk of death) as a
# function of population, MHW severity (continuous temperature), and their
# interaction.
#
# Treat (temperature) is kept CONTINUOUS here because:
#   (a) it is a continuous variable with 12 ordered levels
#   (b) a continuous term is more parsimonious and directly yields a slope
#       (hazard ratio per 1°C increase)
#   (c) Chris asks about a population × MHW severity interaction,
#       which is most naturally expressed as a continuous × factor interaction
#
# If you prefer to treat Treat as a factor, replace `Treat` with `Treat_f`
# throughout Section 5. This will give separate hazard ratios per level but
# will consume more degrees of freedom.

# --- 5a. Fit factorial Cox PH model ------------------------------------------
cox_factorial <- coxph(
  Surv(time, status) ~ SP * Treat,
  data    = survival_data,
  ties    = "efron"      # recommended for continuous time data
)

cat("\n=== Cox Proportional Hazards — factorial (SP × Treat) ===\n")
summary(cox_factorial)

# Tidy output for reporting
cox_tidy <- tidy(cox_factorial, exponentiate = TRUE, conf.int = TRUE)
cat("\n=== Hazard ratios (exponentiated coefficients) ===\n")
print(cox_tidy)

# --- 5b. Likelihood-ratio test: is the interaction term significant? ---------
cox_main_only <- coxph(
  Surv(time, status) ~ SP + Treat,
  data = survival_data,
  ties = "efron"
)

lrt <- anova(cox_main_only, cox_factorial, test = "LRT")
cat("\n=== Likelihood-ratio test: main effects vs. factorial model ===\n")
print(lrt)

# If LRT is not significant, the interaction is not supported and the simpler
# main-effects model is preferred. Report both.

# --- 5c. Check proportional hazards assumption (Schoenfeld residuals) --------
ph_test <- cox.zph(cox_main_only)
cat("\n=== Proportional hazards assumption test (Schoenfeld residuals) ===\n")
print(ph_test)

# Plot Schoenfeld residuals
par(mfrow = c(2, 3))
plot(ph_test)
par(mfrow = c(1, 1))
# If any covariate has p < 0.05 here, the PH assumption may be violated for
# that covariate. In that case consider a stratified Cox model or time-varying
# coefficients.

# --- 5d. Forest plot of hazard ratios ----------------------------------------
library(broom)
library(ggplot2)

cox_tidy <- tidy(cox_main_only, exponentiate = TRUE, conf.int = TRUE) %>%
  mutate(
    term = case_when(
      term == "Treat" ~ "Treat",
      term == "SPCedar" ~ "Cedar",
      term == "SPKwak" ~ "Kwak",
      term == "SPPruth" ~ "Pruth",
      TRUE ~ term
    )
  )

forest_plot = ggplot(cox_tidy, aes(x = estimate, y = term)) +
  geom_point() +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_bw(base_size = 24) + 
  labs(
    x = "Hazard ratio (log scale)",
    y = ""
  ) +
  scale_x_log10()

library(cowplot)

# Add a, b, c, d labels to the 4 facets of the KM plot
# Adjust x/y positions to match your 2x2 facet layout
p_km_labeled <- ggdraw(p_km_combined) +
  draw_label("a", x = 0.12, y = 0.78, fontface = "bold", size = 38) +
  draw_label("b", x = 0.57, y = 0.78, fontface = "bold", size = 38) +
  draw_label("c", x = 0.12, y = 0.42, fontface = "bold", size = 38) +
  draw_label("d", x = 0.57, y = 0.42, fontface = "bold", size = 38)

combined_plot <- plot_grid(
  p_km_labeled,
  forest_plot,
  ncol = 1,
  rel_heights = c(3, 1),
  labels = c("", "e"),   # e labels the forest plot; KM labels already drawn above
  label_size = 38,
  label_x = 0,
  label_y = 0.99
)

ggsave(
  filename = paste0(output_dir, "KM_and_forest_plots.png"),
  plot     = combined_plot,
  height   = 13,
  width    = 12
)


# =============================================================================
# 6. PREDICTED SURVIVAL CURVES FROM COX MODEL
# =============================================================================
# Plot model-predicted survival curves at representative temperatures,
# separately for each population — this shows the interaction graphically.

# Create a grid of covariate values to predict over
new_data <- expand.grid(
  SP    = levels(survival_data$SP),
  Treat = c(14, 17, 19, 21, 22)   # representative severity levels; adjust as needed
)

cox_main_only <- coxph(
  Surv(time, status) ~ SP + Treat,
  data    = survival_data,
  ties    = "efron"      # recommended for continuous time data
)

cox_pred <- survfit(cox_main_only, newdata = new_data)

p_cox_pred <- ggsurvplot(
  cox_pred,
  data         = new_data,
  conf.int     = TRUE,
  legend.title = "Population × Temperature",
  xlab         = "Days from experiment start",
  ylab         = "Predicted survival probability",
  title        = "Cox PH model — predicted survival by population and MHW severity",
  ggtheme      = theme_bw(base_size = 12)
)
print(p_cox_pred)

# =============================================================================
# 8. REPORTING GUIDANCE
# =============================================================================
# In your methods section:
#   "Survival was analysed using a Cox proportional hazards model (Cox, 1972)
#    implemented in R (v4.x) using the survival package (Therneau, 2023). The
#    model included source population (SP), MHW severity (continuous temperature,
#    °C), and their interaction as fixed effects. The proportional hazards
#    assumption was verified using Schoenfeld residuals (cox.zph). Kaplan–Meier
#    curves were plotted for visualisation and pairwise log-rank tests with
#    Bonferroni correction were used to assess pairwise population differences.
#    Snails alive at the end of the experiment were treated as right-censored
#    observations."
#
# Hazard ratio interpretation:
#   HR > 1 → higher risk of death relative to reference level / per unit increase
#   HR < 1 → lower (protective) risk
#   The interaction term SP × Treat tells you whether the slope of mortality
#   risk with temperature differs between populations (i.e., whether warm-adapted
#   populations are more or less sensitive to rising temperature than cold-adapted
#   populations).

# =============================================================================
# Hazard Ratio vs Temperature plot from a mixed-effects Cox PH model
# Y axis : Hazard ratio (log scale)
# X axis : Temperature treatment (continuous, °C)
# Lines  : One per source population (Cedar, Heron, Kwak, Pruth)
# Model  : coxme() with individual snail as random effect
# =============================================================================

library(tidyverse)
library(survival)
library(coxme)       # install.packages("coxme") if needed
library(cowplot)

# --- Assumes survival_data is already prepared (one row per snail) -----------
# Required columns:
#   time     : days to death or end of experiment
#   status   : 1 = died, 0 = censored
#   SP       : factor, source population (Heron as reference level)
#   Treat    : numeric, temperature treatment (°C)
#   ID : unique snail identifier (used as random effect)

# Ensure factor levels are set correctly (Heron = reference)
survival_data <- survival_data %>%
  mutate(
    SP    = factor(SP, levels = c("Heron", "Cedar", "Kwak", "Pruth")),
    Treat = as.numeric(as.character(Treat))
  )

# =============================================================================
# 1. FIT MIXED-EFFECTS COX PH MODEL
# =============================================================================
# Random effect: individual snail (random intercept)
# This accounts for unmeasured individual-level heterogeneity in survival
# beyond the fixed effects of population and temperature.
# Fixed effects: SP (population), Treat (temperature), and their interaction

cox_mixed <- coxme(
  Surv(time, status) ~ SP * Treat + (1 | ID),
  data = survival_data
)
cat("=== Mixed-effects Cox PH model summary ===\n")
print(summary(cox_mixed))

# =============================================================================
# 2. EXTRACT COEFFICIENTS AND VARIANCE-COVARIANCE MATRIX
# =============================================================================

coefs <- fixef(cox_mixed)        # named vector of fixed-effect coefficients
V     <- vcov(cox_mixed)         # variance-covariance matrix of fixed effects

cat("\nFixed-effect coefficient names:\n")
print(names(coefs))
# Typical names (Heron = reference):
#   "Treat", "SPCedar", "SPKwak", "SPPruth",
#   "SPCedar:Treat", "SPKwak:Treat", "SPPruth:Treat"

# =============================================================================
# 3. COMPUTE PREDICTED LOG-HR AND 95% CI ACROSS TEMPERATURE RANGE
# =============================================================================

# Temperature grid for prediction
treat_seq <- seq(
  from       = min(survival_data$Treat, na.rm = TRUE),
  to         = max(survival_data$Treat, na.rm = TRUE),
  length.out = 300
)

ref_pop <- levels(survival_data$SP)[1]   # "Heron"
pops    <- levels(survival_data$SP)       # all four populations

# Helper: build a named contrast vector for a given population at temperature T
# The contrast vector picks out the relevant linear combination of coefficients
make_contrast <- function(pop, T) {
  contrast <- setNames(rep(0, length(coefs)), names(coefs))
  
  # Temperature main effect (applies to all populations)
  contrast["Treat"] <- T
  
  # For non-reference populations, add population main effect and interaction
  if (pop != ref_pop) {
    sp_main <- paste0("SP", pop)            # e.g. "SPCedar"
    sp_int  <- paste0("SP", pop, ":Treat")  # e.g. "SPCedar:Treat"
    
    if (sp_main %in% names(contrast)) contrast[sp_main] <- 1
    if (sp_int  %in% names(contrast)) contrast[sp_int]  <- T
  }
  
  contrast
}

# Compute log HR, SE, and HR for every pop x temperature combination
pred_df <- crossing(SP = pops, Treat = treat_seq) %>%
  rowwise() %>%
  mutate(
    contrast   = list(make_contrast(SP, Treat)),
    log_hr     = as.numeric(unlist(contrast) %*% coefs),
    var_log_hr = as.numeric(t(unlist(contrast)) %*% V %*% unlist(contrast)),
    se_log_hr  = sqrt(pmax(var_log_hr, 0))   # pmax guards against tiny negatives
  ) %>%
  ungroup() %>%
  select(-contrast)

# Anchor all curves to HR = 1 at the minimum temperature:
# This makes the plot show "how much does mortality risk increase relative to
# the lowest treatment?" rather than relative to an arbitrary baseline.
ref_temp <- min(survival_data$Treat, na.rm = TRUE)

pred_df <- pred_df %>%
  group_by(SP) %>%
  mutate(
    log_hr_ref = log_hr[which.min(abs(Treat - ref_temp))],   # log HR at min temp
    log_hr_adj = log_hr - log_hr_ref,                         # centred log HR
    hr         = exp(log_hr_adj),
    hr_lower   = exp(log_hr_adj - 1.96 * se_log_hr),
    hr_upper   = exp(log_hr_adj + 1.96 * se_log_hr)
  ) %>%
  ungroup() %>%
  mutate(SP = factor(SP, levels = pops))

# =============================================================================
# 4. PLOT
# =============================================================================

custom_colors <- c(
  "Heron" = "blue",
  "Cedar" = "red",
  "Kwak"  = "orange",
  "Pruth" = "green4"
)

p_hr <- ggplot(pred_df, aes(x = Treat, y = hr, colour = SP, fill = SP)) +
  
  # 95% CI ribbon
  geom_ribbon(
    aes(ymin = hr_lower, ymax = hr_upper),
    alpha  = 0.12,
    colour = NA
  ) +
  
  # Predicted HR lines
  geom_line(linewidth = 1.3) +
  
  # Reference line at HR = 1
  geom_hline(
    yintercept = 1,
    linetype   = "dashed",
    colour     = "grey40",
    linewidth  = 0.7
  ) +
  
  # Log y-axis so proportional changes are visually equal
  scale_y_log10(
    breaks = c(0.5, 1, 2, 5, 10, 20, 50),
    labels = c("0.5", "1", "2", "5", "10", "20", "50")
  ) +
  
  scale_colour_manual(values = custom_colors, name = "Source population") +
  scale_fill_manual(  values = custom_colors, name = "Source population") +
  
  labs(
    x        = "Temperature treatment (\u00b0C)",
    y        = "Hazard ratio (log scale)",
    title    = "Mortality hazard ratio by population and MHW severity",
    subtitle = paste0(
      "Mixed-effects Cox PH model (individual random effect); ",
      "HR anchored to 1 at ", ref_temp, " \u00b0C"
    ),
    caption  = "Shading: 95% CI (delta method)"
  ) +
  
  theme_bw(base_size = 15) +
  theme(
    legend.position  = "right",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 11, colour = "grey40")
  )

print(p_hr)

# =============================================================================
# 5. SAVE
# =============================================================================

output_dir <- "C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/"

ggsave(
  filename = paste0(output_dir, "cox_HR_vs_temperature.png"),
  plot     = p_hr,
  height   = 7,
  width    = 9
)

cat("\nPlot saved to:", paste0(output_dir, "cox_HR_vs_temperature.png"), "\n")

# =============================================================================
# NOTES
# =============================================================================
# 1. The y-axis uses a log scale so that equal visual distances represent equal
#    multiplicative changes in risk. HR = 4 means four times the instantaneous
#    risk of death compared to the same population at the lowest temperature.
#
# 2. Lines diverging with increasing temperature = significant population x
#    temperature interaction (some populations more sensitive than others).
#    Parallel lines = population main effect only, no interaction.
#
# 3. The random effect (1 | ID) accounts for unmeasured individual
#    heterogeneity. Replace ID with a tank ID if tank is your cluster.
#
# 4. If coxme is unavailable or the random-effect variance is near zero, use
#    coxph() with robust standard errors as a simpler alternative:
#      coxph(Surv(time, status) ~ SP * Treat + cluster(ID), data = survival_data)
#    Then replace fixef() with coef() and vcov() works the same way.
