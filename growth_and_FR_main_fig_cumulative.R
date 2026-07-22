# Script to plot main figure (for growth variables and FR)

# Using just quadratic for several reasons:
# Beaty et al. 2023 showed that quadratic models are overall the best fitting models

# gaussian and beta rely heavily on the constant value added to the RV to transform the 
# data, such that negative values become positive. Depending on the constant value 
# added, the model fits could look very different, so there's some degree of 
# subjectivity here... 

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

# ---------------- match the daily temps to the Hobday thresholds to determine the 'realised' Severity level
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

# from here, you can run all the code until feeding rate, where you need to load additional variables

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

L_rates_avg_with_temps = data.frame(temp = vector(), L_growth = vector(), L_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
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
L_rates_avg_with_temps = rbind(L_rates_avg_with_temps, data.frame(temp = sev_thresh_by_treat$avg_temp, 
                                                                  L_growth = growth_rates_avg$meanL_growth_cum, L_growth_std_dev = growth_rates_avg$sdL_cum, 
                                                                  Treat = sev_thresh_by_treat$Treat,
                                                                  SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                  starting_date = sev_thresh_by_treat$start_date, ending_date = sev_thresh_by_treat$end_date, 
                                                                  difference_days = sev_thresh_by_treat$diff_days))

L_rates_avg_with_temps$L_growth_standardised = L_rates_avg_with_temps$L_growth/L_rates_avg_with_temps$difference_days*42
L_rates_avg_with_temps$L_growth_sd_standardised = L_rates_avg_with_temps$L_growth_std_dev/L_rates_avg_with_temps$difference_days*42

SG_rates_avg_with_temps = data.frame(temp = vector(), SG_growth = vector(), SG_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
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
SG_rates_avg_with_temps = rbind(SG_rates_avg_with_temps, data.frame(temp = sev_thresh_by_treat$avg_temp, 
                                                                    SG_growth = growth_rates_avg$meanSG_growth_cum, SG_growth_std_dev = growth_rates_avg$sdSG_cum, 
                                                                    Treat = sev_thresh_by_treat$Treat,
                                                                    SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                    starting_date = sev_thresh_by_treat$start_date, ending_date = sev_thresh_by_treat$end_date, 
                                                                    difference_days = sev_thresh_by_treat$diff_days))

SG_rates_avg_with_temps$SG_growth_standardised = SG_rates_avg_with_temps$SG_growth/SG_rates_avg_with_temps$difference_days*42
SG_rates_avg_with_temps$SG_growth_sd_standardised = SG_rates_avg_with_temps$SG_growth_std_dev/SG_rates_avg_with_temps$difference_days*42

ShW_rates_avg_with_temps = data.frame(temp = vector(), ShW_growth = vector(), ShW_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
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
ShW_rates_avg_with_temps = rbind(ShW_rates_avg_with_temps, data.frame(temp = sev_thresh_by_treat$avg_temp, 
                                                                      ShW_growth = growth_rates_avg$meanShW_growth_cum, ShW_growth_std_dev = growth_rates_avg$sdShW_cum, 
                                                                      Treat = sev_thresh_by_treat$Treat,
                                                                      SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                      starting_date = sev_thresh_by_treat$start_date, ending_date = sev_thresh_by_treat$end_date, 
                                                                      difference_days = sev_thresh_by_treat$diff_days))

ShW_rates_avg_with_temps$ShW_growth_standardised = ShW_rates_avg_with_temps$ShW_growth/ShW_rates_avg_with_temps$difference_days*42
ShW_rates_avg_with_temps$ShW_growth_sd_standardised = ShW_rates_avg_with_temps$ShW_growth_std_dev/ShW_rates_avg_with_temps$difference_days*42

TiW_rates_avg_with_temps = data.frame(temp = vector(), TiW_growth = vector(), TiW_growth_std_dev = vector(), Treat = vector(), SR = vector(), sev_thresh = vector(),
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
TiW_rates_avg_with_temps = rbind(TiW_rates_avg_with_temps, data.frame(temp = sev_thresh_by_treat$avg_temp, 
                                                                      TiW_growth = growth_rates_avg$meanTiW_growth_cum, TiW_growth_std_dev = growth_rates_avg$sdTiW_cum, 
                                                                      Treat = sev_thresh_by_treat$Treat,
                                                                      SR = growth_rates_avg$SR, sev_thresh = all_sevs2,
                                                                      starting_date = sev_thresh_by_treat$start_date, ending_date = sev_thresh_by_treat$end_date, 
                                                                      difference_days = sev_thresh_by_treat$diff_days))

TiW_rates_avg_with_temps$TiW_growth_standardised = TiW_rates_avg_with_temps$TiW_growth/TiW_rates_avg_with_temps$difference_days*42
TiW_rates_avg_with_temps$TiW_growth_sd_standardised = TiW_rates_avg_with_temps$TiW_growth_std_dev/TiW_rates_avg_with_temps$difference_days*42

# Shell length ----------------- 

# Create an empty list to store plots
plots_SL <- list()

topt_sev_results = data.frame(param = vector(), op_sev = vector(), conf_lower = vector(), 
                              conf_upper = vector(), SR = vector(), model = vector(),
                              RV = vector())

# CC Shell length 
# if you come across an error here, close R and re-load packages
# error is due to package conflict with dplyr
CC_rv <- L_rates_avg_with_temps %>% 
  filter(SR == "Central Coast") %>% 
  select(sev_thresh, L_growth_standardised, SR) %>% 
  rename("rate" = "L_growth_standardised")
CC_rv$sev_thresh = as.numeric(as.character(CC_rv$sev_thresh))

d_fits <- CC_rv %>%
  nest(data = c(sev_thresh, rate)) %>%
  mutate(
    beta = map(data, ~{
      params <- c("a","b","c","d", "e")  # parameters in your formula
      starts <- get_start_vals(CC_rv$sev_thresh, CC_rv$rate, model_name = "beta_2012")
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
        rate ~ beta_2012(sev_thresh, a, b, c, d, e),
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
    
    gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                        data = .x,
                                        iter = c(4,4,4),
                                        start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                        start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                        lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
    
    quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                         data = .x,
                                         iter = c(4,4,4),
                                         start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                         start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                         lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), CC_rv) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Severity level',
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
cc_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), CC_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Severity level',
       y = 'Shell length growth (mm)',
       title = 'Central Coast') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(CC_rv, aes(Treat, L_growth_standardised)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
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
                                  "(temp = sev_thresh, ", 
                                  paste(spec$params, collapse = ", "), ")") )

#CC_rv: Fit data
if (best_model == "quadratic" | best_model == "gaussian") {
  d_fit <- nest(CC_rv, data = c(sev_thresh, rate)) %>% 
    mutate( fit = map(data, ~nls_multstart( 
      formula = fit_formula, data = .x, iter = spec$iter, 
      start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) - 10, 
      start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) + 10, 
      lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      supp_errors = "Y", convergence_count = FALSE )), 
      # create new temperature data 
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      # predict over that data 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
} else {
  beta_params <- c("a", "b", "c", "d", "e") 
  beta_start = get_start_vals(CC_rv$sev_thresh, CC_rv$rate, model_name = "beta_2012")
  # Force valid bounds
  lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
  uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
  starts <- beta_start[beta_params]
  lowers <- lowers[beta_params]
  uppers <- uppers[beta_params]
  
  d_fit <- nest(CC_rv, data = c(sev_thresh, rate)) %>% 
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
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
}

# unnest predictions
d_preds_CC <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Central Coast")

# plot data and predictions
# ggplot() +
#   geom_line(aes(sev_thresh, .fitted), d_preds_CC, col = 'blue') +
#   geom_point(aes(sev_thresh, rate), CC_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
#        y = 'Feeding rate',
#        title = 'Central Coast')

if (best_model == "quadratic") {
  fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = sev_thresh, a, b, c), 
                     data = CC_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     weights = rep(1, times = nrow(CC_rv)))
} else if (best_model == "beta") {
  beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
  fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = sev_thresh, a, b, c, d, e), 
                     data = CC_rv,
                     start = beta_start,
                     lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                     uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                     weights = rep(1, times = nrow(CC_rv)))
} else if (best_model == "gaussian") {
  fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = sev_thresh, rmax, topt, a), 
                     data = CC_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
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
    temp_seq <- seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100) 
    # Extract parameter values for this iteration 
    params <- as.list(.[param_names][1, ]) 
    # Build argument list: temp plus parameters 
    args <- c(list(temp = temp_seq), params) 
    # Call chosen_fun with correct arguments 
    data.frame(sev_thresh = temp_seq, 
               pred = do.call(chosen_fun, args)) 
  }) %>% 
  ungroup()

# calculate bootstrapped confidence intervals
boot1_conf_preds_CC <- group_by(boot1_preds, sev_thresh) %>%
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
                                                      RV = "SL"))

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                      op_sev = ci_params_select_CC_fr$estimate[2], 
                                                      conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                      conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                      SR = ci_params_select_CC_fr$SR[2],
                                                      model = ci_params_select_CC_fr$model[2],
                                                      RV = "SL"))


# SoG Shell length

SoG_rv <- L_rates_avg_with_temps %>% 
  filter(SR == "Strait of Georgia") %>% 
  select(sev_thresh, L_growth_standardised, SR) %>% 
  rename("rate" = "L_growth_standardised")
SoG_rv$sev_thresh = as.numeric(as.character(SoG_rv$sev_thresh))

d_fits <- SoG_rv %>%
  nest(data = c(sev_thresh, rate)) %>%
  mutate(
    beta = map(data, ~{
      params <- c("a","b","c","d", "e")  # parameters in your formula
      starts <- get_start_vals(SoG_rv$sev_thresh, SoG_rv$rate, model_name = "beta_2012")
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
        rate ~ beta_2012(sev_thresh, a, b, c, d, e),
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
    
    gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                        data = .x,
                                        iter = c(4,4,4),
                                        start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                        start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                        lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
    
    quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                         data = .x,
                                         iter = c(4,4,4),
                                         start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                         start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                         lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), SoG_rv) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Severity level',
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
SoG_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), SoG_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Severity level',
       y = 'Shell length growth (mm)',
       title = 'Strait of Georgia') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(SoG_rv, aes(Treat, L_growth_standardised)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
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
                                  "(temp = sev_thresh, ", 
                                  paste(spec$params, collapse = ", "), ")") )

#SoG_rv: Fit data
if (best_model == "quadratic" | best_model == "gaussian") {
  d_fit <- nest(SoG_rv, data = c(sev_thresh, rate)) %>% 
    mutate( fit = map(data, ~nls_multstart( 
      formula = fit_formula, data = .x, iter = spec$iter, 
      start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) - 10, 
      start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) + 10, 
      lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      supp_errors = "Y", convergence_count = FALSE )), 
      # create new temperature data 
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      # predict over that data 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
} else {
  beta_params <- c("a", "b", "c", "d", "e") 
  beta_start = get_start_vals(SoG_rv$sev_thresh, SoG_rv$rate, model_name = "beta_2012")
  # Force valid bounds
  lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
  uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
  starts <- beta_start[beta_params]
  lowers <- lowers[beta_params]
  uppers <- uppers[beta_params]
  
  d_fit <- nest(SoG_rv, data = c(sev_thresh, rate)) %>% 
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
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
}

# unnest predictions
d_preds_SoG <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Strait of Georgia")

# plot data and predictions
# ggplot() +
#   geom_line(aes(sev_thresh, .fitted), d_preds_SoG, col = 'blue') +
#   geom_point(aes(sev_thresh, rate), SoG_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
#        y = 'Feeding rate',
#        title = 'Strait of Georgia')

if (best_model == "quadratic") {
  fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = sev_thresh, a, b, c), 
                     data = SoG_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     weights = rep(1, times = nrow(SoG_rv)))
} else if (best_model == "beta") {
  beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
  fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = sev_thresh, a, b, c, d, e), 
                     data = SoG_rv,
                     start = beta_start,
                     lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                     uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                     weights = rep(1, times = nrow(SoG_rv)))
} else if (best_model == "gaussian") {
  fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = sev_thresh, rmax, topt, a), 
                     data = SoG_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
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
    temp_seq <- seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100) 
    # Extract parameter values for this iteration 
    params <- as.list(.[param_names][1, ]) 
    # Build argument list: temp plus parameters 
    args <- c(list(temp = temp_seq), params) 
    # Call chosen_fun with correct arguments 
    data.frame(sev_thresh = temp_seq, 
               pred = do.call(chosen_fun, args)) 
  }) %>% 
  ungroup()

# calculate bootstrapped confidence intervals
boot1_conf_preds_SoG <- group_by(boot1_preds, sev_thresh) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup() %>% 
  mutate(SR = "Strait of Georgia")

# plot bootstrapped CIs
p = ggplot() +
  geom_line(aes(sev_thresh, .fitted), d_preds_SoG, col = 'blue') +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
  geom_point(aes(sev_thresh, rate), SoG_rv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 17) +
  labs(x = 'Severity level',
       y = 'Shell length growth (mm)')

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

# plot --------
p = ggplot() +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_CC, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_SoG, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("darkblue", "darkred")) +
  labs(x = 'Severity level',
       y = "Units: mm",
       col = "Source Region",
       fill = "Source Region",
       title = "Shell length growth") + 
  theme_cowplot(35) + 
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3,4,5,6,7,8)) +
  expand_limits(x = c(-2.3, 9.3)) +
  scale_y_continuous(breaks = c(0,1,2,3,4,5,6,7,8)) +
  expand_limits(y = c(-1,7)) +
  geom_hline(aes(yintercept = 0), linetype = 2,
             linewidth = 2) +
  annotate("segment",
           x = ci_params_select_SoG_fr$estimate[1],
           xend = ci_params_select_SoG_fr$estimate[1],
           y = -Inf,
           yend = 2.8,
           linetype = 2,
           linewidth = 2) +
  annotate("segment",
         x = ci_params_select_CC_fr$estimate[1],
         xend = ci_params_select_CC_fr$estimate[1],
         y = -Inf,
         yend = 5.9,
         linetype = 2,
         linewidth = 2) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
p
# Store the plot in the list
plots_SL[[1]] <- p

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                      op_sev = ci_params_select_SoG_fr$estimate[1], 
                                                      conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                      conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                      SR = ci_params_select_SoG_fr$SR[1],
                                                      model = ci_params_select_SoG_fr$model[1],
                                                      RV = "SL"))

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                      op_sev = ci_params_select_SoG_fr$estimate[2], 
                                                      conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                      conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                      SR = ci_params_select_SoG_fr$SR[2],
                                                      model = ci_params_select_SoG_fr$model[2],
                                                      RV = "SL"))

# statistical significance tests
cc_se = (topt_sev_results$conf_upper[2]-topt_sev_results$conf_lower[2])/(2*1.96)
sog_se = (topt_sev_results$conf_upper[4]-topt_sev_results$conf_lower[4])/(2*1.96)
diff <- topt_sev_results$op_sev[2] - topt_sev_results$op_sev[4]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value

# statistical significance tests - alpha opt
cc_se = (topt_sev_results$conf_upper[1]-topt_sev_results$conf_lower[1])/(2*1.96)
sog_se = (topt_sev_results$conf_upper[3]-topt_sev_results$conf_lower[3])/(2*1.96)
diff <- topt_sev_results$op_sev[1] - topt_sev_results$op_sev[3]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value # 0.05817723 

# Shell weight ----------------- 
plots_ShW <- list()

# CC Shell weight 

CC_rv <- ShW_rates_avg_with_temps %>% 
  filter(SR == "Central Coast") %>% 
  select(sev_thresh, ShW_growth_standardised, SR) %>% 
  rename("rate" = "ShW_growth_standardised")
CC_rv$sev_thresh = as.numeric(as.character(CC_rv$sev_thresh))

d_fits <- CC_rv %>%
  nest(data = c(sev_thresh, rate)) %>%
  mutate(
    beta = map(data, ~{
      params <- c("a","b","c","d", "e")  # parameters in your formula
      starts <- c(a = 1, b = -2, c = 25, d = 1, e = 1)
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
        rate ~ beta_2012(sev_thresh, a, b, c, d, e),
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
    
    gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                        data = .x,
                                        iter = c(4,4,4),
                                        start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                        start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                        lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
    
    quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                         data = .x,
                                         iter = c(4,4,4),
                                         start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                         start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                         lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), CC_rv) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Severity level',
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
cc_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), CC_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Severity level',
       y = 'Shell weight growth (g)',
       title = 'Central Coast') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(CC_rv, aes(Treat, ShW_growth_standardised)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
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
                                  "(temp = sev_thresh, ", 
                                  paste(spec$params, collapse = ", "), ")") )

#CC_rv: Fit data
if (best_model == "quadratic" | best_model == "gaussian") {
  d_fit <- nest(CC_rv, data = c(sev_thresh, rate)) %>% 
    mutate( fit = map(data, ~nls_multstart( 
      formula = fit_formula, data = .x, iter = spec$iter, 
      start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) - 10, 
      start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) + 10, 
      lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      supp_errors = "Y", convergence_count = FALSE )), 
      # create new temperature data 
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      # predict over that data 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
} else {
  beta_params <- c("a", "b", "c", "d", "e") 
  beta_start = get_start_vals(CC_rv$sev_thresh, CC_rv$rate, model_name = "beta_2012")
  # Force valid bounds
  lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
  uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
  starts <- beta_start[beta_params]
  lowers <- lowers[beta_params]
  uppers <- uppers[beta_params]
  
  d_fit <- nest(CC_rv, data = c(sev_thresh, rate)) %>% 
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
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
}

# unnest predictions
d_preds_CC <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Central Coast")

# plot data and predictions
# ggplot() +
#   geom_line(aes(sev_thresh, .fitted), d_preds_CC, col = 'blue') +
#   geom_point(aes(sev_thresh, rate), CC_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
#        y = 'Feeding rate',
#        title = 'Central Coast')

if (best_model == "quadratic") {
  fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = sev_thresh, a, b, c), 
                     data = CC_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     weights = rep(1, times = nrow(CC_rv)))
} else if (best_model == "beta") {
  beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
  fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = sev_thresh, a, b, c, d, e), 
                     data = CC_rv,
                     start = beta_start,
                     lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                     uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                     weights = rep(1, times = nrow(CC_rv)))
} else if (best_model == "gaussian") {
  fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = sev_thresh, rmax, topt, a), 
                     data = CC_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
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
    temp_seq <- seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100) 
    # Extract parameter values for this iteration 
    params <- as.list(.[param_names][1, ]) 
    # Build argument list: temp plus parameters 
    args <- c(list(temp = temp_seq), params) 
    # Call chosen_fun with correct arguments 
    data.frame(sev_thresh = temp_seq, 
               pred = do.call(chosen_fun, args)) 
  }) %>% 
  ungroup()

# calculate bootstrapped confidence intervals
boot1_conf_preds_CC <- group_by(boot1_preds, sev_thresh) %>%
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
                                                      RV = "ShW"))

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                      op_sev = ci_params_select_CC_fr$estimate[2], 
                                                      conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                      conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                      SR = ci_params_select_CC_fr$SR[2],
                                                      model = ci_params_select_CC_fr$model[2],
                                                      RV = "ShW"))


# SoG Shell weight
SoG_rv <- ShW_rates_avg_with_temps %>% 
  filter(SR == "Strait of Georgia") %>% 
  select(sev_thresh, ShW_growth_standardised, SR) %>% 
  rename("rate" = "ShW_growth_standardised")
SoG_rv$sev_thresh = as.numeric(as.character(SoG_rv$sev_thresh))

d_fits <- SoG_rv %>%
  nest(data = c(sev_thresh, rate)) %>%
  mutate(
    beta = map(data, ~{
      params <- c("a","b","c","d", "e")  # parameters in your formula
      starts <- get_start_vals(SoG_rv$sev_thresh, SoG_rv$rate, model_name = "beta_2012")
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
        rate ~ beta_2012(sev_thresh, a, b, c, d, e),
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
    
    gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                        data = .x,
                                        iter = c(4,4,4),
                                        start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                        start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                        lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
    
    quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                         data = .x,
                                         iter = c(4,4,4),
                                         start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                         start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                         lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), SoG_rv) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Severity level',
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
SoG_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), SoG_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Severity level',
       y = 'Shell weight growth (g)',
       title = 'Strait of Georgia') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(SoG_rv, aes(Treat, ShW_growth_standardised)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
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
                                  "(temp = sev_thresh, ", 
                                  paste(spec$params, collapse = ", "), ")") )

#SoG_rv: Fit data
if (best_model == "quadratic" | best_model == "gaussian") {
  d_fit <- nest(SoG_rv, data = c(sev_thresh, rate)) %>% 
    mutate( fit = map(data, ~nls_multstart( 
      formula = fit_formula, data = .x, iter = spec$iter, 
      start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) - 10, 
      start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) + 10, 
      lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      supp_errors = "Y", convergence_count = FALSE )), 
      # create new temperature data 
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      # predict over that data 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
} else {
  beta_params <- c("a", "b", "c", "d", "e") 
  beta_start = get_start_vals(SoG_rv$sev_thresh, SoG_rv$rate, model_name = "beta_2012")
  # Force valid bounds
  lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
  uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
  starts <- beta_start[beta_params]
  lowers <- lowers[beta_params]
  uppers <- uppers[beta_params]
  
  d_fit <- nest(SoG_rv, data = c(sev_thresh, rate)) %>% 
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
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
}

# unnest predictions
d_preds_SoG <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Strait of Georgia")

# plot data and predictions
# ggplot() +
#   geom_line(aes(sev_thresh, .fitted), d_preds_SoG, col = 'blue') +
#   geom_point(aes(sev_thresh, rate), SoG_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
#        y = 'Feeding rate',
#        title = 'Strait of Georgia')

if (best_model == "quadratic") {
  fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = sev_thresh, a, b, c), 
                     data = SoG_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     weights = rep(1, times = nrow(SoG_rv)))
} else if (best_model == "beta") {
  beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
  fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = sev_thresh, a, b, c, d, e), 
                     data = SoG_rv,
                     start = beta_start,
                     lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                     uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                     weights = rep(1, times = nrow(SoG_rv)))
} else if (best_model == "gaussian") {
  fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = sev_thresh, rmax, topt, a), 
                     data = SoG_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
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
    temp_seq <- seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100) 
    # Extract parameter values for this iteration 
    params <- as.list(.[param_names][1, ]) 
    # Build argument list: temp plus parameters 
    args <- c(list(temp = temp_seq), params) 
    # Call chosen_fun with correct arguments 
    data.frame(sev_thresh = temp_seq, 
               pred = do.call(chosen_fun, args)) 
  }) %>% 
  ungroup()

# calculate bootstrapped confidence intervals
boot1_conf_preds_SoG <- group_by(boot1_preds, sev_thresh) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup() %>% 
  mutate(SR = "Strait of Georgia")

# plot bootstrapped CIs
p = ggplot() +
  geom_line(aes(sev_thresh, .fitted), d_preds_SoG, col = 'blue') +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
  geom_point(aes(sev_thresh, rate), SoG_rv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 17) +
  labs(x = 'Severity level',
       y = 'Shell weight growth (g)')

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

# plot --------
p = ggplot() +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_CC, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_SoG, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("darkblue", "darkred")) +
  labs(x = 'Severity level',
       y = " Units: g",
       col = "Source Region",
       fill = "Source Region",
       title = "Shell weight growth") + 
  theme_cowplot(35) + 
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3,4,5,6,7,8)) +
  expand_limits(x = c(-2.3, 7.3)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)) +
  expand_limits(y = c(-0.1,0.9)) +
  geom_hline(aes(yintercept = 0), linetype = 2,
             linewidth = 2) +
  annotate("segment",
           x = ci_params_select_SoG_fr$estimate[1],
           xend = ci_params_select_SoG_fr$estimate[1],
           y = -Inf,
           yend = 0.47,
           linetype = 2,
           linewidth = 2) +
  annotate("segment",
           x = ci_params_select_CC_fr$estimate[1],
           xend = ci_params_select_CC_fr$estimate[1],
           y = -Inf,
           yend = 0.62,
           linetype = 2,
           linewidth = 2) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
p
# Store the plot in the list
plots_ShW[[5]] <- p


topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                      op_sev = ci_params_select_SoG_fr$estimate[1], 
                                                      conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                      conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                      SR = ci_params_select_SoG_fr$SR[1],
                                                      model = ci_params_select_SoG_fr$model[1],
                                                      RV = "ShW"))

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                      op_sev = ci_params_select_SoG_fr$estimate[2], 
                                                      conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                      conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                      SR = ci_params_select_SoG_fr$SR[2],
                                                      model = ci_params_select_SoG_fr$model[2],
                                                      RV = "ShW"))

# statistical significance tests
cc_se = (topt_sev_results$conf_upper[2]-topt_sev_results$conf_lower[2])/(2*1.96)
sog_se = (topt_sev_results$conf_upper[4]-topt_sev_results$conf_lower[4])/(2*1.96)
diff <- topt_sev_results$op_sev[2] - topt_sev_results$op_sev[4]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value

# statistical significance tests - alpha opt
cc_se = (topt_sev_results$conf_upper[1]-topt_sev_results$conf_lower[1])/(2*1.96)
sog_se = (topt_sev_results$conf_upper[3]-topt_sev_results$conf_lower[3])/(2*1.96)
diff <- topt_sev_results$op_sev[1] - topt_sev_results$op_sev[3]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value # 0.02371422 

# Tissue weight only ----------------- 
plots_TiW <- list()

# CC Tissue weight 

CC_rv <- TiW_rates_avg_with_temps %>% 
  filter(SR == "Central Coast") %>% 
  select(sev_thresh, TiW_growth_standardised, SR) %>% 
  rename("rate" = "TiW_growth_standardised")
CC_rv$sev_thresh = as.numeric(as.character(CC_rv$sev_thresh))

d_fits <- CC_rv %>%
  nest(data = c(sev_thresh, rate)) %>%
  mutate(
    beta = map(data, ~{
      params <- c("a","b","c","d", "e")  # parameters in your formula
      starts <- get_start_vals(CC_rv$sev_thresh, CC_rv$rate, model_name = "beta_2012")
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
        rate ~ beta_2012(sev_thresh, a, b, c, d, e),
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
    
    gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                        data = .x,
                                        iter = c(4,4,4),
                                        start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                        start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                        lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
    
    quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                         data = .x,
                                         iter = c(4,4,4),
                                         start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                         start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                         lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), CC_rv) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Severity level',
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
cc_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), CC_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Severity level',
       y = 'Tissue weight growth (g)',
       title = 'Central Coast') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(CC_rv, aes(Treat, TiW_growth_standardised)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
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
                                  "(temp = sev_thresh, ", 
                                  paste(spec$params, collapse = ", "), ")") )

#CC_rv: Fit data
if (best_model == "quadratic" | best_model == "gaussian") {
  d_fit <- nest(CC_rv, data = c(sev_thresh, rate)) %>% 
    mutate( fit = map(data, ~nls_multstart( 
      formula = fit_formula, data = .x, iter = spec$iter, 
      start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) - 10, 
      start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) + 10, 
      lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      supp_errors = "Y", convergence_count = FALSE )), 
      # create new temperature data 
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      # predict over that data 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
} else {
  beta_params <- c("a", "b", "c", "d", "e") 
  beta_start = get_start_vals(CC_rv$sev_thresh, CC_rv$rate, model_name = "beta_2012")
  # Force valid bounds
  lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
  uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
  starts <- beta_start[beta_params]
  lowers <- lowers[beta_params]
  uppers <- uppers[beta_params]
  
  d_fit <- nest(CC_rv, data = c(sev_thresh, rate)) %>% 
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
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
}

# unnest predictions
d_preds_CC <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Central Coast")

# plot data and predictions
# ggplot() +
#   geom_line(aes(sev_thresh, .fitted), d_preds_CC, col = 'blue') +
#   geom_point(aes(sev_thresh, rate), CC_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
#        y = 'Feeding rate',
#        title = 'Central Coast')

if (best_model == "quadratic") {
  fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = sev_thresh, a, b, c), 
                     data = CC_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     weights = rep(1, times = nrow(CC_rv)))
} else if (best_model == "beta") {
  beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
  fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = sev_thresh, a, b, c, d, e), 
                     data = CC_rv,
                     start = beta_start,
                     lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                     uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                     weights = rep(1, times = nrow(CC_rv)))
} else if (best_model == "gaussian") {
  fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = sev_thresh, rmax, topt, a), 
                     data = CC_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = chosen_model_name),
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
    temp_seq <- seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100) 
    # Extract parameter values for this iteration 
    params <- as.list(.[param_names][1, ]) 
    # Build argument list: temp plus parameters 
    args <- c(list(temp = temp_seq), params) 
    # Call chosen_fun with correct arguments 
    data.frame(sev_thresh = temp_seq, 
               pred = do.call(chosen_fun, args)) 
  }) %>% 
  ungroup()

# calculate bootstrapped confidence intervals
boot1_conf_preds_CC <- group_by(boot1_preds, sev_thresh) %>%
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
                                                      RV = "TiW"))

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                      op_sev = ci_params_select_CC_fr$estimate[2], 
                                                      conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                      conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                      SR = ci_params_select_CC_fr$SR[2],
                                                      model = ci_params_select_CC_fr$model[2],
                                                      RV = "TiW"))


# SoG Tissue weight

SoG_rv <- TiW_rates_avg_with_temps %>% 
  filter(SR == "Strait of Georgia") %>% 
  select(sev_thresh, TiW_growth_standardised, SR) %>% 
  rename("rate" = "TiW_growth_standardised")
SoG_rv$sev_thresh = as.numeric(as.character(SoG_rv$sev_thresh))

d_fits <- SoG_rv %>%
  nest(data = c(sev_thresh, rate)) %>%
  mutate(
    beta = map(data, ~{
      params <- c("a","b","c","d", "e")  # parameters in your formula
      starts <- get_start_vals(SoG_rv$sev_thresh, SoG_rv$rate, model_name = "beta_2012")
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
        rate ~ beta_2012(sev_thresh, a, b, c, d, e),
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
    
    gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                        data = .x,
                                        iter = c(4,4,4),
                                        start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                        start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                        lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                        supp_errors = 'Y',
                                        convergence_count = FALSE)),
    
    quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                         data = .x,
                                         iter = c(4,4,4),
                                         start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                         start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                         lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                         supp_errors = 'Y',
                                         convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', beta:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), SoG_rv) +
#   theme_bw(base_size = 12) +
#   theme(legend.position = 'none') +
#   labs(x = 'Severity level',
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
SoG_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), SoG_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Severity level',
       y = 'Tissue weight growth (g)',
       title = 'Strait of Georgia') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(SoG_rv, aes(Treat, TiW_growth_standardised)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
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
                                  "(temp = sev_thresh, ", 
                                  paste(spec$params, collapse = ", "), ")") )

#SoG_rv: Fit data
if (best_model == "quadratic" | best_model == "gaussian") {
  d_fit <- nest(SoG_rv, data = c(sev_thresh, rate)) %>% 
    mutate( fit = map(data, ~nls_multstart( 
      formula = fit_formula, data = .x, iter = spec$iter, 
      start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) - 10, 
      start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = chosen_model_name) + 10, 
      lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = chosen_model_name), 
      supp_errors = "Y", convergence_count = FALSE )), 
      # create new temperature data 
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      # predict over that data 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
} else {
  beta_params <- c("a", "b", "c", "d", "e") 
  beta_start = get_start_vals(SoG_rv$sev_thresh, SoG_rv$rate, model_name = "beta_2012")
  # Force valid bounds
  lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100)
  uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100)
  starts <- beta_start[beta_params]
  lowers <- lowers[beta_params]
  uppers <- uppers[beta_params]
  
  d_fit <- nest(SoG_rv, data = c(sev_thresh, rate)) %>% 
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
      new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))), 
      preds = map2(fit, new_data, ~augment(.x, newdata = .y)) )
}

# unnest predictions
d_preds_SoG <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Strait of Georgia")

# plot data and predictions
# ggplot() +
#   geom_line(aes(sev_thresh, .fitted), d_preds_SoG, col = 'blue') +
#   geom_point(aes(sev_thresh, rate), SoG_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Severity level (ºC)',
#        y = 'Feeding rate',
#        title = 'Strait of Georgia')

if (best_model == "quadratic") {
  fit_nlsLM <- nlsLM(rate ~ quadratic_2008(temp = sev_thresh, a, b, c), 
                     data = SoG_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     weights = rep(1, times = nrow(SoG_rv)))
} else if (best_model == "beta") {
  beta_start = c(a = 1, b = -2, c = 25, d = 1, e = 1)
  fit_nlsLM <- nlsLM(rate ~ beta_2012(temp = sev_thresh, a, b, c, d, e), 
                     data = SoG_rv,
                     start = beta_start,
                     lowers <- c(a = -10, b = -100, c = -200, d = -100, e = -100),
                     uppers <- c(a = 50, b = 100, c = 300, d = 100, e = 100),
                     weights = rep(1, times = nrow(SoG_rv)))
} else if (best_model == "gaussian") {
  fit_nlsLM <- nlsLM(rate ~ gaussian_1987(temp = sev_thresh, rmax, topt, a), 
                     data = SoG_rv,
                     start = coef(d_fit[["fit"]][[1]]),
                     lower = get_lower_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
                     upper = get_upper_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = chosen_model_name),
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
    temp_seq <- seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100) 
    # Extract parameter values for this iteration 
    params <- as.list(.[param_names][1, ]) 
    # Build argument list: temp plus parameters 
    args <- c(list(temp = temp_seq), params) 
    # Call chosen_fun with correct arguments 
    data.frame(sev_thresh = temp_seq, 
               pred = do.call(chosen_fun, args)) 
  }) %>% 
  ungroup()

# calculate bootstrapped confidence intervals
boot1_conf_preds_SoG <- group_by(boot1_preds, sev_thresh) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup() %>% 
  mutate(SR = "Strait of Georgia")

# plot bootstrapped CIs
p = ggplot() +
  geom_line(aes(sev_thresh, .fitted), d_preds_SoG, col = 'blue') +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
  geom_point(aes(sev_thresh, rate), SoG_rv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 17) +
  labs(x = 'Severity level',
       y = 'Tissue weight growth (g)')

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

# plot --------
p = ggplot() +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_CC, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_SoG, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("darkblue", "darkred")) +
  labs(x = expression("Severity level (" * alpha[i] * ")"),
       y = " Units: g",
       col = "Source Region",
       fill = "Source Region",
       title = "Tissue weight growth") + 
  theme_cowplot(35) + 
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3,4,5,6,7,8)) +
  expand_limits(x = c(-2.3, 7.3)) +
  scale_y_continuous(breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)) +
  expand_limits(y = c(-0.15,1.2)) +
  geom_hline(aes(yintercept = 0), linetype = 2,
             linewidth = 2) +
  annotate("segment",
           x = ci_params_select_SoG_fr$estimate[1],
           xend = ci_params_select_SoG_fr$estimate[1],
           y = -Inf,
           yend = 0.36,
           linetype = 2,
           linewidth = 2) +
  annotate("segment",
           x = ci_params_select_CC_fr$estimate[1],
           xend = ci_params_select_CC_fr$estimate[1],
           y = -Inf,
           yend = 0.67,
           linetype = 2,
           linewidth = 2) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
p
# Store the plot in the list
plots_TiW[[5]] <- p

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                      op_sev = ci_params_select_SoG_fr$estimate[1], 
                                                      conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                      conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                      SR = ci_params_select_SoG_fr$SR[1],
                                                      model = ci_params_select_SoG_fr$model[1],
                                                      RV = "TiW"))

topt_sev_results = rbind(topt_sev_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                      op_sev = ci_params_select_SoG_fr$estimate[2], 
                                                      conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                      conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                      SR = ci_params_select_SoG_fr$SR[2],
                                                      model = ci_params_select_SoG_fr$model[2],
                                                      RV = "TiW"))

# statistical significance tests
cc_se = (topt_sev_results$conf_upper[2]-topt_sev_results$conf_lower[2])/(2*1.96)
sog_se = (topt_sev_results$conf_upper[4]-topt_sev_results$conf_lower[4])/(2*1.96)
diff <- topt_sev_results$op_sev[2] - topt_sev_results$op_sev[4]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value

# statistical significance tests - alpha opt
cc_se = (topt_sev_results$conf_upper[1]-topt_sev_results$conf_lower[1])/(2*1.96)
sog_se = (topt_sev_results$conf_upper[3]-topt_sev_results$conf_lower[3])/(2*1.96)
diff <- topt_sev_results$op_sev[1] - topt_sev_results$op_sev[3]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value # 0.2325454 

# feeding rate -------------
# run lines 371 to 995 in C:\Users\dlcyli\OneDrive\Development of thesis\Nucella experiments\Data\feeding_rate v2.R to load some variables
feeding_base <- read.csv("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/per_weight_feeding_rate.csv")

growth_clean <- growth_base %>% 
  rename(L = Length, SG = Linear_shell_growth, TW = Total_weight, SW = Submerged_weight) %>% 
  mutate(ShW = ifelse(SP == "Cedar", cedar_reg(SW), 
                      ifelse(SP == "Heron", heron_reg(SW),
                             ifelse(SP == "Kwak", kwak_reg(SW),
                                    ifelse(SP == "Pruth", pruth_reg(SW), NA)))),
         SR = as.factor(ifelse(SP == "Cedar" | SP == "Heron", "Strait of Georgia", "Central Coast")),
         SP = as.factor(SP),
         TiW = TW-ShW)

feeding_clean <- feeding_base %>% 
  mutate(SR = as.factor(ifelse(SP == "C" | SP == "H", "Strait of Georgia", "Central Coast")),
         SP = as.factor(SP),
         Stage = as.factor(Stage), 
         wet_weight = ifelse(Length == 0.00 , 0, -2.889 + 0.1434*Length + 0.00309*(Length-32.52)^2)) # from Gooding and Harley (2015)

feeding_clean_total <- feeding_clean %>% 
  unite(unique_ID, c(Stage, SP, Treatment), sep = "_", remove = FALSE) %>% 
  group_by(Stage, SP, Treatment, Date) %>% 
  summarize(totalW_mussels_eaten = sum(wet_weight, na.rm = TRUE), n = n()) %>% 
  ungroup()

feeding_clean_total$mean_T = NA
feeding_clean_total$days = NA
for (i in 1:length(feeding_clean_total$Stage)) {
  if (feeding_clean_total$Stage[i] == "wk2") {
    treatment_temp_of_interest = all_daily_temps[which(all_daily_temps$Treat == feeding_clean_total$Treatment[i]),]
    treatment_temp_of_interest = treatment_temp_of_interest %>% 
      filter(date >= as.Date("2025-06-20") & date <= as.Date("2025-06-27"))
    feeding_clean_total$mean_T[i] = mean(treatment_temp_of_interest$daily_avg_temp)
    feeding_clean_total$days[i] = as.Date("2025-06-27") - as.Date("2025-06-20") 
  } else{
    current_date = feeding_clean_total$Date[i]
    prev_week = paste0("wk",as.numeric(substr(feeding_clean_total$Stage[i], 3,3))-1)
    previous_date = feeding_clean_total$Date[which(feeding_clean_total$Stage == prev_week &
                                                     feeding_clean_total$SP == feeding_clean_total$SP[i] &
                                                     feeding_clean_total$Treatment == feeding_clean_total$Treatment[i])]
    if (length(previous_date) == 0) {
      # case where no record is found for wk2 T15.5 cos no length was taken for this one
      previous_date = "27-Jun-25"
    }
    treatment_temp_of_interest = all_daily_temps[which(all_daily_temps$Treat == feeding_clean_total$Treatment[i]),]
    treatment_temp_of_interest = treatment_temp_of_interest %>% 
      filter(date >= as.Date(previous_date, format = "%d-%b-%y") & date <= as.Date(current_date, format = "%d-%b-%y"))
    feeding_clean_total$mean_T[i] = mean(treatment_temp_of_interest$daily_avg_temp)
    feeding_clean_total$days[i] = as.Date(current_date, format = "%d-%b-%y") - as.Date(previous_date, format = "%d-%b-%y")
  }
}

feeding_clean_total$cc_sev_thresh = NA
feeding_clean_total$sog_sev_thresh = NA
for (i in 1:length(feeding_clean_total$cc_sev_thresh)) {
  temp_df = paste0("T", feeding_clean_total$Treatment[i], "_logger_data_weekly")
  period_of_interest = substr(feeding_clean_total$Stage[i], 3, 3)
  
  # Retrieve the object from the global environment
  temp_df <- get(temp_df, envir = .GlobalEnv)
  
  feeding_clean_total$cc_sev_thresh[i] = temp_df$cc_sev_thresh[which(temp_df$group == period_of_interest)]
  feeding_clean_total$sog_sev_thresh[i] = temp_df$sog_sev_thresh[which(temp_df$group == period_of_interest)]
}

feeding_clean_total$weekly_totalW_mussels_eaten = feeding_clean_total$totalW_mussels_eaten/feeding_clean_total$days*7
feeding_clean_total$sum_TW_snails = NA
feeding_clean_total$sum_TiW_snails = NA
for (i in 1:length(feeding_clean_total$Stage)) {
  # for weeks 2, 4 and 6, we have recorded the weights, so we can match them directly
  if (feeding_clean_total$Stage[i] == "wk2" | feeding_clean_total$Stage[i] == "wk4" | feeding_clean_total$Stage[i] == "wk6") {
    feeding_clean_total$sum_TW_snails[i] = sum(growth_base$Total_weight[which(substr(growth_base$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                                growth_base$Stage == feeding_clean_total$Stage[i] &
                                                                                growth_base$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
    feeding_clean_total$sum_TiW_snails[i] = sum(growth_clean$TiW[which(substr(growth_clean$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                         growth_clean$Stage == feeding_clean_total$Stage[i] &
                                                                         growth_clean$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
    # however, for weeks 3 and 5, we can just match the snail weights to the weeks when they were last recorded
  } else if (feeding_clean_total$Stage[i] == "wk3") {
    feeding_clean_total$sum_TW_snails[i] = sum(growth_base$Total_weight[which(substr(growth_base$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                                growth_base$Stage == "wk2" &
                                                                                growth_base$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
    feeding_clean_total$sum_TiW_snails[i] = sum(growth_clean$TiW[which(substr(growth_clean$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                         growth_clean$Stage == "wk2" &
                                                                         growth_clean$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
  } else if (feeding_clean_total$Stage[i] == "wk5") {
    feeding_clean_total$sum_TW_snails[i] = sum(growth_base$Total_weight[which(substr(growth_base$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                                growth_base$Stage == "wk4" &
                                                                                growth_base$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
    feeding_clean_total$sum_TiW_snails[i] = sum(growth_clean$TiW[which(substr(growth_clean$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                         growth_clean$Stage == "wk4" &
                                                                         growth_clean$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
  }
}

feeding_clean_total$meanW_mussels_per_snail_W = feeding_clean_total$weekly_totalW_mussels_eaten/feeding_clean_total$sum_TW_snails
feeding_clean_total$meanW_mussels_per_snail_TiW = feeding_clean_total$weekly_totalW_mussels_eaten/feeding_clean_total$sum_TiW_snails

plots <- list()

alpha_opt_max_results = data.frame(param = vector(), sev_thresh = vector(), conf_lower = vector(), 
                                   conf_upper = vector(), SR = vector(), model = vector(), 
                                   RV = vector())

for (i in 1:5) {
  
  if (exists("CC_rv")) {
    # Central Coast
    CC_rv_temp <- feeding_clean_total %>% 
      filter((SP == "K" & Stage == unique(feeding_clean_total$Stage)[i]) |
               (SP == "P" & Stage == unique(feeding_clean_total$Stage)[i])) %>% 
      select(cc_sev_thresh, meanW_mussels_per_snail_TiW) %>% 
      rename("rate" = "meanW_mussels_per_snail_TiW")
    CC_rv_temp$sev_thresh = as.numeric(as.character(CC_rv_temp$cc_sev_thresh))
    CC_rv_temp$SR = "Central Coast"
    CC_rv_temp = CC_rv_temp[,2:4]
    CC_rv = rbind(CC_rv, CC_rv_temp)
  } else {
    # Central Coast
    CC_rv <- feeding_clean_total %>% 
      filter((SP == "K" & Stage == unique(feeding_clean_total$Stage)[i]) |
               (SP == "P" & Stage == unique(feeding_clean_total$Stage)[i])) %>% 
      select(cc_sev_thresh, meanW_mussels_per_snail_TiW) %>% 
      rename("rate" = "meanW_mussels_per_snail_TiW")
    CC_rv$sev_thresh = as.numeric(as.character(CC_rv$cc_sev_thresh))
    CC_rv$SR = "Central Coast"
    CC_rv = CC_rv[,2:4]
  }
  
  if (exists("SoG_rv")) {
    # Strait of Georgia
    SoG_rv_temp <- feeding_clean_total %>% 
      filter((SP == "C" & Stage == unique(feeding_clean_total$Stage)[i]) |
               (SP == "H" & Stage == unique(feeding_clean_total$Stage)[i])) %>% 
      select(sog_sev_thresh, meanW_mussels_per_snail_TiW) %>% 
      rename("rate" = "meanW_mussels_per_snail_TiW")
    SoG_rv_temp$sev_thresh = as.numeric(as.character(SoG_rv_temp$sog_sev_thresh))
    SoG_rv_temp$SR = "Strait of Georgia"
    SoG_rv_temp = SoG_rv_temp[,2:4]
    SoG_rv = rbind(SoG_rv, SoG_rv_temp)
  } else {
    # Strait of Georgia
    SoG_rv <- feeding_clean_total %>% 
      filter((SP == "C" & Stage == unique(feeding_clean_total$Stage)[i]) |
               (SP == "H" & Stage == unique(feeding_clean_total$Stage)[i])) %>% 
      select(sog_sev_thresh, meanW_mussels_per_snail_TiW) %>% 
      rename("rate" = "meanW_mussels_per_snail_TiW")
    SoG_rv$sev_thresh = as.numeric(as.character(SoG_rv$sog_sev_thresh))
    SoG_rv$SR = "Strait of Georgia"
    SoG_rv = SoG_rv[,2:4]
  }
  
}

# fit five chosen model formulations in rTPC
d_fits <- nest(CC_rv, data = c(sev_thresh, rate)) %>%
  mutate(briere = map(data, ~nls_multstart(rate~briere2_1999(sev_thresh, tmin, tmax, a, b),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'briere2_1999') - 10,
                                           start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'briere2_1999') + 10,
                                           lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'briere2_1999'),
                                           upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'briere2_1999'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', briere:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), CC_rv) +
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

d_ic

# filter for best model
best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)

# best model is set to quadratic everywhere
best_model = "quadratic"

# get colour code
col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# plot
cc_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), CC_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Feeding rate',
       title = 'Central Coast') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(CC_rv, aes(sev_thresh, rate)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rate')

#CC_rv: Quadratic: Fit data
# fit with Gaussian model
d_fit <- nest(CC_rv, data = c(sev_thresh, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))),
         # predict over that data,
         preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds_CC <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Central Coast")

# plot data and predictions
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
#   geom_point(aes(temp, rate), CC_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (ºC)',
#        y = 'Feeding rate',
#        title = 'Central Coast')

#CC_rv: Quadratic: refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                               data = CC_rv,
                               start = coef(d_fit$quadratic[[1]]),
                               lower = get_lower_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = 'quadratic_2008'),
                               upper = get_upper_lims(CC_rv$sev_thresh, CC_rv$rate, model_name = 'quadratic_2008'),
                               weights = rep(1, times = nrow(CC_rv)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# look at the data
# head(boot1$t)

# hist(boot1, layout = c(2,2))

#CC_rv: quadratic: Now plot the bootstrapped models
#create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(sev_thresh = seq(min(CC_rv$sev_thresh), max(CC_rv$sev_thresh), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = sev_thresh, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds_CC <- group_by(boot1_preds, sev_thresh) %>%
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
         RV = "fr",
         model = "quadratic")

alpha_opt_max_results = rbind(alpha_opt_max_results, data.frame(param = ci_params_select_CC_fr$param[1],
                                                                sev_thresh = ci_params_select_CC_fr$estimate[1], 
                                                                conf_lower = ci_params_select_CC_fr$conf_lower[1], 
                                                                conf_upper = ci_params_select_CC_fr$conf_upper[1],
                                                                SR = ci_params_select_CC_fr$SR[1],
                                                                model = ci_params_select_CC_fr$model[1],
                                                                RV = "FR"))

alpha_opt_max_results = rbind(alpha_opt_max_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                                sev_thresh = ci_params_select_CC_fr$estimate[2], 
                                                                conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                                conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                                SR = ci_params_select_CC_fr$SR[2],
                                                                model = ci_params_select_CC_fr$model[2],
                                                                RV = "FR"))

# fit five chosen model formulations in rTPC
d_fits <- nest(SoG_rv, data = c(sev_thresh, rate)) %>%
  mutate(briere = map(data, ~nls_multstart(rate~briere2_1999(sev_thresh, tmin, tmax, a, b),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'briere2_1999') - 10,
                                           start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'briere2_1999') + 10,
                                           lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'briere2_1999'),
                                           upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'briere2_1999'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = sev_thresh, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', briere:quadratic)

# get predictions using augment
newdata <- tibble(sev_thresh = seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, sev_thresh < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
# ggplot(d_preds, aes(sev_thresh, .fitted)) +
#   geom_line(aes(col = model_name)) +
#   geom_label_repel(aes(sev_thresh, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
#   geom_point(aes(sev_thresh, rate), SoG_rv) +
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

d_ic

# filter for best model
best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)

# best model is set to quadratic everywhere
best_model = "quadratic"

# get colour code
col_best_mod = RColorBrewer::brewer.pal(n = 6, name = "Dark2")[6]

# plot
SoG_best_fr <- ggplot(d_preds, aes(sev_thresh, .fitted)) +
  geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
  geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
  geom_label_repel(aes(sev_thresh, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
  geom_point(aes(sev_thresh, rate), SoG_rv) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Feeding rate',
       title = 'Strait of Georgia') +
  geom_hline(aes(yintercept = 0), linetype = 2) 

#Visualize the data
# ggplot(SoG_rv, aes(Treat, L_growth)) +
#   geom_point() +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (ºC)',
#        y = 'Rate')

#SoG_rv: Quadratic: Fit data
# fit with Gaussian model
d_fit <- nest(SoG_rv, data = c(sev_thresh, rate)) %>%
  mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') - 10,
                                              start_upper = get_start_vals(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008') + 10,
                                              lower = get_lower_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$sev_thresh, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         # create new temperature data
         new_data = map(data, ~tibble(sev_thresh = seq(min(.x$sev_thresh), max(.x$sev_thresh), length.out = 100))),
         # predict over that data,
         preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))

# unnest predictions
d_preds_SoG <- select(d_fit, preds) %>%
  unnest(preds) %>% 
  mutate(SR = "Strait of Georgia")

# plot data and predictions
# ggplot() +
#   geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
#   geom_point(aes(temp, rate), SoG_rv, size = 2, alpha = 0.5) +
#   theme_bw(base_size = 12) +
#   labs(x = 'Temperature (ºC)',
#        y = 'Feeding rate',
#        title = 'Strait of Georgia')

#SoG_rv: Quadratic: refit model using nlsLM
fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = sev_thresh, a, b, c),
                               data = SoG_rv,
                               start = coef(d_fit$quadratic[[1]]),
                               lower = get_lower_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = 'quadratic_2008'),
                               upper = get_upper_lims(SoG_rv$sev_thresh, SoG_rv$rate, model_name = 'quadratic_2008'),
                               weights = rep(1, times = nrow(SoG_rv)))

# bootstrap using case resampling
boot1 <- Boot(fit_nlsLM, method = 'case')

# look at the data
# head(boot1$t)

# hist(boot1, layout = c(2,2))

#SoG_rv: quadratic: Now plot the bootstrapped models
#create predictions of each bootstrapped model
boot1_preds <- boot1$t %>%
  as.data.frame() %>%
  drop_na() %>%
  mutate(iter = 1:n()) %>%
  group_by_all() %>%
  do(data.frame(sev_thresh = seq(min(SoG_rv$sev_thresh), max(SoG_rv$sev_thresh), length.out = 100))) %>%
  ungroup() %>%
  mutate(pred = quadratic_2008(temp = sev_thresh, a, b, c))

# calculate bootstrapped confidence intervals
boot1_conf_preds_SoG <- group_by(boot1_preds, sev_thresh) %>%
  summarise(conf_lower = quantile(pred, 0.025),
            conf_upper = quantile(pred, 0.975)) %>%
  ungroup() %>% 
  mutate(SR = "Strait of Georgia")

# plot bootstrapped CIs
p = ggplot() +
  geom_line(aes(sev_thresh, .fitted), d_preds_SoG, col = 'blue') +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
  geom_point(aes(sev_thresh, rate), SoG_rv, size = 2, alpha = 0.5) +
  theme_bw(base_size = 17) +
  labs(x = 'Temperature (ºC)',
       y = 'Feeding rate',
       title = paste0('TPC during ', unique(feeding_clean_total$Stage)[i]))

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
         RV = "fr",
         model = "quadratic")

# plot --------
p = ggplot() +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_CC, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun=mean, geom="point", size = 5) +
  stat_summary(aes(y = rate, x = sev_thresh, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 1.5) +
  geom_line(aes(sev_thresh, .fitted, col = SR), d_preds_SoG, linewidth = 2) +
  geom_ribbon(aes(sev_thresh, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
  scale_colour_manual(values = c("blue", "red")) +
  scale_fill_manual(values = c("darkblue", "darkred")) +
  labs(x = expression("Severity level (" * alpha[i] * ")"),
       y = expression("g of " * italic(M.~trossulus) *
                            " / g of " * italic(N.~lamellosa) * " / week"),
       col = "Source Region",
       fill = "Source Region",
       title = "Per capita weekly feeding rate") + 
  theme_cowplot(35) + 
  scale_x_continuous(breaks = c(-3,-2,-1,0,1,2,3,4,5,6,7,8)) +
  expand_limits(x = c(-2.3, 7.3)) +
  scale_y_continuous(breaks = c(0,1,2,3)) +
  expand_limits(y = c(-0.5,3)) +
  geom_hline(aes(yintercept = 0), linetype = 2,
             linewidth = 2) +
  annotate("segment",
           x = ci_params_select_SoG_fr$estimate[1],
           xend = ci_params_select_SoG_fr$estimate[1],
           y = -Inf,
           yend = 1.4,
           linetype = 2,
           linewidth = 2) +
  annotate("segment",
           x = ci_params_select_CC_fr$estimate[1],
           xend = ci_params_select_CC_fr$estimate[1],
           y = -Inf,
           yend = 1.4,
           linetype = 2,
           linewidth = 2) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
p
# Store the plot in the list
plots[[1]] <- p

alpha_opt_max_results = rbind(alpha_opt_max_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                                sev_thresh = ci_params_select_SoG_fr$estimate[1], 
                                                                conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                                conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                                SR = ci_params_select_SoG_fr$SR[1],
                                                                model = ci_params_select_SoG_fr$model[1],
                                                                RV = "FR"))

alpha_opt_max_results = rbind(alpha_opt_max_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                                sev_thresh = ci_params_select_SoG_fr$estimate[2], 
                                                                conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                                conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                                SR = ci_params_select_SoG_fr$SR[2],
                                                                model = ci_params_select_SoG_fr$model[2],
                                                                RV = "FR"))

# statistical significance tests 
cc_se = (alpha_opt_max_results$conf_upper[2]-alpha_opt_max_results$conf_lower[2])/(2*1.96)
sog_se = (alpha_opt_max_results$conf_upper[4]-alpha_opt_max_results$conf_lower[4])/(2*1.96)
diff <- alpha_opt_max_results$sev_thresh[2] - alpha_opt_max_results$sev_thresh[4]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value

# statistical significance tests - alpha opt
cc_se = (alpha_opt_max_results$conf_upper[1]-alpha_opt_max_results$conf_lower[1])/(2*1.96)
sog_se = (alpha_opt_max_results$conf_upper[3]-alpha_opt_max_results$conf_lower[3])/(2*1.96)
diff <- alpha_opt_max_results$sev_thresh[1] - alpha_opt_max_results$sev_thresh[3]
se_diff <- sqrt(cc_se^2 + sog_se^2)

z <- diff / se_diff
p_value <- 2 * (1 - pnorm(abs(z)))

z
p_value # 0 

# compute mean and CI for sev threshs - CC
vals = c(topt_sev_results$op_sev[2],topt_sev_results$op_sev[6],topt_sev_results$op_sev[10],alpha_opt_max_results$sev_thresh[2])
mean_val <- mean(vals)
n <- length(vals)
se <- sd(vals) / sqrt(n)
t_crit <- qt(0.975, df = n - 1)

lower <- mean_val - t_crit * se
upper <- mean_val + t_crit * se
t_crit * se

# compute mean and CI for sev threshs - SoG
vals = c(topt_sev_results$op_sev[4],topt_sev_results$op_sev[8],topt_sev_results$op_sev[12],alpha_opt_max_results$sev_thresh[4])
mean_val <- mean(vals)
n <- length(vals)
se <- sd(vals) / sqrt(n)
t_crit <- qt(0.975, df = n - 1)

lower <- mean_val - t_crit * se
upper <- mean_val + t_crit * se
t_crit * se

# compute mean and CI for sev threshs - both
vals = c(topt_sev_results$op_sev[2],topt_sev_results$op_sev[6],topt_sev_results$op_sev[10],alpha_opt_max_results$sev_thresh[2],
         topt_sev_results$op_sev[4],topt_sev_results$op_sev[8],topt_sev_results$op_sev[12],alpha_opt_max_results$sev_thresh[4])
mean_val <- mean(vals)
n <- length(vals)
se <- sd(vals) / sqrt(n)
t_crit <- qt(0.975, df = n - 1)

lower <- mean_val - t_crit * se
upper <- mean_val + t_crit * se
t_crit * se

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/Figures")
ggsave(filename="main_fig_v4.png", height=20, width=27, 
       plot=grid.arrange(plots_SL[[1]],plots_ShW[[5]],
                         plots[[1]],plots_TiW[[5]],
                         ncol = 2))
