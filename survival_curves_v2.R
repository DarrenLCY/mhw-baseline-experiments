# Script to plot Kaplan-Meier survival curves

#Load packages----
pkgs <- c("tidyverse", "lubridate", "car", "visreg", "cowplot", "survminer", "survival",
          "emmeans", "lme4", "waldo")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

#Load packages----
pkgs <- c("rTPC", "nls.multstart", "broom", "tidyverse", "cowplot", "ggrepel")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

rm(list = ls())
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/")

#Load csvs ----
growth_base <- read.csv("growth_variables.csv")
growth_base <- read.csv("growth_variables_v2.csv") # v2 is where i correct the IDs
# we agreed not to use the corrected one
# growth_base_corrected <- read.csv("growth_variables_corrected.csv")
# compare(growth_base, growth_base_corrected)
# compare(growth_base, growth_base_corrected, max_diffs = Inf)

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

colnames(growth_base)
unique(growth_base$Died) # "" - alive, "1"- dead, "N/A" - already dead in prior weeks

# convert Date.dead
growth_base <- growth_base %>%
  mutate(Date.dead = as.Date(Date.dead, format = "%d-%b-%y"))

# define experiment start + end
start_date <- as.Date("2025-06-13")
end_date   <- as.Date("2025-07-26")

# create survival dataset (ONE ROW PER SNAIL)
length(unique(growth_base$ID))
surv_df <- growth_base %>%
  select(ID, SP, Treat, Died, Date.dead) %>%
  distinct() %>%
  mutate(
    # time = days from start
    time = ifelse(Died == "1",
                  as.numeric(Date.dead - start_date),
                  as.numeric(end_date - start_date)),
    
    # status = 1 = died, 0 = survived
    status = ifelse(Died == "1", 1, 0),
    
    # population (rename cleanly)
    population = SP,
    
    # MHW treatment (numeric or factor)
    MHW = Treat
  )
unique(surv_df$Treat)
length(which(surv_df$Treat == "22"))

# try to remove NA
surv_df = surv_df[which(surv_df$Died != "N/A"),]

# check lengths
nrow(surv_df)
length(surv_df$time)
length(surv_df$status)

library(survival)
library(survminer)

surv_obj <- Surv(time = surv_df$time, event = surv_df$status)


km_fit <- survfit(surv_obj ~ population, data = surv_df)

ggsurvplot(
  km_fit,
  data = surv_df,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  ggtheme = theme_bw()
)

ggsurvplot(
  survfit(surv_obj ~ population, data = surv_df),
  data = surv_df,
  facet.by = "MHW",
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  ggtheme = theme_bw()
)


# treat MHW properly
surv_df$MHW <- as.factor(surv_df$MHW)

cox_model <- coxph(surv_obj ~ population * MHW, data = surv_df)

summary(cox_model)

surv_df <- surv_df %>%
  mutate(group = interaction(population, MHW))

fit <- survfit(surv_obj ~ group, data = surv_df)

ggsurvplot(
  fit,
  data = surv_df,
  pval = TRUE,
  ggtheme = theme_bw()
)

exp(coef(cox_model))
cox.zph(cox_model)

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub")
ggsave(filename="survival_curves_v2.png", height=8, width=15, plot=p)
