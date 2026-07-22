# Script to analyze CTmax and TSM changes for feeding rate
# first part of the code is identical to feeding_rate v2.R script
# Actual analyses start on: # Central Coast CTmax changes ----------

#Load packages----
pkgs <- c("tidyverse", "lubridate", "car", "visreg", "cowplot", "survminer", "survival",
          "emmeans", "lme4")
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
feeding_base <- read.csv("per_weight_feeding_rate.csv")

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
         TiW = TW-ShW)

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

# -------------- now that we have the daily temps, we can compute weekly average temps
# later, you have to individually compute the weekly temps on a case by case basis 
T11_logger_data_weekly <- T11_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "11.0")

T12.5_logger_data_weekly <- T12.5_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "12.5")

T14_logger_data_weekly <- T14_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "14.0")

T15.5_logger_data_weekly <- T15.5_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "15.5")

T17_logger_data_weekly <- T17_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "17.0")

T18_logger_data_weekly <- T18_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "18.0")

T19_logger_data_weekly <- T19_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "19.0")

T20_logger_data_weekly <- T20_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "20.0")

T20.5_logger_data_weekly <- T20.5_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "20.5")

T21_logger_data_weekly <- T21_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "21.0")

T21.5_logger_data_weekly <- T21.5_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "21.5")

T22_logger_data_weekly <- T22_logger_data_daily %>%
  mutate(group = (row_number() - 1) %/% 7 + 1) %>%
  mutate(group = if_else(row_number() > n() - 5, 6, group)) %>% # include the extra days for the last period
  group_by(group) %>%
  summarise(
    start_date = min(date),
    end_date = max(date),
    avg_temp = mean(daily_avg_temp, na.rm = TRUE)
  ) %>%
  mutate(Treat = "22.0")

# ---------------- match the daily temps to the Hobday thresholds to determine the 'realised' severity threshold
# we will match them for 3 different periods (weeks 0-2, 2-4 and 4-6), due to the variable realised temperatures

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments")

# OISST data from UCAR
OISST_nc = nc_open("oisst_of_snails_beaty.nc")
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

#let's see how it looks like first... Run experimental_temps.R until you have the final plot
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

# ------ realised severity threshs for SoG

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

# 6-week heatwave, from 13 June to 26 July - T11
T11_logger_data_weekly$cc_sev_thresh = NA
T11_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev11 %>% filter(t >= T11_logger_data_weekly$start_date[i] & t <= T11_logger_data_weekly$end_date[i])
  subset_realised_temp = T11_logger_data_daily %>% filter(date >= T11_logger_data_weekly$start_date[i] & date <= T11_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T11_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev11 %>% filter(t >= T11_logger_data_weekly$start_date[i] & t <= T11_logger_data_weekly$end_date[i])
  subset_realised_temp = T11_logger_data_daily %>% filter(date >= T11_logger_data_weekly$start_date[i] & date <= T11_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T11_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T12.5
T12.5_logger_data_weekly$cc_sev_thresh = NA
T12.5_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev12.5 %>% filter(t >= T12.5_logger_data_weekly$start_date[i] & t <= T12.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T12.5_logger_data_daily %>% filter(date >= T12.5_logger_data_weekly$start_date[i] & date <= T12.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T12.5_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev12.5 %>% filter(t >= T12.5_logger_data_weekly$start_date[i] & t <= T12.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T12.5_logger_data_daily %>% filter(date >= T12.5_logger_data_weekly$start_date[i] & date <= T12.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T12.5_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T14
T14_logger_data_weekly$cc_sev_thresh = NA
T14_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev14 %>% filter(t >= T14_logger_data_weekly$start_date[i] & t <= T14_logger_data_weekly$end_date[i])
  subset_realised_temp = T14_logger_data_daily %>% filter(date >= T14_logger_data_weekly$start_date[i] & date <= T14_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T14_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev14 %>% filter(t >= T14_logger_data_weekly$start_date[i] & t <= T14_logger_data_weekly$end_date[i])
  subset_realised_temp = T14_logger_data_daily %>% filter(date >= T14_logger_data_weekly$start_date[i] & date <= T14_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T14_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T15.5
T15.5_logger_data_weekly$cc_sev_thresh = NA
T15.5_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev15.5 %>% filter(t >= T15.5_logger_data_weekly$start_date[i] & t <= T15.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T15.5_logger_data_daily %>% filter(date >= T15.5_logger_data_weekly$start_date[i] & date <= T15.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T15.5_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev15.5 %>% filter(t >= T15.5_logger_data_weekly$start_date[i] & t <= T15.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T15.5_logger_data_daily %>% filter(date >= T15.5_logger_data_weekly$start_date[i] & date <= T15.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T15.5_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T17
T17_logger_data_weekly$cc_sev_thresh = NA
T17_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev17 %>% filter(t >= T17_logger_data_weekly$start_date[i] & t <= T17_logger_data_weekly$end_date[i])
  subset_realised_temp = T17_logger_data_daily %>% filter(date >= T17_logger_data_weekly$start_date[i] & date <= T17_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T17_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev17 %>% filter(t >= T17_logger_data_weekly$start_date[i] & t <= T17_logger_data_weekly$end_date[i])
  subset_realised_temp = T17_logger_data_daily %>% filter(date >= T17_logger_data_weekly$start_date[i] & date <= T17_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T17_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T18
T18_logger_data_weekly$cc_sev_thresh = NA
T18_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev18 %>% filter(t >= T18_logger_data_weekly$start_date[i] & t <= T18_logger_data_weekly$end_date[i])
  subset_realised_temp = T18_logger_data_daily %>% filter(date >= T18_logger_data_weekly$start_date[i] & date <= T18_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T18_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev18 %>% filter(t >= T18_logger_data_weekly$start_date[i] & t <= T18_logger_data_weekly$end_date[i])
  subset_realised_temp = T18_logger_data_daily %>% filter(date >= T18_logger_data_weekly$start_date[i] & date <= T18_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T18_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T19
T19_logger_data_weekly$cc_sev_thresh = NA
T19_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev19 %>% filter(t >= T19_logger_data_weekly$start_date[i] & t <= T19_logger_data_weekly$end_date[i])
  subset_realised_temp = T19_logger_data_daily %>% filter(date >= T19_logger_data_weekly$start_date[i] & date <= T19_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T19_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev19 %>% filter(t >= T19_logger_data_weekly$start_date[i] & t <= T19_logger_data_weekly$end_date[i])
  subset_realised_temp = T19_logger_data_daily %>% filter(date >= T19_logger_data_weekly$start_date[i] & date <= T19_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T19_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T20
T20_logger_data_weekly$cc_sev_thresh = NA
T20_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev20 %>% filter(t >= T20_logger_data_weekly$start_date[i] & t <= T20_logger_data_weekly$end_date[i])
  subset_realised_temp = T20_logger_data_daily %>% filter(date >= T20_logger_data_weekly$start_date[i] & date <= T20_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T20_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev20 %>% filter(t >= T20_logger_data_weekly$start_date[i] & t <= T20_logger_data_weekly$end_date[i])
  subset_realised_temp = T20_logger_data_daily %>% filter(date >= T20_logger_data_weekly$start_date[i] & date <= T20_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T20_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T20.5
T20.5_logger_data_weekly$cc_sev_thresh = NA
T20.5_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev20.5 %>% filter(t >= T20.5_logger_data_weekly$start_date[i] & t <= T20.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T20.5_logger_data_daily %>% filter(date >= T20.5_logger_data_weekly$start_date[i] & date <= T20.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T20.5_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev20.5 %>% filter(t >= T20.5_logger_data_weekly$start_date[i] & t <= T20.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T20.5_logger_data_daily %>% filter(date >= T20.5_logger_data_weekly$start_date[i] & date <= T20.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T20.5_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T21
T21_logger_data_weekly$cc_sev_thresh = NA
T21_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev21 %>% filter(t >= T21_logger_data_weekly$start_date[i] & t <= T21_logger_data_weekly$end_date[i])
  subset_realised_temp = T21_logger_data_daily %>% filter(date >= T21_logger_data_weekly$start_date[i] & date <= T21_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T21_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev21 %>% filter(t >= T21_logger_data_weekly$start_date[i] & t <= T21_logger_data_weekly$end_date[i])
  subset_realised_temp = T21_logger_data_daily %>% filter(date >= T21_logger_data_weekly$start_date[i] & date <= T21_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T21_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T21.5
T21.5_logger_data_weekly$cc_sev_thresh = NA
T21.5_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev21.5 %>% filter(t >= T21.5_logger_data_weekly$start_date[i] & t <= T21.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T21.5_logger_data_daily %>% filter(date >= T21.5_logger_data_weekly$start_date[i] & date <= T21.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T21.5_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev21.5 %>% filter(t >= T21.5_logger_data_weekly$start_date[i] & t <= T21.5_logger_data_weekly$end_date[i])
  subset_realised_temp = T21.5_logger_data_daily %>% filter(date >= T21.5_logger_data_weekly$start_date[i] & date <= T21.5_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T21.5_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# 6-week heatwave, from 13 June to 26 July - T22
T22_logger_data_weekly$cc_sev_thresh = NA
T22_logger_data_weekly$sog_sev_thresh = NA

all_sevs = seq(-3,10,by=0.1)
for (i in 1:6) {
  # SoG 
  subset_target_temp = sog_sev22 %>% filter(t >= T22_logger_data_weekly$start_date[i] & t <= T22_logger_data_weekly$end_date[i])
  subset_realised_temp = T22_logger_data_daily %>% filter(date >= T22_logger_data_weekly$start_date[i] & date <= T22_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T22_logger_data_weekly$sog_sev_thresh[i] = chosen_sev
  
  # CC
  subset_target_temp = cc_sev22 %>% filter(t >= T22_logger_data_weekly$start_date[i] & t <= T22_logger_data_weekly$end_date[i])
  subset_realised_temp = T22_logger_data_daily %>% filter(date >= T22_logger_data_weekly$start_date[i] & date <= T22_logger_data_weekly$end_date[i])
  current_diff = 100
  for (sev in all_sevs) {
    sev_of_interest = subset_target_temp$seas + subset_target_temp$thresh_less_seas*sev
    if (abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp)) < current_diff) {
      chosen_sev = sev
      current_diff = abs(sum(sev_of_interest-subset_realised_temp$daily_avg_temp))
    }
  }
  T22_logger_data_weekly$cc_sev_thresh[i] = chosen_sev
}

# combine all weekly temps into one df
all_weekly_temps = rbind(T11_logger_data_weekly, T12.5_logger_data_weekly, T14_logger_data_weekly, T15.5_logger_data_weekly,
                         T17_logger_data_weekly, T18_logger_data_weekly, T19_logger_data_weekly, T20_logger_data_weekly,
                         T20.5_logger_data_weekly, T21_logger_data_weekly, T21.5_logger_data_weekly, T22_logger_data_weekly)

# combine all daily temps into one df
all_daily_temps = rbind(T11_logger_data_daily, T12.5_logger_data_daily, T14_logger_data_daily, T15.5_logger_data_daily,
                        T17_logger_data_daily, T18_logger_data_daily, T19_logger_data_daily, T20_logger_data_daily,
                        T20.5_logger_data_daily, T21_logger_data_daily, T21.5_logger_data_daily, T22_logger_data_daily)

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
  if (feeding_clean_total$Stage[i] == "wk2" | feeding_clean_total$Stage[i] == "wk4" | feeding_clean_total$Stage[i] == "wk6") {
    feeding_clean_total$sum_TW_snails[i] = sum(growth_base$Total_weight[which(substr(growth_base$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                                growth_base$Stage == feeding_clean_total$Stage[i] &
                                                                                growth_base$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
    feeding_clean_total$sum_TiW_snails[i] = sum(growth_clean$TiW[which(substr(growth_clean$SP, 1, 1) == feeding_clean_total$SP[i] &
                                                                         growth_clean$Stage == feeding_clean_total$Stage[i] &
                                                                         growth_clean$Treat == feeding_clean_total$Treatment[i])], na.rm = TRUE)
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

# absolute temps ------------------------------------- 

# Create an empty list to store plots
plots <- list()

topt_ctmax_results = data.frame(param = vector(), temp = vector(), conf_lower = vector(), 
                                conf_upper = vector(), SR = vector(), model = vector(), period = vector(),
                                RV = vector())

for (i in 1:length(unique(feeding_clean_total$Stage))) {
  
  # Central Coast
  CC_rv <- feeding_clean_total %>% 
    filter((SP == "K" & Stage == unique(feeding_clean_total$Stage)[i]) |
             (SP == "P" & Stage == unique(feeding_clean_total$Stage)[i])) %>% 
    select(mean_T, meanW_mussels_per_snail_TiW) %>% 
    rename("rate" = "meanW_mussels_per_snail_TiW")
  CC_rv$mean_T = as.numeric(as.character(CC_rv$mean_T))
  CC_rv$SR = "Central Coast"
  
  # fit five chosen model formulations in rTPC
  d_fits <- nest(CC_rv, data = c(mean_T, rate)) %>%
    mutate(briere = map(data, ~nls_multstart(rate~briere2_1999(mean_T, tmin, tmax, a, b),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'briere2_1999') - 10,
                                             start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'briere2_1999') + 10,
                                             lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'briere2_1999'),
                                             upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'briere2_1999'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
           gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = mean_T, rmax, topt, a),
                                               data = .x,
                                               iter = c(4,4,4),
                                               start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'gaussian_1987') - 10,
                                               start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'gaussian_1987') + 10,
                                               lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'gaussian_1987'),
                                               upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'gaussian_1987'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
           quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = mean_T, a, b, c),
                                                data = .x,
                                                iter = c(4,4,4),
                                                start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                                start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                                lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', briere:quadratic)
  
  # get predictions using augment
  newdata <- tibble(mean_T = seq(min(CC_rv$mean_T), max(CC_rv$mean_T), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, mean_T < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(mean_T, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(mean_T, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(mean_T, rate), CC_rv) +
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
  cc_best_fr <- ggplot(d_preds, aes(mean_T, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(mean_T, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(mean_T, rate), CC_rv) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Feeding rate',
         title = 'Central Coast') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(CC_rv, aes(mean_T, rate)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Rate')
  
  #CC_rv: Quadratic: Fit data
  # fit with Gaussian model
  d_fit <- nest(CC_rv, data = c(mean_T, rate)) %>%
    mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = mean_T, a, b, c),
                                                data = .x,
                                                iter = c(4,4,4),
                                                start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') - 10,
                                                start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') + 10,
                                                lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
           # create new temperature data
           new_data = map(data, ~tibble(mean_T = seq(min(.x$mean_T), max(.x$mean_T), length.out = 100))),
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
  fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = mean_T, a, b, c),
                                 data = CC_rv,
                                 start = coef(d_fit$quadratic[[1]]),
                                 lower = get_lower_lims(CC_rv$mean_T, CC_rv$rate, model_name = 'quadratic_2008'),
                                 upper = get_upper_lims(CC_rv$mean_T, CC_rv$rate, model_name = 'quadratic_2008'),
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
    do(data.frame(mean_T = seq(min(CC_rv$mean_T), max(CC_rv$mean_T), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = quadratic_2008(temp = mean_T, a, b, c))
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_CC <- group_by(boot1_preds, mean_T) %>%
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
  
  topt_ctmax_results = rbind(topt_ctmax_results, data.frame(param = ci_params_select_CC_fr$param[1],
                                                            temp = ci_params_select_CC_fr$estimate[1], 
                                                            conf_lower = ci_params_select_CC_fr$conf_lower[1], 
                                                            conf_upper = ci_params_select_CC_fr$conf_upper[1],
                                                            SR = ci_params_select_CC_fr$SR[1],
                                                            model = ci_params_select_CC_fr$model[1],
                                                            period = i+1, RV = "FR"))
  
  topt_ctmax_results = rbind(topt_ctmax_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                            temp = ci_params_select_CC_fr$estimate[2], 
                                                            conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                            conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                            SR = ci_params_select_CC_fr$SR[2],
                                                            model = ci_params_select_CC_fr$model[2],
                                                            period = i+1, RV = "FR"))
  
  # Strait of Georgia
  SoG_rv <- feeding_clean_total %>% 
    filter((SP == "C" & Stage == unique(feeding_clean_total$Stage)[i]) |
             (SP == "H" & Stage == unique(feeding_clean_total$Stage)[i])) %>% 
    select(mean_T, meanW_mussels_per_snail_TiW) %>% 
    rename("rate" = "meanW_mussels_per_snail_TiW")
  SoG_rv$mean_T = as.numeric(as.character(SoG_rv$mean_T))
  SoG_rv$SR = "Strait of Georgia"
  
  # fit five chosen model formulations in rTPC
  d_fits <- nest(SoG_rv, data = c(mean_T, rate)) %>%
    mutate(briere = map(data, ~nls_multstart(rate~briere2_1999(mean_T, tmin, tmax, a, b),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'briere2_1999') - 10,
                                             start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'briere2_1999') + 10,
                                             lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'briere2_1999'),
                                             upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'briere2_1999'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
           gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = mean_T, rmax, topt, a),
                                               data = .x,
                                               iter = c(4,4,4),
                                               start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'gaussian_1987') - 10,
                                               start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'gaussian_1987') + 10,
                                               lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'gaussian_1987'),
                                               upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'gaussian_1987'),
                                               supp_errors = 'Y',
                                               convergence_count = FALSE)),
           quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = mean_T, a, b, c),
                                                data = .x,
                                                iter = c(4,4,4),
                                                start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                                start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                                lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', briere:quadratic)
  
  # get predictions using augment
  newdata <- tibble(mean_T = seq(min(SoG_rv$mean_T), max(SoG_rv$mean_T), length.out = 100))
  d_preds <- d_stack %>%
    mutate(., preds = map(fit, augment, newdata = newdata)) %>%
    select(-fit) %>%
    unnest(preds)
  
  # take a random point from each model for labelling
  d_labs <- filter(d_preds, mean_T < 30) %>%
    group_by(., model_name) %>%
    sample_n(., 1) %>%
    ungroup()
  
  # plot
  # ggplot(d_preds, aes(mean_T, .fitted)) +
  #   geom_line(aes(col = model_name)) +
  #   geom_label_repel(aes(mean_T, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  #   geom_point(aes(mean_T, rate), SoG_rv) +
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
  SoG_best_fr <- ggplot(d_preds, aes(mean_T, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(mean_T, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(mean_T, rate), SoG_rv) +
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
  d_fit <- nest(SoG_rv, data = c(mean_T, rate)) %>%
    mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = mean_T, a, b, c),
                                                data = .x,
                                                iter = c(4,4,4),
                                                start_lower = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') - 10,
                                                start_upper = get_start_vals(.x$mean_T, .x$rate, model_name = 'quadratic_2008') + 10,
                                                lower = get_lower_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$mean_T, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
           # create new temperature data
           new_data = map(data, ~tibble(mean_T = seq(min(.x$mean_T), max(.x$mean_T), length.out = 100))),
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
  fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = mean_T, a, b, c),
                                 data = SoG_rv,
                                 start = coef(d_fit$quadratic[[1]]),
                                 lower = get_lower_lims(SoG_rv$mean_T, SoG_rv$rate, model_name = 'quadratic_2008'),
                                 upper = get_upper_lims(SoG_rv$mean_T, SoG_rv$rate, model_name = 'quadratic_2008'),
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
    do(data.frame(mean_T = seq(min(SoG_rv$mean_T), max(SoG_rv$mean_T), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = quadratic_2008(temp = mean_T, a, b, c))
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_SoG <- group_by(boot1_preds, mean_T) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot bootstrapped CIs
  p = ggplot() +
    geom_line(aes(mean_T, .fitted), d_preds_SoG, col = 'blue') +
    geom_ribbon(aes(mean_T, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_SoG, fill = 'blue', alpha = 0.3) +
    geom_point(aes(mean_T, rate), SoG_rv, size = 2, alpha = 0.5) +
    theme_bw(base_size = 17) +
    labs(x = 'Temperature (ºC)',
         y = 'Feeding rate',
         title = paste0('TPC during ', unique(feeding_clean_total$Stage)[i]))
  
  # plot --------
  p = ggplot() +
    stat_summary(aes(y = rate, x = mean_T, col = SR), data = CC_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = mean_T, col = SR), data = CC_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(mean_T, .fitted, col = SR), d_preds_CC, linewidth = 1) +
    geom_ribbon(aes(mean_T, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
    stat_summary(aes(y = rate, x = mean_T, col = SR), data = SoG_rv, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = mean_T, col = SR), data = SoG_rv, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(mean_T, .fitted, col = SR), d_preds_SoG, linewidth = 1) +
    geom_ribbon(aes(mean_T, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
    scale_colour_manual(values = c("skyblue", "coral")) +
    scale_fill_manual(values = c("skyblue3", "coral3")) +
    labs(x = 'Temperature (ºC)',
         y = if (i == 1 | i == 4) "Per capita weekly feeding rate" else NULL,
         col = "Source Region",
         fill = "Source Region",
         title = paste0('TPC during week ', i+1)) + 
    theme_cowplot(16) + 
    scale_x_continuous(breaks = c(10,12,14,16,18,20,22)) +
    expand_limits(x = c(10,24)) +
    scale_y_continuous(breaks = c(0,1,2,3)) +
    expand_limits(y = c(-0.5,3)) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text  = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  # Store the plot in the list
  plots[[i]] <- p
  
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
  
  topt_ctmax_results = rbind(topt_ctmax_results, data.frame(param = ci_params_select_SoG_fr$param[1],
                                                            temp = ci_params_select_SoG_fr$estimate[1], 
                                                            conf_lower = ci_params_select_SoG_fr$conf_lower[1], 
                                                            conf_upper = ci_params_select_SoG_fr$conf_upper[1],
                                                            SR = ci_params_select_SoG_fr$SR[1],
                                                            model = ci_params_select_SoG_fr$model[1],
                                                            period = i+1, RV = "FR"))
  
  topt_ctmax_results = rbind(topt_ctmax_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                            temp = ci_params_select_SoG_fr$estimate[2], 
                                                            conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                            conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                            SR = ci_params_select_SoG_fr$SR[2],
                                                            model = ci_params_select_SoG_fr$model[2],
                                                            period = i+1, RV = "FR"))
  
}

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/Figures")
ggsave(filename="TPCs_temp_FR.png", height=10, width=15, 
       plot=grid.arrange(plots[[1]],plots[[2]],plots[[3]],
                         plots[[4]],plots[[5]],ncol = 3))

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/")
# write.csv(topt_ctmax_results, "Topt_and_CTmax_FR.csv")
topt_ctmax_results = read.csv("Topt_and_CTmax_FR.csv")

# Central Coast CTmax changes ----------

get_mode <- function(x) {
  uniq_x <- unique(x)
  uniq_x[which.max(tabulate(match(x, uniq_x)))]
}

ctmax_results_CC = topt_ctmax_results[which(topt_ctmax_results$SR == "Central Coast" &
                                              topt_ctmax_results$param == "ctmax"),]
ctmax_results_CC = ctmax_results_CC[2:5,] # filter out first period cos the TPCs don't look good

# get modes for the different periods
get_mode(feeding_clean_total$Date[which(feeding_clean_total$Stage == "wk2")]) # 27-Jun-25
get_mode(feeding_clean_total$Date[which(feeding_clean_total$Stage == "wk3")]) # 4-Jul-25
get_mode(feeding_clean_total$Date[which(feeding_clean_total$Stage == "wk4")]) # 12-Jul-25
get_mode(feeding_clean_total$Date[which(feeding_clean_total$Stage == "wk5")]) # 18-Jul-25
get_mode(feeding_clean_total$Date[which(feeding_clean_total$Stage == "wk6")]) # 26-Jul-25

ctmax_results_CC$date = as.Date("4-Jul-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_CC$date[1] = as.Date("4-Jul-25", format = "%d-%b-%y")
ctmax_results_CC$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_CC$date[3] = as.Date("18-Jul-25", format = "%d-%b-%y")
ctmax_results_CC$date[4] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_CC, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/")
# write.csv(ctmax_results_CC, "FR_CTmax_changes_CC.csv")
ctmax_results_CC = read.csv("FR_CTmax_changes_CC.csv")
ctmax_results_CC$date = as.Date(ctmax_results_CC$date)

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
start_date <- as.Date("2018-07-04")
end_date <- as.Date("2018-07-26")
ts2 = ts %>%
  filter(t >= start_date & t <= end_date)

# change the dates from 2018 to 2025 in ts2 (doesn't make any difference cos it's just the clim which is the same every year anyway)
ts2$t <- seq(
  from = as.Date("2025-07-04"),
  to   = as.Date("2025-07-26"),
  by   = "day"
)

# change in climatology
max(ts2$seas)-min(ts2$seas)
max(ts2$doy)-min(ts2$doy)

ctmax_plots_CC = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+11), color="red", size = 3) +
  geom_smooth(aes(y=ctmax_results_CC$temp, x=ctmax_results_CC$date), method = "lm", se = TRUE, color = "blue", linewidth = 3) +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 35)) +
  labs(x = "Date") + ggtitle("Central Coast")

# ggsave(filename="FR_CTmax_changes_CC.png", height=8, width=15, 
#        plot=ctmax_plots_CC)

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_CC$date,
  temp  = ctmax_results_CC$temp
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

# regression plot CTmax CC --------------
ctmax_results_CC <- ctmax_results_CC %>% 
  mutate(ctmax_diff = c(NA,diff(ctmax_results_CC$temp)))
ts3 <- ts2 %>%
  filter(t %in% ctmax_results_CC$date)
ts3 <- ts3 %>% 
  mutate(seas_diff = c(NA,diff(ts3$seas)))
CC_reg_plot_ctmax = ggplot() +
  geom_point(aes(y=ctmax_results_CC$ctmax_diff, x=ts3$seas_diff), size=4) +
  geom_smooth(aes(y=ctmax_results_CC$ctmax_diff, x=ts3$seas_diff), method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", size=1) +
  labs(x = "Change in seasonal temperatures (°C)", y = "Change in CTmax (°C)")+
  theme_bw() + 
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) + ggtitle("Central Coast")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  seas_diff = ts3$seas_diff,
  ctmax_diff  = ctmax_results_CC$ctmax_diff
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(ctmax_diff ~ seas_diff, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["seas_diff"]]
pval  <- sm$coefficients["seas_diff", "Pr(>|t|)"]
ci    <- confint(m, "seas_diff", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))

library(metafor)
# Calculate standard errors from confidence intervals
ctmax_results_CC$se <- (ctmax_results_CC$conf_upper - ctmax_results_CC$conf_lower) / (2 * 1.96)

# Perform meta-analysis
meta_result <- rma(yi = temp, sei = se, data = ctmax_results_CC)

# Summary of the meta-analysis
summary(meta_result)

# Central Coast TSM changes ----------

# load lighthouse data
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")
CC_lighthouse_weekly_95th_temp = readRDS("CC_lighthouse_weekly_95th_temp_FR.Rds")

CC_TSM_df <- ctmax_results_CC %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(CC_lighthouse_weekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

CC_tsm_plots=ggplot(CC_TSM_df, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

CC_tsm_plots = CC_tsm_plots + 
  geom_line(data = ts2, aes(x = t, y = seas-4), color="red", size = 3) +
  geom_smooth(aes(y=CC_TSM_df$TSM, x=CC_TSM_df$date), method = "lm", se = TRUE, color = "blue", linewidth = 3) +
  scale_y_continuous(name = "TSM (°C)", sec.axis = sec_axis( trans=~.+4, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        text = element_text(size = 35)) +
  labs(x = "Date") 

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = CC_TSM_df$date,
  TSM  = CC_TSM_df$TSM
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

# regression plot TSM CC --------------
CC_TSM_df <- CC_TSM_df %>% 
  mutate(TSM_diff = c(NA,diff(CC_TSM_df$TSM)))
ts3 <- ts2 %>%
  filter(t %in% CC_TSM_df$date)
ts3 <- ts3 %>% 
  mutate(seas_diff = c(NA,diff(ts3$seas)))
CC_reg_plot_tsm = ggplot() +
  geom_point(aes(y=CC_TSM_df$TSM_diff, x=ts3$seas_diff), size=4) +
  geom_smooth(aes(y=CC_TSM_df$TSM_diff, x=ts3$seas_diff), method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", size=1) +
  labs(x = "Change in seasonal temperatures (°C)", y = "Change in TSM (°C)")+
  theme_bw() + 
  theme(text = element_text(size = 20),
        axis.title.y = element_blank())

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  seas_diff = ts3$seas_diff,
  TSM_diff  = CC_TSM_df$TSM_diff
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM_diff ~ seas_diff, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["seas_diff"]]
pval  <- sm$coefficients["seas_diff", "Pr(>|t|)"]
ci    <- confint(m, "seas_diff", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))

# Strait of Georgia CTmax changes ----------

ctmax_results_SoG = topt_ctmax_results[which(topt_ctmax_results$SR == "Strait of Georgia" &
                                              topt_ctmax_results$param == "ctmax"),]
ctmax_results_SoG = ctmax_results_SoG[2:5,] # filter out first period cos the TPCs don't look good

ctmax_results_SoG$date = as.Date("4-Jul-25", format = "%d-%b-%y") # you need this to make the whole column Date first (else it doesn't work)
ctmax_results_SoG$date[1] = as.Date("4-Jul-25", format = "%d-%b-%y")
ctmax_results_SoG$date[2] = as.Date("12-Jul-25", format = "%d-%b-%y")
ctmax_results_SoG$date[3] = as.Date("18-Jul-25", format = "%d-%b-%y")
ctmax_results_SoG$date[4] = as.Date("26-Jul-25", format = "%d-%b-%y")

ctmax_plots = ggplot(ctmax_results_SoG, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/")
# write.csv(ctmax_results_SoG, "FR_CTmax_changes_SoG.csv")
ctmax_results_SoG = read.csv("FR_CTmax_changes_SoG.csv")
ctmax_results_SoG$date = as.Date(ctmax_results_SoG$date)

# add seasonally varying climatology to the plot
OISST_of_interest = data.frame(t = oisst_dates,
                               temp = sst[,2])

ts = ts2clm(OISST_of_interest, climatologyPeriod = c("1983-01-01", "2012-12-31"))
start_date <- as.Date("2018-07-04")
end_date <- as.Date("2018-07-26")
ts2 = ts %>%
  filter(t >= start_date & t <= end_date)

# change the dates from 2018 to 2025 in ts2 (doesn't make any difference cos it's just the clim which is the same every year anyway)
ts2$t <- seq(
  from = as.Date("2025-07-04"),
  to   = as.Date("2025-07-26"),
  by   = "day"
)

# change in climatology
max(ts2$seas)-min(ts2$seas)
max(ts2$doy)-min(ts2$doy)

ctmax_plots_SoG = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+11), color="red", size = 3) +
  geom_smooth(aes(y=ctmax_results_SoG$temp, x=ctmax_results_SoG$date), method = "lm", se = TRUE, color = "blue", linewidth = 3) +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size = 35)) +
  labs(x = "Date") + ggtitle("Strait of Georgia")

# ggsave(filename="FR_CTmax_changes_SoG.png", height=8, width=15, 
#        plot=ctmax_plots_SoG)

# save both in one plot
# ggsave(filename="FR_CTmax_changes.png", height=10, width=15, 
#        plot=grid.arrange(ctmax_plots_CC,
#                          ctmax_plots_SoG,nrow=2))

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = ctmax_results_SoG$date,
  temp  = ctmax_results_SoG$temp
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

# regression plot CTmax SoG --------------
ctmax_results_SoG <- ctmax_results_SoG %>% 
  mutate(ctmax_diff = c(NA,diff(ctmax_results_SoG$temp)))
ts4 <- ts2 %>%
  filter(t %in% ctmax_results_SoG$date)
ts4 <- ts4 %>% 
  mutate(seas_diff = c(NA,diff(ts4$seas)))
SoG_reg_plot_ctmax = ggplot() +
  geom_point(aes(y=ctmax_results_SoG$ctmax_diff, x=ts4$seas_diff), size=4) +
  geom_smooth(aes(y=ctmax_results_SoG$ctmax_diff, x=ts4$seas_diff), method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", size=1) +
  labs(x = "Change in seasonal temperatures (°C)", y = expression("Change in CT" [max] * " (°C)"))+
  theme_bw() + 
  theme(text = element_text(size = 20),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "coral3",
                                  size = 35)) + ggtitle("Strait of Georgia")

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  seas_diff = ts4$seas_diff,
  ctmax_diff  = ctmax_results_SoG$ctmax_diff
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(ctmax_diff ~ seas_diff, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["seas_diff"]]
pval  <- sm$coefficients["seas_diff", "Pr(>|t|)"]
ci    <- confint(m, "seas_diff", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))

library(metafor)
# Calculate standard errors from confidence intervals
ctmax_results_SoG$se <- (ctmax_results_SoG$conf_upper - ctmax_results_SoG$conf_lower) / (2 * 1.96)

# Perform meta-analysis
meta_result <- rma(yi = temp, sei = se, data = ctmax_results_SoG)

# Summary of the meta-analysis
summary(meta_result)

# Strait of Georgia TSM changes ----------

# load lighthouse data
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")
SoG_lighthouse_weekly_95th_temp = readRDS("SoG_lighthouse_weekly_95th_temp_FR.Rds")

SoG_TSM_df <- ctmax_results_SoG %>% 
  select(c(date, temp, conf_lower, conf_upper)) %>% 
  cbind(SoG_lighthouse_weekly_95th_temp[,6]) %>% 
  mutate(TSM = temp - temp95th,
         date = as.Date(date),
         TSM_lower = conf_lower - temp95th,
         TSM_upper = conf_upper - temp95th)

SoG_tsm_plots=ggplot(SoG_TSM_df, aes(date, TSM)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = TSM_lower, ymax = TSM_upper)) +
  theme_bw(base_size = 20)

SoG_tsm_plots = SoG_tsm_plots + 
  geom_line(data = ts2, aes(x = t, y = seas-8), color="red", size = 3) +
  geom_smooth(aes(y=SoG_TSM_df$TSM, x=SoG_TSM_df$date), method = "lm", se = TRUE, color = "blue", linewidth = 3) +
  scale_y_continuous(name = "TSM (°C)", sec.axis = sec_axis( trans=~.+8, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        text = element_text(size = 35)) +
  labs(x = "Date") 

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  date = SoG_TSM_df$date,
  TSM  = SoG_TSM_df$TSM
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

# regression plot TSM SoG --------------
SoG_TSM_df <- SoG_TSM_df %>% 
  mutate(TSM_diff = c(NA,diff(SoG_TSM_df$TSM)))
ts4 <- ts2 %>%
  filter(t %in% SoG_TSM_df$date)
ts4 <- ts4 %>% 
  mutate(seas_diff = c(NA,diff(ts4$seas)))
SoG_reg_plot_tsm = ggplot() +
  geom_point(aes(y=SoG_TSM_df$TSM_diff, x=ts4$seas_diff), size=4) +
  geom_smooth(aes(y=SoG_TSM_df$TSM_diff, x=ts4$seas_diff), method = "lm", se = TRUE, color = "blue") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "red", size=1) +
  labs(x = "Change in seasonal temperatures (°C)", y = "Change in TSM (°C)")+
  theme_bw() + 
  theme(text = element_text(size = 20)) # Change all text size

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  seas_diff = ts4$seas_diff,
  TSM_diff  = SoG_TSM_df$TSM_diff
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(TSM_diff ~ seas_diff, data = dat)
sm <- summary(m)
# Extract slope (beta1), p-value, CI, R^2
beta1 <- coef(m)[["seas_diff"]]
pval  <- sm$coefficients["seas_diff", "Pr(>|t|)"]
ci    <- confint(m, "seas_diff", level = 0.95)
r2    <- sm$r.squared
cat(sprintf(
  "Slope = %.3f (95%% CI: %.3f to %.3f), p = %.3g, R^2 = %.3f\n",
  beta1, ci[1], ci[2], pval, r2
))
# Interpret direction and significance (two-sided test)
direction <- if (beta1 > 0) "positive" else if (beta1 < 0) "negative" else "zero"
sig_text  <- if (pval < 0.05) "statistically significant" else "not statistically significant"
cat(sprintf("Conclusion: The slope is %s and %s at α = 0.05.\n", direction, sig_text))

# save all plots
library(gridExtra)
ggsave(filename="FR_CTmax_TSM.png", height=13, width=15,
       plot=grid.arrange(ctmax_plots_SoG, 
                         ctmax_plots_CC,
                         SoG_tsm_plots,
                         CC_tsm_plots, nrow=2))

ggsave(filename="FR_CTmax_TSM_vs_seasonal_temp.png", height=13, width=15,
       plot=grid.arrange(SoG_reg_plot_ctmax, 
                         CC_reg_plot_ctmax,
                         SoG_reg_plot_tsm,
                         CC_reg_plot_tsm, nrow=2))
