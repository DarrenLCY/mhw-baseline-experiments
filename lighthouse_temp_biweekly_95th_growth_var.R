#Visualizing the long-term temperature data at Calvert & Nanaimo
#Data downloaded from https://open.canada.ca/data/en/dataset/719955f2-bf8e-44f7-bc26-6bd623e82884
#Last edited Sept 2022

#Load in necessary packages----
pkgs <- c("janitor", "dplyr", "tidyverse", "ggplot2", "zoo", "lubridate", "cowplot")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")

#Load csvs & summarize temp data----
#from Egg Island (closest lighthouse to Calvert) & Departure Bay (closest to Nanaimo sites)
egg <- read.csv("Egg_Island_-_Daily_Sea_Surface_Temperature_and_Salinity_1970-2021.csv")
departure <- read.csv("Departure_Bay_PBS_-_Daily_Sea_Surface_Temperature_and_Salinity_1914-2021.csv")

#Clean dataset and filter for 2012-2020 data
egg_1 <- egg[ , c(1,3)] %>% 
  row_to_names(row_number = 1) %>% 
  rename("date" = "DATE (YYYY-MM-DD)", "temp" = "TEMPERATURE ( C )") %>% 
  mutate(date = ymd(date),
         temp = as.numeric(temp),
         station = "Egg Island, Central Coast",
         md_format = format(date, "%m-%d")) %>% 
  filter(date > "1970-01-01" & date < "2020-01-01", 
         temp != 999.9) %>% 
  filter(md_format > "06-13" & md_format < "07-26") %>% 
  mutate(week = ifelse((md_format >= "06-13" & md_format <= "06-29"), "week2",
                       ifelse((md_format >= "06-29" & md_format <= "07-12"), "week4",
                       ifelse((md_format >= "07-12" & md_format <= "07-26"), "week6", NA))))

departure_1 <- departure[ , c(1,3)] %>% 
  row_to_names(row_number = 1) %>% 
  rename("date" = "DATE (YYYY-MM-DD)", "temp" = "TEMPERATURE ( C )") %>% 
  mutate(date = ymd(date),
         temp = as.numeric(temp),
         station = "Departure Bay, Strait of Georgia",
         md_format = format(date, "%m-%d")) %>% 
  filter(date > "1970-01-01" & date < "2020-01-01", 
         temp != 999.9) %>% 
  filter(md_format > "06-13" & md_format < "07-26") %>% 
  mutate(week = ifelse((md_format >= "06-13" & md_format <= "06-29"), "week2",
                       ifelse((md_format >= "06-29" & md_format <= "07-12"), "week4",
                              ifelse((md_format >= "07-12" & md_format <= "07-26"), "week6", NA))))

# you want the 95th percentile of all the temps in the week prior to the FR measurement, from 2012 to 2020. 
# for instance, for the TPC on 2018-07-25, we'll filter for all temps from 2018-07-19 to 2018-07-25 inclusive

egg_2 <- egg_1 %>% 
  mutate(week = as.factor(week)) %>% 
  group_by(week) %>% 
  summarize(avgtemp = mean(temp), maxtemp = max(temp), sdtemp = sd(temp), temp90th = quantile(temp, 0.90),
            temp95th = quantile(temp, .95))

departure_2 <- departure_1 %>% 
  mutate(week = as.factor(week)) %>% 
  group_by(week) %>% 
  summarize(avgtemp = mean(temp), maxtemp = max(temp), sdtemp = sd(temp), temp90th = quantile(temp, 0.90),
            temp95th = quantile(temp, .95))

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")
saveRDS(egg_2, "CC_lighthouse_biweekly_95th_temp_growth_var.Rds")
saveRDS(departure_2, "SoG_lighthouse_biweekly_95th_temp_growth_var.Rds")
