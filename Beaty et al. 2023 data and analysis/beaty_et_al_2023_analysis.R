# code adapted from Beaty et al. 2023

#Load packages----
pkgs <- c("tidyverse", "lubridate", "car", "visreg", "cowplot", "survminer", "survival",
          "emmeans", "lme4")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

rm(list = ls())
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/Beaty et al. 2023 data and analysis")

#Load csvs ----
#March 2nd: I saved the meso_collated from the Mesocosm folder into the CH 2 Analysis folder as a csv. If you make any changes to the base df, replace the one in git folder
#March 12th: CO24 was outlier, typo in raw datasheet, changed fmor 3.42 to 2.42 in csv
#March 4th: I saved the survival datasheet
#March 1th: I saved feeding rate datasheet from the 'Weekly feeding rate_Oct 23.xlsx'in the Meso > feeding rate data sheets folder
meso_base <- read.csv("meso_collated.csv")
meso_feed <- read.csv("meso_percap_feeding_rate.csv")

#Clean growth, feeding & survival dataframes----
#Add a Source Region column
#Ensure that all treatments line up with the correct tank
#Classify Tank 16 (which was supposed to be a 12) as a 15 degree treatment because you could never really bring that tank's temp down due to equipment issues
#Remove the Died & Notes column
#Rename column headers
#Remove the tanks that experienced chiller failures & 100% mortality within first experimental days: 3, 5, 6 :( :(
#Remove the lowered pH treatment tanks: 17, 20, 21, 24, 18, 19, 22, 23
#Estimate shell weight (ShW) & tissue weight (TiW) based on the following submerged regressions for each population where x is SW (submerged weight):
#Pruth	y = 1.5889x + 0.1392
#Kwakshua	y = 1.5958x + 0.0646
#Cedar	y = 1.61x + 0.0266
#Heron	y = 1.6104x - .1292
#Calculate TiW (tissue weight) based on Shell weight and Total weight
#Remove any rows with NAs in the L, T, SG, TiW & ShW columns (these died, but I kept them in my original datasheet)
#Separate the Treat into Temp & pH columns

pruth_reg <- function(x){
  ShW <- 1.5889*x + 0.1392
  return(ShW)
}
kwak_reg <- function(x){
  ShW <- 1.5958*x + 0.0646
  return(ShW)
}
cedar_reg <- function(x){
  ShW <- 1.61*x + 0.0266
  return(ShW)
}
heron_reg <- function(x){
  ShW <- 1.6104*x - .1292
  return(ShW)
}

tanks_remove <- c(3, 5, 6, 17, 20, 21, 24, 18, 19, 22, 23)

meso_clean <- meso_base %>% 
  filter(!Tank %in% tanks_remove) %>% 
  filter(!is.na(Tank)) %>% 
  rename(L = Length, Th = Thickness, SG = Linear_shell_growth, TW = Total_weight, SW = Submerged_weight) %>% 
  mutate(ShW = ifelse(SP == "Cedar", cedar_reg(SW), 
                      ifelse(SP == "Heron", heron_reg(SW),
                             ifelse(SP == "Kwak", kwak_reg(SW),
                                    ifelse(SP == "Pruth", pruth_reg(SW), NA)))),
         SR = as.factor(ifelse(SP == "Cedar" | SP == "Heron", "Strait of Georgia", "Central Coast")),
         SP = as.factor(SP),
         Stage = as.factor(Stage),
         Treat = as.factor(ifelse(Tank == 9 | Tank == 12, "12",
                                  ifelse(Tank == 4 | Tank == 7 | Tank == 11 | Tank == 14 | Tank == 16, "15",
                                         ifelse(Tank == 1 | Tank == 10 | Tank == 13, "19",
                                                ifelse(Tank == 2 | Tank == 8 | Tank == 15, "22", NA))))), 
         TiW = TW-ShW) %>% 
  filter_at(vars(L,Th, SG, ShW, TiW), any_vars(!is.na(.))) %>% 
  select(Stage, SR, SP, Tank, Treat, ID, L, Th, ShW, TiW, SG)

#You have to add the TiW & ShW values for each ID to the init ID so that you only have 1 row/snail (otherwise it messes up the sample sizes later
meso_clean_i <- meso_clean %>% 
  filter(Stage == "Init") %>% 
  select(!c(TiW, ShW))
meso_clean_m <- meso_clean %>% 
  filter(Stage == "Mid") %>% 
  select(ID, TiW, ShW)
meso_clean_f <- meso_clean %>% 
  filter(Stage == "Final") %>% 
  rename(L_final = L, Th_final = Th, ShW_final = ShW, TiW_final = TiW) %>% 
  select(Stage, ID, L_final, Th_final, ShW_final, TiW_final, SG)

meso_clean_1 <- left_join(meso_clean_i, meso_clean_m, by = "ID") %>% 
  left_join(meso_clean_f, by = "ID") %>% 
  droplevels

#For some reason CB44 was duplicated twice and has 2 incomplete rows --> remove rows 64-66
meso_clean_1 <- meso_clean_1[-c(64:66), ]

#Now we have a dataframe with rows for all the right IDs, but we need to switch it into a long format keeping the IDs. 
#SO: split it into 2 dataframes again, and then rbind them
meso_clean_2 <- meso_clean_1 %>% 
  select(SR, SP, Tank, Treat, ID, Stage.y, L_final, Th_final, ShW_final, TiW_final, SG.y) %>% 
  rename(L = L_final, Th = Th_final, ShW = ShW_final, TiW = TiW_final, Stage = Stage.y, SG = SG.y) %>%
  mutate(Alive = ifelse(is.na(Stage), 0, 1),
         Stage = ifelse(is.na(Stage), "Final", "Final")) %>% 
  select(Stage, SR, SP, Tank, Treat, ID, L, Th, ShW, TiW, SG, Alive)

meso_clean <- meso_clean_1 %>% 
  select(Stage.x, SR, SP, Tank, Treat, ID, L, Th, ShW, TiW, SG.x) %>% 
  rename(Stage = Stage.x, SG = SG.x) %>% 
  mutate(Alive = 1) %>% 
  rbind(meso_clean_2)

#Clean the feeding rate data
meso_food_clean <-meso_feed %>% 
  filter(!Tank %in% tanks_remove) %>% 
  filter(!is.na(Tank)) %>% 
  select(!c(Notes, Treatment)) %>% 
  mutate(SP = as.factor(ifelse(SP == "C", "Cedar",
                               ifelse(SP == "H", "Heron", 
                                      ifelse(SP == "K", "Kwak",
                                             ifelse(SP == "P", "Pruth", NA))))),
         Date = as.Date(Date, format = "%m-%d-%Y"),
         tank_sp = paste(Tank, SP, sep = "_"),
         Treat = as.factor(ifelse(Tank == 9 | Tank == 12, "12",
                                  ifelse(Tank == 4 | Tank == 7 | Tank == 11 | Tank == 14 | Tank == 16, "15",
                                         ifelse(Tank == 1 | Tank == 10 | Tank == 13, "19",
                                                ifelse(Tank == 2 | Tank == 8 | Tank == 15, "22", NA)))))) %>% 
  mutate(SR = ifelse(SP == "Cedar" | SP == "Heron", "Strait of Georgia", "Central Coast")) %>% 
  select(Date, Tank, SR, SP, tank_sp, Treat, Per_cap)

#Load packages----
pkgs <- c("rTPC", "nls.multstart", "broom", "tidyverse", "cowplot", "ggrepel")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)
library(minpack.lm)
library(gridExtra)

# FR only ----------------- 

# Create an empty list to store plots
plots <- list()

topt_ctmax_results = data.frame(param = vector(), temp = vector(), conf_lower = vector(), 
                                conf_upper = vector(), SR = vector(), model = vector(), period = vector(),
                                RV = vector())

for (i in 1:length(unique(meso_food_clean$Date))) {
  CC_fr <- meso_food_clean %>% 
    filter(SR == "Central Coast" & Date == unique(meso_food_clean$Date)[i]) %>% 
    select(Treat, Per_cap, SR) %>% 
    rename("temp" = "Treat", "rate" = "Per_cap")
  CC_fr$temp = as.numeric(as.character(CC_fr$temp))
  
  SoG_fr <- meso_food_clean %>% 
    filter(SR == "Strait of Georgia" & Date == unique(meso_food_clean$Date)[i]) %>% 
    select(Treat, Per_cap, SR) %>% 
    rename("temp" = "Treat", "rate" = "Per_cap")
  SoG_fr$temp = as.numeric(as.character(SoG_fr$temp))
  
  # fit five chosen model formulations in rTPC
  d_fits <- nest(CC_fr, data = c(temp, rate)) %>%
    mutate(briere = map(data, ~nls_multstart(rate~briere2_1999(temp, tmin, tmax, a, b),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
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
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', briere:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(CC_fr$temp), max(CC_fr$temp), length.out = 100))
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
  #   geom_point(aes(temp, rate), CC_fr) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Central Coast') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #CC_fr: Now start the AIC process
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
  cc_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), CC_fr) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Feeding rate',
         title = 'Central Coast') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(CC_fr, aes(temp, rate)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Rate')
  
  #CC_fr: Quadratic: Fit data
  # fit with Gaussian model
  d_fit <- nest(CC_fr, data = c(temp, rate)) %>%
    mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                                data = .x,
                                                iter = c(4,4,4),
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
           # create new temperature data
           new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
           # predict over that data,
           preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))
  
  # unnest predictions
  d_preds_CC <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Central Coast")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
  #   geom_point(aes(temp, rate), CC_fr, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Central Coast')
  
  #CC_fr: Quadratic: refit model using nlsLM
  fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                                 data = CC_fr,
                                 start = coef(d_fit$quadratic[[1]]),
                                 lower = get_lower_lims(CC_fr$temp, CC_fr$rate, model_name = 'quadratic_2008'),
                                 upper = get_upper_lims(CC_fr$temp, CC_fr$rate, model_name = 'quadratic_2008'),
                                 weights = rep(1, times = nrow(CC_fr)))
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  #CC_fr: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(CC_fr$temp), max(CC_fr$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = quadratic_2008(temp = temp, a, b, c))
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_CC <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Central Coast")
  
  #CC_fr: quadratic: Estimate parameters & CI intervals 
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
                                                            period = i, RV = "FR"))
  
  topt_ctmax_results = rbind(topt_ctmax_results, data.frame(param = ci_params_select_CC_fr$param[2],
                                                            temp = ci_params_select_CC_fr$estimate[2], 
                                                            conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                            conf_upper = ci_params_select_CC_fr$conf_upper[2],
                                                            SR = ci_params_select_CC_fr$SR[2],
                                                            model = ci_params_select_CC_fr$model[2],
                                                            period = i, RV = "FR"))
  
  # fit five chosen model formulations in rTPC
  d_fits <- nest(SoG_fr, data = c(temp, rate)) %>%
    mutate(briere = map(data, ~nls_multstart(rate~briere2_1999(temp, tmin, tmax, a, b),
                                             data = .x,
                                             iter = c(4,4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'briere2_1999') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'briere2_1999'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
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
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)))
  
  # stack models
  d_stack <- select(d_fits, -data) %>%
    pivot_longer(., names_to = 'model_name', values_to = 'fit', briere:quadratic)
  
  # get predictions using augment
  newdata <- tibble(temp = seq(min(SoG_fr$temp), max(SoG_fr$temp), length.out = 100))
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
  #   geom_point(aes(temp, rate), SoG_fr) +
  #   theme_bw(base_size = 12) +
  #   theme(legend.position = 'none') +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding',
  #        title = 'Strait of Georgia') +
  #   geom_hline(aes(yintercept = 0), linetype = 2) +
  #   scale_color_brewer(type = 'qual', palette = 2)
  
  #SoG_fr: Now start the AIC process
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
  SoG_best_fr <- ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(group = model_name), col = 'grey50', alpha = 0.5) +
    geom_line(data = filter(d_preds, model_name == best_model), col = col_best_mod) +
    geom_label_repel(aes(temp, .fitted, label = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', data = filter(d_labs, model_name == best_model), col = col_best_mod) +
    geom_point(aes(temp, rate), SoG_fr) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Feeding rate',
         title = 'Strait of Georgia') +
    geom_hline(aes(yintercept = 0), linetype = 2) 
  
  #Visualize the data
  # ggplot(SoG_fr, aes(Treat, L_growth)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Rate')
  
  #SoG_fr: Quadratic: Fit data
  # fit with Gaussian model
  d_fit <- nest(SoG_fr, data = c(temp, rate)) %>%
    mutate(quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                                data = .x,
                                                iter = c(4,4,4),
                                                start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 10,
                                                start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 10,
                                                lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                                supp_errors = 'Y',
                                                convergence_count = FALSE)),
           # create new temperature data
           new_data = map(data, ~tibble(temp = seq(min(.x$temp), max(.x$temp), length.out = 100))),
           # predict over that data,
           preds =  map2(quadratic, new_data, ~augment(.x, newdata = .y)))
  
  # unnest predictions
  d_preds_SoG <- select(d_fit, preds) %>%
    unnest(preds) %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot data and predictions
  # ggplot() +
  #   geom_line(aes(temp, .fitted), d_preds_SoG, col = 'blue') +
  #   geom_point(aes(temp, rate), SoG_fr, size = 2, alpha = 0.5) +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Feeding rate',
  #        title = 'Strait of Georgia')
  
  #SoG_fr: Quadratic: refit model using nlsLM
  fit_nlsLM <- minpack.lm::nlsLM(rate~quadratic_2008(temp = temp, a, b, c),
                                 data = SoG_fr,
                                 start = coef(d_fit$quadratic[[1]]),
                                 lower = get_lower_lims(SoG_fr$temp, SoG_fr$rate, model_name = 'quadratic_2008'),
                                 upper = get_upper_lims(SoG_fr$temp, SoG_fr$rate, model_name = 'quadratic_2008'),
                                 weights = rep(1, times = nrow(SoG_fr)))
  
  # bootstrap using case resampling
  boot1 <- Boot(fit_nlsLM, method = 'case')
  
  # look at the data
  # head(boot1$t)
  
  # hist(boot1, layout = c(2,2))
  
  #SoG_fr: quadratic: Now plot the bootstrapped models
  #create predictions of each bootstrapped model
  boot1_preds <- boot1$t %>%
    as.data.frame() %>%
    drop_na() %>%
    mutate(iter = 1:n()) %>%
    group_by_all() %>%
    do(data.frame(temp = seq(min(SoG_fr$temp), max(SoG_fr$temp), length.out = 100))) %>%
    ungroup() %>%
    mutate(pred = quadratic_2008(temp = temp, a, b, c))
  
  # calculate bootstrapped confidence intervals
  boot1_conf_preds_SoG <- group_by(boot1_preds, temp) %>%
    summarise(conf_lower = quantile(pred, 0.025),
              conf_upper = quantile(pred, 0.975)) %>%
    ungroup() %>% 
    mutate(SR = "Strait of Georgia")
  
  # plot --------
  p = ggplot() +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_fr, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = CC_fr, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_CC, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_CC,  alpha = 0.3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_fr, fun=mean, geom="point", size = 3) +
    stat_summary(aes(y = rate, x = temp, col = SR), data = SoG_fr, fun.data = "mean_se", geom = "errorbar", width = 0.2, size = 0.5) +
    geom_line(aes(temp, .fitted, col = SR), d_preds_SoG, linewidth = 1) +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper, fill = SR), boot1_conf_preds_SoG,  alpha = 0.3) +
    scale_colour_manual(values = c("skyblue", "coral")) +
    scale_fill_manual(values = c("skyblue3", "coral3")) +
    labs(x = 'Temperature (ºC)',
         y = if (i == 1 | i == 5) "Per capita weekly feeding rate" else NULL,
         col = "Source Region",
         fill = "Source Region",
         title = paste0('TPC on ', unique(meso_food_clean$Date)[i])) + 
    theme_cowplot(16) + 
    scale_x_continuous(breaks = c(10,12,14,16,18,20,22)) +
    expand_limits(x = c(10,24)) +
    scale_y_continuous(breaks = c(0,1,2)) +
    expand_limits(y = c(-0.5,2)) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text  = element_text(size = 16),
          axis.title = element_text(size = 18),
          axis.title.x = if (i == 1 | i == 2 | i == 3 | i == 4) element_blank() else NULL,)
  
  # Store the plot in the list
  plots[[i]] <- p
  
  #SoG_fr: quadratic: Estimate parameters & CI intervals
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
                                                            period = i, RV = "FR"))
  
  topt_ctmax_results = rbind(topt_ctmax_results, data.frame(param = ci_params_select_SoG_fr$param[2],
                                                            temp = ci_params_select_SoG_fr$estimate[2], 
                                                            conf_lower = ci_params_select_SoG_fr$conf_lower[2], 
                                                            conf_upper = ci_params_select_SoG_fr$conf_upper[2],
                                                            SR = ci_params_select_SoG_fr$SR[2],
                                                            model = ci_params_select_SoG_fr$model[2],
                                                            period = i, RV = "FR"))
}

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/Figures/Beaty_results")
ggsave(filename="FR_TPCs_Beaty.png", height=8, width=15, 
       plot=grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],
                         ncol = 4))
write.csv(topt_ctmax_results, "Topt_and_CTmax_FR.csv")
topt_ctmax_results = read.csv("Topt_and_CTmax_FR.csv")

# Central Coast CTmax changes ----------

ctmax_results_CC = topt_ctmax_results[which(topt_ctmax_results$SR == "Central Coast" &
                                              topt_ctmax_results$param == "ctmax"),]
ctmax_results_CC = ctmax_results_CC[3:8,] # filter out first 2 periods cos the TPCs don't look good (probably due to the heatwave)
ctmax_results_CC$date = unique(meso_food_clean$Date)[3:8]
ctmax_plots = ggplot(ctmax_results_CC, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Data/Figures/Beaty_results")
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
start_date <- as.Date("2018-07-30")
end_date <- as.Date("2018-09-05")
ts2 = ts %>%
  filter(t >= start_date & t <= end_date)

ctmax_plots_CC = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+13), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_CC$temp, x=ctmax_results_CC$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-13, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "skyblue3",
                                  size = 35)) +
  labs(x = "Date") + ggtitle("Central Coast")

max(ts2$seas)-min(ts2$seas)

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  ctmax = ctmax_results_CC$temp,
  date  = ctmax_results_CC$date
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(ctmax ~ date, data = dat)
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

# ggsave(filename="FR_CTmax_changes_CC.png", height=8, width=15, 
#        plot=ctmax_plots_CC)

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
  geom_vline(xintercept = 0, linetype = "dotted", color = "red", size=1) +
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
CC_lighthouse_weekly_95th_temp = readRDS("CC_lighthouse_weekly_95th_temp.Rds")
CC_lighthouse_weekly_95th_temp = CC_lighthouse_weekly_95th_temp[3:8,] # discard weeks 1 and 2 cos of heatwave

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
  geom_line(data = ts2, aes(x = t, y = seas), color="red", size = 2) +
  geom_smooth(aes(y=CC_TSM_df$TSM, x=CC_TSM_df$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = "TSM (°C)", sec.axis = sec_axis( trans=~.*1, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red")) +
  labs(x = "Date") 

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  TSM = CC_TSM_df$TSM,
  date  = CC_TSM_df$date
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
  geom_vline(xintercept = 0, linetype = "dotted", color = "red", size=1) +
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
ctmax_results_SoG = ctmax_results_SoG[3:8,] # filter out first 2 periods cos the TPCs don't look good (probably due to the heatwave)
ctmax_results_SoG$date = unique(meso_food_clean$Date)[3:8]
ctmax_plots = ggplot(ctmax_results_SoG, aes(date, temp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

# write.csv(ctmax_results_SoG, "FR_CTmax_changes_SoG.csv")
ctmax_results_SoG = read.csv("FR_CTmax_changes_SoG.csv")
ctmax_results_SoG$date = as.Date(ctmax_results_SoG$date)

# add seasonally varying climatology to the plot
OISST_of_interest = data.frame(t = oisst_dates,
                               temp = sst[,2])

ts = ts2clm(OISST_of_interest, climatologyPeriod = c("1983-01-01", "2012-12-31"))
start_date <- as.Date("2018-07-30")
end_date <- as.Date("2018-09-05")
ts2 = ts %>%
  filter(t >= start_date & t <= end_date)

ctmax_plots_SoG = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+8), color="red", size = 2) +
  geom_smooth(aes(y=ctmax_results_SoG$temp, x=ctmax_results_SoG$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = expression(CT[max] ~ "(°C)"), sec.axis = sec_axis( trans=~.-8, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "coral3",
                                  size = 35)) +
  labs(x = "Date") + ggtitle("Strait of Georgia")

max(ts2$seas)-min(ts2$seas)

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  ctmax = ctmax_results_SoG$temp,
  date  = ctmax_results_SoG$date
)
# 2) Fit linear model: TSM change ~ seasonal temp change
m  <- lm(ctmax ~ date, data = dat)
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

# ggsave(filename="FR_CTmax_changes_SoG.png", height=8, width=15, 
#        plot=ctmax_plots_SoG)

# save both in one plot
# ggsave(filename="FR_CTmax_changes.png", height=10, width=15, 
#        plot=grid.arrange(ctmax_plots_CC,
#                          ctmax_plots_SoG,nrow=2))

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
  geom_vline(xintercept = 0, linetype = "dotted", color = "red", size=1) +
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
SoG_lighthouse_weekly_95th_temp = readRDS("SoG_lighthouse_weekly_95th_temp.Rds")
SoG_lighthouse_weekly_95th_temp = SoG_lighthouse_weekly_95th_temp[3:8,] # discard weeks 1 and 2 cos of heatwave

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
  geom_line(data = ts2, aes(x = t, y = seas-11), color="red", size = 2) +
  geom_smooth(aes(y=SoG_TSM_df$TSM, x=SoG_TSM_df$date), method = "lm", se = TRUE, color = "blue") +
  scale_y_continuous(name = "TSM (°C)", sec.axis = sec_axis( trans=~.+11, name="Temperature (°C)")) +
  theme(axis.title.y.right = element_text(color = "red"),
        axis.text.y.right  = element_text(color = "red")) +
  labs(x = "Date") 

# test whether the slope is significantly positive or negative
# 1) Build a clean data frame (keeps x and y aligned)
dat <- data.frame(
  TSM = SoG_TSM_df$TSM,
  date  = SoG_TSM_df$date
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
  geom_vline(xintercept = 0, linetype = "dotted", color = "red", size=1) +
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
