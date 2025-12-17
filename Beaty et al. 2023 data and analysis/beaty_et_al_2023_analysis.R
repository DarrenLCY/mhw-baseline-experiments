#This script will fit TPCs to my mesocosm snails data!! 
#Begin by running the 'snails_meso.R' script until line 320 (DLSH: changed from 290 here). You want the dataframe with columns for averaged RV w/in each tank
#This is 'meso_lm_block', which has the change in growth metrics averaged w/in tank

#Start with creating TPCs for tissue weight
#Now try with the other growth metrics (Shell length, shell weight, per capita feeding rate)

#Load packages----
pkgs <- c("rTPC", "nls.multstart", "broom", "tidyverse", "cowplot", "ggrepel")
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

# ----------------- FR only - get meso_food_clean variable only from other script

# ---------- Central Coast
# Create an empty list to store plots
plots <- list()

ctmax_results = data.frame(ctmax = vector(), conf_lower = vector(), conf_upper = vector())

for (i in 1:length(unique(meso_food_clean$Date))) {
  CC_fr <- meso_food_clean %>% 
    filter(SR == "Central Coast" & Date == unique(meso_food_clean$Date)[i]) %>% 
    select(Treat, Per_cap) %>% 
    rename("temp" = "Treat", "rate" = "Per_cap")
  CC_fr$temp = as.numeric(as.character(CC_fr$temp))
  
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
  ggplot(d_preds, aes(temp, .fitted)) +
    geom_line(aes(col = model_name)) +
    geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
    geom_point(aes(temp, rate), CC_fr) +
    theme_bw(base_size = 12) +
    theme(legend.position = 'none') +
    labs(x = 'Temperature (ºC)',
         y = 'Feeding',
         title = 'Central Coast') +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    scale_color_brewer(type = 'qual', palette = 2)
  
  #CC_fr: Now start the AIC process----
  d_ic <- d_stack %>%
    mutate(., info = map(fit, glance),
           AICc =  map_dbl(fit, MuMIn::AICc)) %>%
    select(-fit) %>%
    unnest(info) %>%
    select(model_name, sigma, AIC, AICc, BIC, df.residual)
  
  d_ic
  
  # filter for best model
  best_model = filter(d_ic, AICc == min(AICc)) %>% pull(model_name)
  best_model
  
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
  
  #Visualize the data----
  # ggplot(CC_fr, aes(Treat, Per_cap)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Rate')
  
  #CC_fr: Quadratic: Fit data----
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
  
  #CC_fr: Quadratic: refit model using nlsLM----
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
  
  #CC_fr: quadratic: Now plot the bootstrapped models----
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
  
  # plot bootstrapped CIs
  p = ggplot() +
    geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_CC, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), CC_fr, size = 2, alpha = 0.5) +
    theme_bw(base_size = 17) +
    labs(x = 'Temperature (ºC)',
         y = 'Feeding rate',
         title = paste0('TPC on ', unique(meso_food_clean$Date)[i]))
  
  # Store the plot in the list
  plots[[i]] <- p
  
  if (i > 2) {
    #CC_fr: quadratic: Estimate parameters & CI intervals ----
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
    
    ctmax_results = rbind(ctmax_results, data.frame(ctmax = ci_params_select_CC_fr$estimate[2], 
                                                    conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                    conf_upper = ci_params_select_CC_fr$conf_upper[2]))
    
    # ggplot(ci_params_select_CC_fr, aes(param, estimate)) +
    #   geom_point(size = 4) +
    #   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
    #   theme_bw() +
    #   facet_wrap(~param, scales = 'free') +
    #   scale_x_discrete('') +
    #   labs(title = 'quadratic - CC')
  }
}

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Beaty et al")
library(gridExtra)
ggsave(filename="FR_TPCs_CC.png", height=8, width=15, 
       plot=grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],
                         ncol = 4))

ctmax_results$date = unique(meso_food_clean$Date)[3:8]
ctmax_plots = ggplot(ctmax_results, aes(date, ctmax)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ggsave(filename="FR_CTmax_changes_CC.png", height=8, width=15, 
       plot=ctmax_plots)
saveRDS(ctmax_results, "FR_CTmax_changes_CC.Rds")
ctmax_results = readRDS("FR_CTmax_changes_CC.Rds")

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

# 12degC treatment
ts = ts2clm(OISST_of_interest, climatologyPeriod = c("1983-01-01", "2012-12-31"))
start_date <- as.Date("2018-07-30")
end_date <- as.Date("2018-09-05")
ts2 = ts %>%
  filter(t >= start_date & t <= end_date)

ctmax_plots = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+13), color="red", size = 2) +
  scale_y_continuous(name = "CTmax (degC)", sec.axis = sec_axis( trans=~.-13, name="Temperature (degC)")) +
  theme(axis.title.y.right = element_text(color = "red"))

ggsave(filename="FR_CTmax_changes_CC.png", height=8, width=15, 
       plot=ctmax_plots)

library(metafor)
# Calculate standard errors from confidence intervals
ctmax_results$se <- (ctmax_results$conf_upper - ctmax_results$conf_lower) / (2 * 1.96)

# Perform meta-analysis
meta_result <- rma(yi = ctmax, sei = se, data = ctmax_results)

# Summary of the meta-analysis
summary(meta_result)

# ---------- Strait of Georgia
# Create an empty list to store plots
plots <- list()

ctmax_results = data.frame(ctmax = vector(), conf_lower = vector(), conf_upper = vector())

for (i in 1:length(unique(meso_food_clean$Date))) {
  CC_fr <- meso_food_clean %>% 
    filter(SR == "Strait of Georgia" & Date == unique(meso_food_clean$Date)[i]) %>% 
    select(Treat, Per_cap) %>% 
    rename("temp" = "Treat", "rate" = "Per_cap")
  CC_fr$temp = as.numeric(as.character(CC_fr$temp))
  
  #Visualize the data----
  # ggplot(CC_fr, aes(Treat, Per_cap)) +
  #   geom_point() +
  #   theme_bw(base_size = 12) +
  #   labs(x = 'Temperature (ºC)',
  #        y = 'Rate')
  
  #CC_fr: Quadratic: Fit data----
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
  
  #CC_fr: Quadratic: refit model using nlsLM----
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
  
  #CC_fr: quadratic: Now plot the bootstrapped models----
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
  
  # plot bootstrapped CIs
  p = ggplot() +
    geom_line(aes(temp, .fitted), d_preds_CC, col = 'blue') +
    geom_ribbon(aes(temp, ymin = conf_lower, ymax = conf_upper), boot1_conf_preds_CC, fill = 'blue', alpha = 0.3) +
    geom_point(aes(temp, rate), CC_fr, size = 2, alpha = 0.5) +
    theme_bw(base_size = 17) +
    labs(x = 'Temperature (ºC)',
         y = 'Feeding rate',
         title = paste0('TPC on ', unique(meso_food_clean$Date)[i]))
  
  # Store the plot in the list
  plots[[i]] <- p
  
  if (i > 2) {
    #CC_fr: quadratic: Estimate parameters & CI intervals ----
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
    
    ctmax_results = rbind(ctmax_results, data.frame(ctmax = ci_params_select_CC_fr$estimate[2], 
                                                    conf_lower = ci_params_select_CC_fr$conf_lower[2], 
                                                    conf_upper = ci_params_select_CC_fr$conf_upper[2]))
    
    # ggplot(ci_params_select_CC_fr, aes(param, estimate)) +
    #   geom_point(size = 4) +
    #   geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
    #   theme_bw() +
    #   facet_wrap(~param, scales = 'free') +
    #   scale_x_discrete('') +
    #   labs(title = 'quadratic - CC')
  }
}

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/Beaty et al")
library(gridExtra)
ggsave(filename="FR_TPCs_SoG.png", height=8, width=15, 
       plot=grid.arrange(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],plots[[6]],plots[[7]],plots[[8]],
                         ncol = 4))

ctmax_results$date = unique(meso_food_clean$Date)[3:8]
ctmax_plots=ggplot(ctmax_results, aes(date, ctmax)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = conf_lower, ymax = conf_upper)) +
  theme_bw(base_size = 20)

ggsave(filename="FR_CTmax_changes_SoG.png", height=8, width=15, 
       plot=ctmax_plots)
saveRDS(ctmax_results, "FR_CTmax_changes_SoG.Rds")
ctmax_results = readRDS("FR_CTmax_changes_SoG.Rds")

# add seasonally varying climatology to the plot
OISST_of_interest = data.frame(t = oisst_dates,
                               temp = sst[,2])

# 12degC treatment
ts = ts2clm(OISST_of_interest, climatologyPeriod = c("1983-01-01", "2012-12-31"))
start_date <- as.Date("2018-07-30")
end_date <- as.Date("2018-09-05")
ts2 = ts %>%
  filter(t >= start_date & t <= end_date)

ctmax_plots = ctmax_plots + 
  geom_line(data = ts2, aes(x = t, y = seas+8), color="red", size = 2) +
  scale_y_continuous(name = "CTmax (degC)", sec.axis = sec_axis( trans=~.-8, name="Temperature (degC)")) +
  theme(axis.title.y.right = element_text(color = "red"))

ggsave(filename="FR_CTmax_changes_SoG.png", height=8, width=15, 
       plot=ctmax_plots)

library(metafor)
# Calculate standard errors from confidence intervals
ctmax_results$se <- (ctmax_results$conf_upper - ctmax_results$conf_lower) / (2 * 1.96)

# Perform meta-analysis
meta_result <- rma(yi = ctmax, sei = se, data = ctmax_results)

# Summary of the meta-analysis
summary(meta_result)





