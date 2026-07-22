rm(list = ls())

library(dplyr) # For basic data manipulation
library(ggplot2) # For visualising data
library(heatwaveR) # For detecting MHWs
library(tidync) # For easily dealing with NetCDF data
library(gridExtra)
library(geosphere)
library(grid)
library(maps)
library(lubridate)
library(ncdf4)
library(tidyr)
white_spaces = "                      " # 22 
white_spaces2 = "                   " # for the plots in second column: 20 

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/")
cum_int_nc <- nc_open("BC_2021MHWs_local_adaptation_model.nc")
cum_int_NCM_nc <- nc_open("BC_2021MHWs_niche_conservatism_model.nc")
cum_int_FT_nc <- nc_open("BC_2021MHWs_22deg_fixed_thresh.nc")
mean_SST_nc <- nc_open("BC_mean_SST_field.nc")

# double mhws[lon,lat]   (Contiguous storage)
mhw_cum_int <- ncvar_get(cum_int_nc, "mhws")
mhw_cum_int_NCM <- ncvar_get(cum_int_NCM_nc, "mhws")
mhw_cum_int_FT <- ncvar_get(cum_int_FT_nc, "mhws")
lat = ncvar_get(cum_int_nc, varid = "lat")
lon = ncvar_get(cum_int_nc, varid = "lon")
lon_lat <- expand.grid(lon, lat)

mhw_cum_int_df <- data.frame(cbind(lon_lat, as.vector(mhw_cum_int[,])))
colnames(mhw_cum_int_df)[1] = "lon"
colnames(mhw_cum_int_df)[2] = "lat"
colnames(mhw_cum_int_df)[3] = "cum_intensity"

mhw_cum_int_NCM_df <- data.frame(cbind(lon_lat, as.vector(mhw_cum_int_NCM[,])))
colnames(mhw_cum_int_NCM_df)[1] = "lon"
colnames(mhw_cum_int_NCM_df)[2] = "lat"
colnames(mhw_cum_int_NCM_df)[3] = "cum_intensity"

mhw_cum_int_FT_df <- data.frame(cbind(lon_lat, as.vector(mhw_cum_int_FT[,])))
colnames(mhw_cum_int_FT_df)[1] = "lon"
colnames(mhw_cum_int_FT_df)[2] = "lat"
colnames(mhw_cum_int_FT_df)[3] = "cum_intensity"

#########################
### World Map for MHW ###
#########################

theme_JMS <- function () {
  theme_bw(base_size = 12) + theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 13, face = "plain"),
    legend.key.size = unit(0.3, 'cm'),
    legend.text = element_text(size = 15, face = "plain"),
    legend.title = element_text(size = 13, face = "bold")
  )
}

### Recenter ####

center <- 200 # positive values only - US centered view is 260

# shift coordinates to recenter loggers
mhw_cum_int_df$long.recenter <-  ifelse(mhw_cum_int_df$lon < center - 180 , 
                                          mhw_cum_int_df$lon + 360, mhw_cum_int_df$lon) 
mhw_cum_int_NCM_df$long.recenter <-  ifelse(mhw_cum_int_NCM_df$lon < center - 180 , 
                                            mhw_cum_int_NCM_df$lon + 360, mhw_cum_int_NCM_df$lon) 
mhw_cum_int_FT_df$long.recenter <-  ifelse(mhw_cum_int_FT_df$lon < center - 180 , 
                                           mhw_cum_int_FT_df$lon + 360, mhw_cum_int_FT_df$lon) 

# shift coordinates to recenter world
library(ggmap)
world <- map_data("world")
world$long.recenter <-  ifelse(world$long  < center - 180 , world$long + 360, world$long)

### Function to regroup split lines and polygons
# takes dataframe, column with long and unique group variable, returns df with added column named group.regroup
RegroupElements <- function(df, longcol, idcol){  
  g <- rep(1, length(df[,longcol]))
  if (diff(range(df[,longcol])) > 300) {          # check if longitude within group differs more than 300 deg, ie if element was split
    d <- df[,longcol] > mean(range(df[,longcol])) # we use the mean to help us separate the extreme values
    g[!d] <- 1     # some marker for parts that stay in place (we cheat here a little, as we do not take into account concave polygons)
    g[d] <- 2      # parts that are moved
  }
  g <-  paste(df[, idcol], g, sep=".") # attach to id to create unique group variable for the dataset
  df$group.regroup <- g
  df
}

### Function to close regrouped polygons
# takes dataframe, checks if 1st and last longitude value are the same, if not, inserts first as last and reassigns order variable
ClosePolygons <- function(df, longcol, ordercol){
  if (df[1,longcol] != df[nrow(df),longcol]) {
    tmp <- df[1,]
    df <- rbind(df,tmp)
  }
  o <- c(1: nrow(df))  # rassign the order variable
  df[,ordercol] <- o
  df
}

library(plyr)
# now regroup
worldmap.rg <- ddply(world, .(group), RegroupElements, "long.recenter", "group")

# close polys
worldmap.cp <- ddply(worldmap.rg, .(group.regroup), ClosePolygons, "long.recenter", "order")  # use the new grouping var

# 1001 Steps Park (49° 01' 53" N 122° 52' 33" W) -> 49.031389, -122.875833
# Waterloo (49° 16' 22" N 123° 10' 40" W) -> 49.272778, -123.177778
coords_df = data.frame(lat = c(49.031389, 49.272778),
                       lon = c(-122.875833, -123.177778))
coords_df$long.recenter <-  ifelse(coords_df$lon < center - 180 , 
                                   coords_df$lon + 360, coords_df$lon) 

map_mhws <- function(df, panel_no, plot_title, legend_or_not) {
  MAP<-ggplot() + theme_JMS() + 
    geom_raster(data = df, 
                aes(x = long.recenter, y = lat, 
                    fill = cum_intensity), interpolate = TRUE) +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Warm-season\ncumulative\nintensity (°C-day)") +
    geom_polygon(aes(long.recenter,lat,group=group.regroup), size = 0.1, 
                 fill="lightgray", colour = "lightgray", data=worldmap.cp) + theme_bw() +
    scale_size(guide = "none") + 
    theme(panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          size = 0.5, linetype = "solid"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title.x=element_blank(), 
          axis.title.y=element_blank(), 
          plot.title = element_text(hjust = 0.5),
          legend.position = if (isTRUE(legend_or_not)) "right" else "none",
          text = element_text(size = 25),
          plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt"),
          legend.key.size = unit(1.5, "cm"))+ 
    scale_y_continuous(expand = c(0,0), limits = c(-60, 90))+
    scale_x_continuous(expand = c(0,0), limits = c(20, 380)) +
    ggtitle(bquote(~bold(~.(panel_no)) ~.(white_spaces) 
                   ~bold(~.(plot_title)) ~.(white_spaces))) +
    coord_cartesian(
      xlim = c(231, 238),   # converted -129 to -122
      ylim = c(47, 53)
    ) +
    geom_point(data = coords_df[1,], aes(x = long.recenter, y = lat), colour = "orange", size = 7) +
    geom_point(data = coords_df[2,], aes(x = long.recenter, y = lat), colour = "blue", size = 7)
  MAP
}

### MHW map
mhw_cum_int_df$cum_intensity[which(mhw_cum_int_df$cum_intensity == "NaN")] = NA
mhws_LAM = map_mhws(mhw_cum_int_df, "", "Statistically determined thresholds", TRUE) 
mhw_cum_int_df_filtered = mhw_cum_int_df[which(!is.na(mhw_cum_int_df$cum_intensity)),]
max(mhw_cum_int_df_filtered$cum_intensity) # 260.4554

mhw_cum_int_NCM_df$cum_intensity[which(mhw_cum_int_NCM_df$cum_intensity == "NaN")] = NA
mhws_NCM = map_mhws(mhw_cum_int_NCM_df, "", "Biologically determined thresholds", TRUE) 
mhw_cum_int_NCM_df_filtered = mhw_cum_int_NCM_df[which(!is.na(mhw_cum_int_NCM_df$cum_intensity)),]

mhw_cum_int_FT_df$cum_intensity[which(mhw_cum_int_FT_df$cum_intensity == "NaN")] = NA
mhws_FT = map_mhws(mhw_cum_int_FT_df, "", "Fixed threshold 22degC", FALSE) 
mhw_cum_int_FT_df_filtered = mhw_cum_int_FT_df[which(!is.na(mhw_cum_int_FT_df$cum_intensity)),]

mhw_cum_int_diff_df = mhw_cum_int_df
mhw_cum_int_diff_df$cum_intensity = mhw_cum_int_df$cum_intensity - mhw_cum_int_NCM_df$cum_intensity
mhws_diff = map_mhws(mhw_cum_int_diff_df, "", "LAM - NCM", FALSE) 

# save
setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub/field_mortality_2021_heat_dome")
ggsave(filename="LSM_vs_NCM_v2.png", height=5, width=15, 
       plot=grid.arrange(mhws_LAM, mhws_NCM, nrow=1))

# add time series of tide, and air and ocean temp --------------

# define end and start dates - 1001 Steps
end_date <- as.Date("2021-07-22")
start_date <- end_date - 29

air_temp_1001_steps = read.csv("WHITE ROCK CAMPBELL SCIENTIFIC_climate-hourly.csv")
# make sure date column is Date type
air_temp_1001_steps$LOCAL_DATE <- as.POSIXct(
  air_temp_1001_steps$LOCAL_DATE,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"   # or your actual timezone
)
# filter
filtered_air_temp_1001_steps <- air_temp_1001_steps %>%
  filter(LOCAL_DATE >= start_date & LOCAL_DATE <= end_date)

air_temp_waterloo = read.csv("VANCOUVER HARBOUR CS_climate-hourly.csv")
# make sure date column is Date type
air_temp_waterloo$LOCAL_DATE <- as.POSIXct(
  air_temp_waterloo$LOCAL_DATE,
  format = "%Y-%m-%d %H:%M:%S",
  tz = "UTC"   # or your actual timezone
)
# filter
filtered_air_temp_waterloo <- air_temp_waterloo %>%
  filter(LOCAL_DATE >= start_date & LOCAL_DATE <= end_date)

ocean_temp = read.csv("ECCC_MSC_BUOYS_6823_4047_a50f.csv")
ocean_temp = ocean_temp[2:length(ocean_temp$time),]
# make sure date column is Date type
ocean_temp$time <- as.POSIXct(
  ocean_temp$time,
  format = "%Y-%m-%dT%H:%M:%S",
  tz = "UTC"   # or your actual timezone
)
# filter
filtered_ocean_temp <- ocean_temp %>%
  filter(time >= start_date & time <= end_date)

# create full sequence
full_time <- data.frame(
  time = seq(
    from = as.POSIXct("2021-06-23 00:00:00"),
    to   = as.POSIXct("2021-07-22 00:00:00"),
    by   = "10 min"
  )
)

# join and fill missing with NA
filtered_ocean_temp <- full_time %>%
  left_join(filtered_ocean_temp, by = "time")
filtered_ocean_temp = filtered_ocean_temp[61:length(filtered_ocean_temp$time),]

tide_levels_waterloo = read.csv("Kitsilano Tides/07707_data.csv")
tide_levels_waterloo = tide_levels_waterloo[7:length(tide_levels_waterloo$Station.Name...Nom.de.la.station),]
colnames(tide_levels_waterloo)[1] = "date"
colnames(tide_levels_waterloo)[2] = "water_level"
# make sure date column is Date type
tide_levels_waterloo$date <- as.POSIXct(
  tide_levels_waterloo$date,
  format = "%Y/%m/%d %H:%M",
  tz = "UTC"   # or your actual timezone
)
# filter
filtered_tide_levels_waterloo <- tide_levels_waterloo %>%
  filter(date >= start_date & date <= end_date)

tide_levels_1001_steps = read.csv("Sand Heads Tides/07594_data.csv")
tide_levels_1001_steps = tide_levels_1001_steps[7:length(tide_levels_1001_steps$Station.Name...Nom.de.la.station),]
colnames(tide_levels_1001_steps)[1] = "date"
colnames(tide_levels_1001_steps)[2] = "water_level"
# make sure date column is Date type
tide_levels_1001_steps$date <- as.POSIXct(
  tide_levels_1001_steps$date,
  format = "%Y/%m/%d %H:%M",
  tz = "UTC"   # or your actual timezone
)
# filter
filtered_tide_levels_1001_steps <- tide_levels_1001_steps %>%
  filter(date >= start_date & date <= end_date)

# plots ---------------- 

# Tide
filtered_tide_levels_1001_steps <- filtered_tide_levels_1001_steps %>%
  mutate(water_level = as.numeric(as.character(water_level)))
filtered_tide_levels_waterloo <- filtered_tide_levels_waterloo %>%
  mutate(water_level = as.numeric(as.character(water_level)))

# Ocean temp
filtered_ocean_temp <- filtered_ocean_temp %>%
  mutate(avg_sea_sfc_temp_pst10mts = as.numeric(as.character(avg_sea_sfc_temp_pst10mts)))

# Air temp
filtered_air_temp_1001_steps <- filtered_air_temp_1001_steps %>%
  mutate(TEMP = as.numeric(as.character(TEMP)))
filtered_air_temp_waterloo <- filtered_air_temp_waterloo %>%
  mutate(TEMP = as.numeric(as.character(TEMP)))

# standardise datetime columns
filtered_tide_levels_1001_steps <- filtered_tide_levels_1001_steps %>%
  mutate(datetime = as.POSIXct(date))

filtered_ocean_temp <- filtered_ocean_temp %>%
  mutate(datetime = as.POSIXct(time))

filtered_air_temp_1001_steps <- filtered_air_temp_1001_steps %>%
  mutate(datetime = as.POSIXct(LOCAL_DATE))

filtered_air_temp_waterloo <- filtered_air_temp_waterloo %>%
  mutate(datetime = as.POSIXct(LOCAL_DATE))

# # compute ranges
# tide_range <- range(tide$water_level, na.rm = TRUE)
# temp_range <- range(c(ocean$avg_sea_sfc_temp_pst10mts, air$TEMP), na.rm = TRUE)
# 
# # scaling function
# scale_factor <- diff(tide_range) / diff(temp_range)
# temp_min <- temp_range[1]
# tide_min <- tide_range[1]

# define temperature tick marks (every 2°C)
temp_breaks <- seq(
  floor(temp_range[1]),
  ceiling(temp_range[2]),
  by = 2
)

low_tide_periods <- filtered_tide_levels_1001_steps %>%
  arrange(date) %>%
  mutate(
    below = water_level < 1.75,
    group = cumsum(!below)
  ) %>%
  filter(below) %>%
  group_by(group) %>%
  reframe(
    xmin = date[1],
    xmax = date[n()]
  )

above_22degC_periods <- filtered_air_temp_1001_steps %>%
  arrange(LOCAL_DATE) %>%
  mutate(
    above = TEMP > 22,
    group = cumsum(!above)
  ) %>%
  filter(above) %>%
  group_by(group) %>%
  reframe(
    xmin = LOCAL_DATE[1],
    xmax = LOCAL_DATE[n()]
  )

library(dplyr)

overlap_periods <- above_22degC_periods %>%
  dplyr::rename(xmin1 = xmin, xmax1 = xmax) %>%
  inner_join(
    low_tide_periods %>%
      dplyr::rename(xmin2 = xmin, xmax2 = xmax),
    by = character()   # cross join
  ) %>%
  # keep only overlapping intervals
  filter(xmin1 <= xmax2 & xmax1 >= xmin2) %>%
  # compute intersection
  mutate(
    xmin = pmax(xmin1, xmin2),
    xmax = pmin(xmax1, xmax2)
  ) %>%
  select(xmin, xmax) %>%
  arrange(xmin)

temp_panel = ggplot() +
  # Ocean temp
  geom_line(data = filtered_ocean_temp,
            aes(x = as.POSIXct(time),
                y = avg_sea_sfc_temp_pst10mts,
                colour = "Ocean temp"),
            linewidth = 1) +
  # Air temp
  geom_line(data = filtered_air_temp_1001_steps,
            aes(x = as.POSIXct(LOCAL_DATE),
                y = TEMP,
                colour = "Air temp"),
            linewidth = 1) +
  # Air temp - waterloo
  geom_line(data = filtered_air_temp_waterloo,
            aes(x = as.POSIXct(LOCAL_DATE),
                y = TEMP,
                colour = "Air temp", alpha = 0.1),
            linewidth = 1) +
  geom_hline(
    yintercept = 22,
    colour = "darkgreen",
    linetype = "dashed",
    linewidth = 1
  ) +
  scale_colour_manual(values = c(
    "Ocean temp" = "blue",
    "Air temp" = "red"
  )) +
  theme_minimal() + labs(x = "Date", y = "Temperature (°C)") +
  scale_y_continuous(breaks = temp_breaks) +
  scale_x_datetime(date_breaks = "2 days", date_labels = "%d %b") +
  theme(legend.position = "none",
        text = element_text(size = 35),
        plot.background = element_rect(colour = "black", fill = NA, size = 2)) +
  geom_rect(
    data = overlap_periods,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = "red",
    alpha = 0.15
  ) 

temp_panel

tide_panel = ggplot() +
  # Tide
  geom_line(data = filtered_tide_levels_1001_steps,
            aes(x = as.POSIXct(date), y = water_level),
            linewidth = 1) +
  geom_line(data = filtered_tide_levels_waterloo,
            aes(x = as.POSIXct(date), y = water_level, alpha = 0.5),
            linewidth = 1) +
  geom_hline(
    yintercept = 1,
    colour = "orange",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 1.75,
    colour = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_hline(
    yintercept = 0.5,
    colour = "blue",
    linetype = "dashed",
    linewidth = 1
  ) +
  theme_minimal() + labs(x = "Date", y = "Water level (m)") +
  scale_x_datetime(date_breaks = "2 days", date_labels = "%d %b") +
  theme(legend.position = "none",
        text = element_text(size = 35),
        plot.background = element_rect(colour = "black", fill = NA, size = 2))

tide_panel

# mort data ---------------
mort_data = read.csv("C:/Users/dlcyli/OneDrive/2025/Paper 2/Nucella_mortality_by_shore_level.csv")

mort_data <- mort_data %>%
  mutate(group = paste(Site, Shore_level, sep = " - "))

mort_data$group[which(mort_data$group == "Waterloo - very_low")] = "Waterloo - very low"
mort_data$group[which(mort_data$group == "1001_Steps - low")] = "1001 Steps - low"
mort_data$group[which(mort_data$group == "1001_Steps - mid")] = "1001 Steps - mid"

mort_data$group <- factor(
  mort_data$group,
  levels = c(
    "Waterloo - very low",
    "1001 Steps - low",
    "1001 Steps - mid"
  )
)

mort_plot = ggplot(mort_data, aes(x = group, y = X.recent_dead)) +
  geom_jitter(
    aes(color = group),
    width = 0.1,
    height = 0,
    size = 5
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.5,
    color = "black"
  ) +
  scale_color_manual(values = c(
    "1001 Steps - low"   = "orange",
    "1001 Steps - mid"   = "red",
    "Waterloo - very low" = "blue"
  )) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text = element_text(size = 35),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10, unit = "pt"),
    plot.background = element_rect(colour = "black", fill = NA, size = 2)
  ) +
  ylab("Mortality (%)") +
  xlab("Site")

ggsave(
  filename = "2021_heat_dome.png",
  height = 25,
  width = 20,
  plot = grid.arrange(
    
    # first row (2 panels)
    arrangeGrob(mhws_LAM, mhws_NCM, ncol = 2),
    
    # next rows (1 panel each)
    temp_panel,
    tide_panel,
    mort_plot,
    
    nrow = 4,
    heights = c(1.3,1,1,0.8)
  )
)

