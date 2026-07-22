#' ------------------------------------------------------------------------------
#' Script to plot field mortality for N. lamellosa during 2021 heat dome
#' 
#' @author  Darren Li Shing Hiung, University of Tasmania
#' @date    15 June 2026 
#' @version 1.0
#' @license MIT
#' ------------------------------------------------------------------------------
rm(list = ls())

#Load packages----
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
white_spaces = "                      " # 22 
white_spaces2 = "                   " # for the plots in second column: 20 

mort_data = read.csv("C:/Users/dlcyli/OneDrive/2025/Paper 2/Nucella_mortality_by_shore_level.csv")

setwd("C:/Users/dlcyli/OneDrive/Development of thesis/Nucella experiments/GitHub")
cum_int_NCM_nc <- nc_open("BC_2021MHWs_niche_conservatism_model.nc")

# double mhws[lon,lat]   (Contiguous storage)
mhw_cum_int_NCM <- ncvar_get(cum_int_NCM_nc, "mhws")
lat = ncvar_get(cum_int_NCM_nc, varid = "lat")
lon = ncvar_get(cum_int_NCM_nc, varid = "lon")
lon_lat <- expand.grid(lon, lat)

mhw_cum_int_NCM_df <- data.frame(cbind(lon_lat, as.vector(mhw_cum_int_NCM[,])))
colnames(mhw_cum_int_NCM_df)[1] = "lon"
colnames(mhw_cum_int_NCM_df)[2] = "lat"
colnames(mhw_cum_int_NCM_df)[3] = "cum_intensity"

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
mhw_cum_int_NCM_df$long.recenter <-  ifelse(mhw_cum_int_NCM_df$lon < center - 180 , 
                                            mhw_cum_int_NCM_df$lon + 360, mhw_cum_int_NCM_df$lon) 

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
                    fill = cum_intensity), interpolate = FALSE) +
    scale_fill_gradient(low = "white", high = "red",
                        name = "Integrated\nextreme\nintensity (°C-day)") +
    geom_polygon(aes(long.recenter,lat,group=group.regroup), size = 0.1, 
                 fill="lightgray", colour = "lightgray", data=worldmap.cp) + theme_bw() +
    scale_size(guide = "none") + 
    theme(panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          size = 0.5, linetype = "solid"),
          panel.border = element_rect(colour = "black", fill=NA, size=1),
          axis.title.x=element_blank(), 
          axis.title.y=element_blank(), 
          plot.title = element_blank(),
          legend.position = if (isTRUE(legend_or_not)) "left" else "none",
          legend.key.height = unit(2, "cm"),
          text = element_text(size = 30),
          plot.margin = margin(t = 10, r = 20, b = 10, l = 10, unit = "pt"))+ # adds box around legend and fixes spacing 
    scale_y_continuous(expand = c(0,0), limits = c(-60, 90))+
    scale_x_continuous(expand = c(0,0), limits = c(20, 380)) +
    ggtitle(bquote(~bold(~.(panel_no)) ~.(white_spaces) 
                   ~bold(~.(plot_title)) ~.(white_spaces))) +
    geom_point(data = coords_df[1,], aes(x = long.recenter, y = lat), colour = "orange", size = 5) +
    geom_point(data = coords_df[2,], aes(x = long.recenter, y = lat), colour = "blue", size = 5) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = "a",
      hjust = 2,   # push slightly left
      vjust = 1.7,   # push slightly down
      size = 15,     # adjust font size
      fontface = "bold"
    ) +
    coord_cartesian(
      xlim = c(235, 238), 
      ylim = c(48, 51)
    )
  MAP
}

### MHW map
mhw_cum_int_NCM_df_zoom <- mhw_cum_int_NCM_df %>%
  dplyr::filter(
    long.recenter >= 235 & long.recenter <= 238,
    lat >= 48 & lat <= 51
  )

mhw_cum_int_NCM_df_zoom$cum_intensity[which(mhw_cum_int_NCM_df_zoom$cum_intensity == "NaN")] = NA
mhws_NCM = map_mhws(mhw_cum_int_NCM_df_zoom, "", "Biologically determined thresholds", TRUE) 
max(mhw_cum_int_NCM_df_zoom$cum_intensity, na.rm = T) # 5.9

# box plots -------------
mort_data <- mort_data %>%
  mutate(group = paste(Site, Shore_level, sep = " - "))

mort_data$group[which(mort_data$group == "Waterloo - very_low")] = "Waterloo - very low"
mort_data$group[which(mort_data$group == "1001_Steps - low")] = "1001 Steps - low"
mort_data$group[which(mort_data$group == "1001_Steps - mid")] = "1001 Steps - mid"

box_plots = ggplot(mort_data, aes(x = group, y = X.recent_dead, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.4) +   # hide default outliers (cleaner)
  geom_jitter(
    aes(color = group),
    width = 0.15,        # spread points horizontally
    height = 0,
    size = 4
  ) +
  scale_fill_manual(values = c(
    "1001 Steps - low"   = "orange",
    "1001 Steps - mid"   = "red",
    "Waterloo - very low" = "blue"
  )) +
  scale_color_manual(values = c(
    "1001 Steps - low"   = "orange",
    "1001 Steps - mid"   = "red",
    "Waterloo - very low" = "blue"
  )) +
  theme_bw() +
  theme(legend.position = "none",
        text = element_text(size = 30)) +
  ylab("Mortality (%)") +
  xlab("Site") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "b",
    hjust = 2,   # push slightly left
    vjust = 2,   # push slightly down
    size = 15,     # adjust font size
    fontface = "bold"
  ) 

# save
setwd("C:/Users/dlcyli/OneDrive/2025/Paper 2")
ggsave(filename="field_mortality_2021_heat_dome.png", height=18, width=15, 
       plot=grid.arrange(mhws_NCM, box_plots, nrow=2))

ggplot(mort_data, aes(x = group, y = X.recent_dead)) +
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
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 30)
  ) +
  ylab("Mortality (%)") +
  xlab("Site")

# statistical tests ---------------

comparisons <- list(
  c("1001 Steps - low", "1001 Steps - mid"),
  c("1001 Steps - low", "Waterloo - very low"),
  c("1001 Steps - mid", "Waterloo - very low")
)

library(rstatix)

pairwise_results <- mort_data %>%
  wilcox_test(X.recent_dead ~ group, 
              comparisons = comparisons,
              p.adjust.method = "bonferroni")

pairwise.wilcox.test(
  mort_data$X.recent_dead,
  mort_data$group,
  p.adjust.method = "bonferroni"
)
pairwise_results
