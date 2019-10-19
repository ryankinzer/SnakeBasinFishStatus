# Author: Ryan N. Kinzer
# Purpose: Produce report map based figures for results.
# Created: 9/30/19
# Last Modified: 
# Notes: 

#-----------------------------------------------------------------

# load needed libraries
library(tidyverse)
library(sf)

# load sites

load('./data/ConfigurationFiles/site_config.rda')

config <- configuration %>%
  filter(SiteID %in% site_df$SiteID) %>%
  #filter(SiteID != 'GRA') %>%
  select(SiteID, SiteType, SiteName, RKM, Latitude, Longitude, Node) %>%
  distinct() %>%
  sf::st_as_sf(coords = c('Longitude', 'Latitude'),
               crs = 4326)

# load polygons and lines
load('./data/ConfigurationFiles/Snake_POP_metadata.rda')

stream_sp <- as_Spatial(SR_streams)
sites_sp <- as_Spatial(config)

tmp <- maptools::snapPointsToLines(sites_sp, stream_sp)
tmp2 <- st_as_sf(tmp)

stream_d <- SR_streams %>%
  group_by(GNIS_ID, GNIS_NAME) %>%
  summarise_at(vars(LENGTHKM, AreaSqKM, TotDASqKM), sum) %>%
  filter(!is.na(GNIS_NAME))

large_st <- filter(stream_d, AreaSqKM > 50)

#tmp3 <- st_intersection(stream_d, tmp2)

# pts <- st_nearest_points(config, SR_streams)
# 
# try(st_nearest_feature(config,SR_streams))
# tmp <- st_nearest_points(config, SR_streams[st_nearest_feature(config,SR_streams),], pairwise = TRUE)
# 
# tmp <- SR_streams[st_nearest_feature(config,SR_streams),] 
# GNIS <- na.omit(unique(as.character(tmp$GNIS_NAME)))
# 
# tmp2 <- st_nearest_points(config, tmp)




# tmp2 <- tmp %>%
#   st_sf() %>%
#   st_cast()

#snap_pts <- st_cast(pts, "POINT")[2] 

#pts <- st_snap(stream_d, config)
# tmp <- stream_d %>%
#   st_join(config, left = FALSE)

ggplot() +
  geom_sf(data = SR_st_pop, aes(fill = MPG)) +
  geom_sf(data = large_st, colour = 'blue') +
  geom_sf(fill = 'black', data = config) +
  coord_sf()

library(leaflet)
library(htmlwidgets)
library(htmltools)

# Title for leaflet maps
tag.map.title <- tags$style(HTML("
  .leaflet-control.map-title { 
    transform: translate(-50%,20%);
    position: fixed !important;
    left: 50%;
    text-align: center;
    padding-left: 10px; 
    padding-right: 10px; 
    background: rgba(255,255,255,0.75);
    font-weight: bold;
    font-size: 16px;
  }
"))

# Steelhead
st_title <- tags$div(
  tag.map.title, HTML("ICTRT Steelhead Populations and DABOM Sites")
)  

st_copop <- colorFactor(palette = 'viridis', domain = SR_st_pop$POP_NAME)

steelhead_sites <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addPolygons(data = SR_st_pop, group = "polys",
              popup = ~as.character(TRT_POPID), layerId = SR_st_pop$POP_NAME,
              stroke = TRUE, color = 'black',
              fillColor = ~st_copop(POP_NAME), weight = 1, fillOpacity = .5) %>%
  addPolylines(data = large_st, weight = 1) %>%
  addCircleMarkers(data = config, radius = 5, popup = ~SiteID,
                   popupOptions = popupOptions(noHide = T)) %>%
  addControl(st_title, position = 'topleft', className = "map-title")

saveWidget(steelhead_sites, file="steelhead_sites.html")
mapshot(steelhead_sites, file = 'steelhead_sites.png')

# Chinook
ch_title <- tags$div(
  tag.map.title, HTML("ICTRT Chinook Populations and DABOM Sites")
)

ch_copop <- colorFactor(palette = 'viridis', domain = SR_ch_pop$POP_NAME)

chinook_sites <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addPolygons(data = SR_ch_pop, group = "polys",
              popup = ~as.character(TRT_POPID), layerId = SR_ch_pop$POP_NAME,
              stroke = TRUE, color = 'black',
              fillColor = ~ch_copop(POP_NAME), weight = 1, fillOpacity = .5) %>%
  addPolylines(data = large_st, weight = 1) %>%
  addCircleMarkers(data = config, radius = 5, popup = ~SiteID,
                   popupOptions = popupOptions(noHide = T)) %>%
  addControl(ch_title, position = 'topleft', className = "map-title")

saveWidget(chinook_sites, file="chinook_sites.html")

# %>%
#   addCircles(SR_array_sites$Longitude, SR_array_sites$Latitude, group = 'sites', radius = 100,
#              color = 'blue', 
#              popup = paste("<b>Site ID:</b>", SR_array_sites$SiteID, "<br>"),
#              popupOptions = popupOptions(noHide = T, textsize = "15px"),
#              highlightOptions = highlightOptions(color = "white",
#                                                  weight = 5, bringToFront = F, opacity = 1))



library(mapview)
mapView(SR_st_pop)




# load function and assign group variables to all sites; TRT POPs and GSI
source('./R/assign_POP_GSI.R')

# set species
spp = 'Steelhead'
# set year
yr = 2018

# mapping data
pop_ls <- assign_POP_GSI(species = spp, configuration, site_df)
grp_df <- pop_ls[[1]]
map_df <- pop_ls[[2]]

#site_labs <- distinct(grp_df, SiteID)

# Sex....
# left_join(map_df, modSexDf, by = c('TRT_POPID' = 'TRT')) %>%
#  ggplot() +
#    geom_sf(aes(fill = propF)) +
#    scale_fill_continuous(type = 'viridis') +
# #   geom_sf(data = grp_df) +
#    geom_sf_text(aes(label = round(propF,2)), size = 3)
