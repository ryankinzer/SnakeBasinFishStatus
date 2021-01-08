# Author: Ryan N. Kinzer
# Purpose: Produce report map based figures for results.
# Created: 9/30/19
# Last Modified: 
# Notes: 

#-----------------------------------------------------------------

# load needed libraries
library(tidyverse)
library(sf)
library(scales)
library(ggmap)
library(maps)

# Get states
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))

# load sites and metadata

load('./data/ConfigurationFiles/site_config.rda')
site_meta <- read_csv('./data/ConfigurationFiles/site_metadata.csv')

# load polygons and lines
load('./data/ConfigurationFiles/Snake_POP_metadata.rda')

# Aggregate snake basin

SR_basin <- raster::aggregate(as_Spatial(SR_st_pop))
SR_basin <- st_as_sf(SR_basin)

# select only DABOM sites
config <- configuration %>%
  filter(SiteID %in% site_df$SiteID) %>%
  #filter(SiteID != 'GRA') %>%
  select(SiteID, SiteType, SiteTypeName, SiteName, RKM, Latitude, Longitude) %>%
  distinct() %>%
  left_join(site_meta, by = 'SiteID') %>%
  sf::st_as_sf(coords = c('Longitude', 'Latitude'),
               crs = 4326)

site_types <- unique(config$SiteTypeName)

# dissolve line segments into full streams
stream_d <- SR_streams %>%
  group_by(GNIS_ID, GNIS_NAME) %>%
  summarise_at(vars(LENGTHKM, AreaSqKM, TotDASqKM), sum) %>%
  filter(!is.na(GNIS_NAME))

# filter for maps
large_st <- filter(stream_d, AreaSqKM > 40)

# test_st <- stream_d %>%
#   filter(GNIS_NAME %in% c('Snake River', 'Salmon River', 'Tucannon River', 'Imnaha River', 'Grande Ronde River', 'South Fork Salmon River', 'Middle Fork Salmon River', 'Clearwater River', 'Potlatch River', 'Lapwai Creek', 'Lochsa River', 'Selway River', 'South Fork Clearwater River', 'Wallowa River', 'Wenaha River', 'Big Creek', 'Lemhi River', 'Lolo Creek', 'Asotin Creek', 'Pahsimeroi River', 'East Fork South Fork Salmon River', 'Panther Creek', 'North Fork Salmon River', 'Yankee Fork Salmon River', 'Little Salmon River'))

# Get Google basemap
#register_google(key = "hidden", write = TRUE)

pop_sp <- as_Spatial(SR_st_pop)
b <- sp::bbox(pop_sp) # get bounding box

map_center <- apply(b,1,mean)

p <- ggmap(get_googlemap(center = map_center, #c(lon = -122.335167, lat = 47.608013), #center = map_center,
                         zoom = 7, scale = 2,
                         maptype ='terrain',
                         color = 'color'))

#---------------------------------------------
# Map for Rick 
p <- ggmap(get_googlemap(center = "Lapwai, Idaho",
                         zoom = 10, scale = 2,
                         maptype = 'terrain',
                         color = 'color'))

arrays <- config %>% filter(DABOM_Branch == 'Lapwai')

arrays <- bind_cols(arrays, 
                    st_coordinates(arrays) %>%
                      as_tibble() %>%
                      rename(lng = X, lat = Y))

labs <- arrays %>% group_by(BPA_Funding_OM) %>% summarise(n = n()) %>% transmute(labs = paste0(BPA_Funding_OM, ' (n = ',n,')')) %>% pull(labs)

lapwai_gg <- p +
  #ggplot() +
  #geom_sf(data = states, fill = 'white') +
  geom_sf(data = SR_basin, alpha = .5, inherit.aes = FALSE) +
  #geom_sf(data = SR_st_pop, aes(fill = TRT_POPID), inherit.aes = FALSE) +
  geom_sf(data = large_st, colour = 'blue', inherit.aes = FALSE) +
  geom_point(data = arrays, aes(x = lng, y = lat), size = 3, inherit.aes = FALSE) +
  ggrepel::geom_label_repel(data = arrays, aes(x = lng, y = lat, label = SiteID), size = 6) +
  #scale_fill_viridis_d(alpha = .35, direction = -1, option = "D", guide = FALSE) +
  #scale_colour_manual(values = topo.colors(4)[1:3], labels = labs) +
  #scale_colour_viridis_d(option = "C", end = .8, labels = labs) +
  #scale_colour_brewer(palette = 'PuRd', labels = labs) +
  theme_void() +
  theme(text = element_text(family = 'serif'),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(colour = 'grey35'),
        legend.position = c(.18,.12),
        legend.background = element_rect(fill = alpha('grey75', .5))) +
  labs(title = "In-stream PIT-tag Detection Systems within the Lapwai Creek Watershed",
       subtitle = "")# +
#coord_sf(xlim = c(-112,-120), ylim = c(43,48))
lapwai_gg

ggsave(lapwai_gg, filename = './Figures/lapwai_sites_map.png', height = 9, width = 9)

#-------------------------------------------
# Steelhead POPs and O&M Sites

arrays <- config %>%
  filter(SiteTypeName == site_types[1],
         Operational)

arrays <- bind_cols(arrays, 
                    st_coordinates(arrays) %>%
                      as_tibble() %>%
                      rename(lng = X, lat = Y))

labs <- arrays %>% group_by(BPA_Funding_OM) %>% summarise(n = n()) %>% transmute(labs = paste0(BPA_Funding_OM, ' (n = ',n,')')) %>% pull(labs)

#---------------------------------------------------------------------------
# Sites with labels
#---------------------------------------------------------------------------
steelhead_gg <- p +
  #ggplot() +
  #geom_sf(data = states, fill = 'white') +
  geom_sf(data = SR_basin, alpha = .5, inherit.aes = FALSE) +
  #geom_sf(data = SR_st_pop, aes(fill = TRT_POPID), inherit.aes = FALSE) +
  geom_sf(data = large_st, colour = 'blue', inherit.aes = FALSE) +
  geom_point(data = arrays, aes(x = lng, y = lat), size = 3, inherit.aes = FALSE) +
  ggrepel::geom_label_repel(data = arrays, aes(x = lng, y = lat, label = SiteID), size = 2) +
  #scale_fill_viridis_d(alpha = .35, direction = -1, option = "D", guide = FALSE) +
  #scale_colour_manual(values = topo.colors(4)[1:3], labels = labs) +
  #scale_colour_viridis_d(option = "C", end = .8, labels = labs) +
  #scale_colour_brewer(palette = 'PuRd', labels = labs) +
  theme_void() +
  theme(text = element_text(family = 'serif'),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(colour = 'grey35'),
        legend.position = c(.18,.12),
        legend.background = element_rect(fill = alpha('grey75', .5))) +
  labs(title = "Snake River Basin In-stream PIT-tag Detection Systems",
       subtitle = "Sixty-seven IPTDS are deployed across the basin to monitor abundance,\n life history and productivity of Steelhead and spring/summer Chinook salmon.")# +
#coord_sf(xlim = c(-112,-120), ylim = c(43,48))
steelhead_gg

ggsave(steelhead_gg, filename = './Figures/IPTDS_sites_labels_map.png', height = 9, width = 9)

#-------------------------------------------------------------------
# Funding arrays.
#-------------------------------------------------------------------
steelhead_gg <- p +
  #ggplot() +
  #geom_sf(data = states, fill = 'white') +
  geom_sf(data = SR_basin, alpha = .5, inherit.aes = FALSE) +
  #geom_sf(data = SR_st_pop, aes(fill = TRT_POPID), inherit.aes = FALSE) +
  geom_sf(data = large_st, colour = 'blue', inherit.aes = FALSE) +
  geom_point(data = arrays, aes(x = lng, y = lat,
                            colour = BPA_Funding_OM), size = 3, inherit.aes = FALSE) +
  #ggrepel::geom_label_repel(data = arrays, aes(x = lng, y = lat, label = SiteID), size = 2) +
  #geom_label(data = centers, aes(x = x, y = y, label = TRT_POPID), size = 3, inherit.aes = FALSE) +
  scale_fill_viridis_d(alpha = .35, direction = -1, option = "D", guide = FALSE) +
  #scale_colour_manual(values = topo.colors(4)[1:3], labels = labs) +
  scale_colour_viridis_d(option = "C", end = .8, labels = labs) +
  #scale_colour_brewer(palette = 'PuRd', labels = labs) +
  theme_void() +
  theme(text = element_text(family = 'serif'),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(colour = 'grey35'),
        legend.position = c(.18,.12),
        legend.background = element_rect(fill = alpha('grey75', .5))) +
  labs(title = "Snake River Basin In-stream PIT-tag Detection Systems",
       subtitle = "Sixty-seven IPTDS are deployed across the basin to monitor abundance,\n life history and productivity of Steelhead and spring/summer Chinook salmon.",
       colour = 'Operation and Maintenance')# +
  #coord_sf(xlim = c(-112,-120), ylim = c(43,48))

steelhead_gg

ggsave(steelhead_gg, filename = './Figures/IPTDS_sites_map.png', height = 9, width = 9)

# Steelhead POPs and DABOM sites
source('./R/definePopulations.R')
pop_names <- definePopulations('Steelhead') %>%
  distinct(TRT) %>%
  mutate(est = TRUE)

# get center point for TRT_POPID labels
poly <- as_Spatial(SR_st_pop)
centers <- data.frame(rgeos::gCentroid(poly, byid = TRUE))
centers$TRT_POPID <- SR_st_pop$TRT_POPID


st_centers <- left_join(centers, pop_names, by = c('TRT_POPID' = 'TRT')) %>%
  mutate(est = ifelse(is.na(est), FALSE, est))

config <- bind_cols(config, 
                    st_coordinates(config) %>%
                      as_tibble() %>%
                      rename(lng = X, lat = Y))

labs <- config %>% group_by(SiteTypeName) %>% summarise(n = n()) %>% transmute(labs = paste0(SiteTypeName, ' (n = ',n,')')) %>% pull(labs)


steelhead_dabom <- p +
  geom_sf(data = SR_st_pop, aes(fill = TRT_POPID), inherit.aes = FALSE) +
  geom_sf(data = large_st, colour = 'blue', inherit.aes = FALSE) +
  geom_point(data = config, aes(x = lng, y = lat,
                                colour = SiteTypeName), size = 3, inherit.aes = FALSE) +
  geom_label(data = st_centers[st_centers$est,], aes(x = x, y = y, label = TRT_POPID), size = 3, inherit.aes = FALSE) +
  scale_fill_viridis_d(alpha = .5, direction = -1, option = "D", guide = FALSE) +
  scale_colour_brewer(palette = 'PuRd', labels = labs) +
  theme_void() +
  theme(text = element_text(family = 'serif'),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(colour = 'grey35'),
        legend.position = c(.22,.2),
        legend.background = element_rect(fill = alpha('grey75', .5))) +
  labs(title = "Snake Basin Steelhead Populations",
       subtitle = "Adult escapement is estimated for fish returning to tributaries located within twenty-two populations using the DABOM model and \n 114 distinct PIT-tag observation sites.",
       colour = 'Detection Type')

steelhead_dabom
ggsave(steelhead_dabom,filename = './Figures/steelhead_dabom_sites.png', height = 9, width = 9)

# Chinook POPs and DABOM sites

pop_names <- definePopulations('Chinook') %>%
  distinct(TRT) %>%
  mutate(est = TRUE)

# get center point for TRT_POPID labels
poly <- as_Spatial(SR_ch_pop)
centers <- data.frame(rgeos::gCentroid(poly, byid = TRUE))
centers$TRT_POPID <- SR_ch_pop$TRT_POPID

ch_centers <- left_join(centers, pop_names, by = c('TRT_POPID' = 'TRT')) %>%
  mutate(est = ifelse(is.na(est), FALSE, est))

config <- bind_cols(config, 
                    st_coordinates(config) %>%
                      as_tibble() %>%
                      rename(lng = X, lat = Y))

labs <- config %>% group_by(SiteTypeName) %>% summarise(n = n()) %>% transmute(labs = paste0(SiteTypeName, ' (n = ',n,')')) %>% pull(labs)


chinook_dabom <- p +
  geom_sf(data = SR_ch_pop, aes(fill = TRT_POPID), inherit.aes = FALSE) +
  geom_sf(data = large_st, colour = 'blue', inherit.aes = FALSE) +
  geom_point(data = config, aes(x = lng, y = lat,
                                colour = SiteTypeName), size = 3, inherit.aes = FALSE) +
  geom_label(data = ch_centers[ch_centers$est,], aes(x = x, y = y, label = TRT_POPID), size = 3, inherit.aes = FALSE) +
  scale_fill_viridis_d(alpha = .5, direction = -1, option = "D", guide = FALSE) +
  scale_colour_brewer(palette = 'PuRd', labels = labs) +
  theme_void() +
  theme(text = element_text(family = 'serif'),
        plot.title = element_text(face = 'bold'),
        plot.subtitle = element_text(colour = 'grey35'),
        legend.position = c(.22,.2),
        legend.background = element_rect(fill = alpha('grey75', .5))) +
  labs(title = "Snake Basin Spring/Summer Chinook Salmon Populations",
       subtitle = "Adult escapement is estimated for fish returning to tributaries located within thirty-one populations using the DABOM model and \n 114 distinct PIT-tag observation sites.",
       colour = 'Detection Type')

chinook_dabom
ggsave(chinook_dabom,filename = './Figures/chinook_dabom_sites.png', height = 9, width = 9)

# Plot dynamic maps

library(leaflet)
#library(mapview)
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
  tag.map.title, HTML("Snake Basin Steelhead Populations")
)  

st_copop <- colorFactor(palette = 'viridis', domain = SR_st_pop$POP_NAME)

array <- filter(config, DetectionType == 'In-stream Array')
weirs <- filter(config, DetectionType == 'Weir Scanning')
ladder <- filter(config, DetectionType == 'Ladder Observation')

tmp_config <- filter(config, !is.na(DetectionType))

pal <- colorFactor(
  palette = topo.colors(9)[1:3],
  #reverse = TRUE,
  domain = tmp_config$DetectionType
)

steelhead_sites <- leaflet() %>%
  addProviderTiles(providers$Esri.WorldTopoMap) %>%
  addPolygons(data = SR_st_pop,
              popup = ~as.character(POP_NAME),
              label = ~TRT_POPID,
              layerId = SR_st_pop$POP_NAME,
              group = 'Populations',
              stroke = TRUE, color = 'black',
              fillColor = ~st_copop(POP_NAME), weight = 1, fillOpacity = .5) %>%
  addPolylines(data = large_st, weight = 1) %>%
  addCircleMarkers(data = array, radius = 5,
                   group = 'In-Stream Arrays',
                   popup = paste("<b>Site ID:</b>", array$SiteID, "<br>",
                                 "<b>BPA Funding:</b>", array$BPA_Funding_OM, "<br>",
                                 "<b>O&M:</b>", array$Agency_OM, "<br>"),
                   color = ~pal(array$DetectionType),
                  #clusterOptions = markerClusterOptions(),
                  popupOptions = popupOptions(noHide = T)) %>%
  addCircleMarkers(data = weirs, radius = 5,
                   group = 'Weir Scanning',
                   popup = paste("<b>Site ID:</b>", weirs$SiteID, "<br>",
                                 "<b>BPA Funding:</b>", weirs$BPA_Funding_OM, "<br>",
                                 "<b>O&M:</b>", weirs$Agency_OM, "<br>"),
                   color = ~pal(weirs$DetectionType),
                   #clusterOptions = markerClusterOptions(),
                   popupOptions = popupOptions(noHide = T)) %>%
  addCircleMarkers(data = ladder, radius = 5,
                   group = 'Ladder Observations',
                   popup = paste("<b>Site ID:</b>", ladder$SiteID, "<br>",
                                 "<b>BPA Funding:</b>", ladder$BPA_Funding_OM, "<br>",
                                 "<b>O&M:</b>", ladder$Agency_OM, "<br>"),
                   color = ~pal(ladder$DetectionType),
                   #clusterOptions = markerClusterOptions(),
                   popupOptions = popupOptions(noHide = T)) %>%
  # addLabelOnlyMarkers(data = centers, lng = ~x, lat = ~y,
  #                     label = ~as.character(TRT_POPID),
  #                     labelOptions = labelOptions(noHide = T,
  #                                                  direction = 'top', textOnly = T)) %>%
  addControl(st_title, position = 'topleft', className = "map-title") %>%
  addLegend(data = tmp_config, position = "bottomleft", pal = pal, values = ~DetectionType,
          title = "Detection Type",
          #labFormat = labelFormat(prefix = "$"),
          opacity = 1) %>%
  addLayersControl(
    overlayGroups = c("Populations", "In-Stream Arrays", "Ladder Observations", "Weir Scanning"),
    options = layersControlOptions(collapsed = FALSE)
  )

steelhead_sites

saveWidget(steelhead_sites, file = 'steelhead_sites.html')
#mapshot(steelhead_sites, file = './Figures/steelhead_sites.png')


# %>%
#   addCircles(SR_array_sites$Longitude, SR_array_sites$Latitude, group = 'sites', radius = 100,
#              color = 'blue', 
#              popup = paste("<b>Site ID:</b>", SR_array_sites$SiteID, "<br>"),
#              popupOptions = popupOptions(noHide = T, textsize = "15px"),
#              )


# load function and assign group variables to all sites; TRT POPs and GSI
# source('./R/assign_POP_GSI.R')
# 
# # set species
# spp = 'Steelhead'
# # set year
# yr = 2018
# 
# # mapping data
# pop_ls <- assign_POP_GSI(species = spp, configuration, site_df)
# grp_df <- pop_ls[[1]]
# map_df <- pop_ls[[2]]

#site_labs <- distinct(grp_df, SiteID)

# Sex....
# left_join(map_df, modSexDf, by = c('TRT_POPID' = 'TRT')) %>%
#  ggplot() +
#    geom_sf(aes(fill = propF)) +
#    scale_fill_continuous(type = 'viridis') +
# #   geom_sf(data = grp_df) +
#    geom_sf_text(aes(label = round(propF,2)), size = 3)


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
