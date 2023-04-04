# Author: Ryan N. Kinzer
# Purpose: Download and configure IPTDS infrastructure for tag observation
# processing and the DABOM model.
# Created: 9/27/19
# Last Modified: 4/22/22
# Notes: 
# The output and saved file from this script is used for processing tag
# observations, with TRT population and GSI grouping designations
# and for visualizing infrastructure and mapping tag detections.
# 
# Updates were made in SY2021 runs to include sites downstream of LGR
# and to change node names to include _D, _M and _U.
#-----------------------------------------------------------------

# Install PITcleanr

devtools::install_github('ryankinzer/PITcleanr', ref = 'snake_main')

# load needed libraries
library(tidyverse)
library(PITcleanr)

library(ggraph)
library(tidygraph)
source('./R/buildNetwork_tbl.R')

# download meta data for all ptagis INT and MRR sites
ptagis_sites <- buildConfig()

# customize some nodes because of name changes across the years and combining of sites into single nodes
configuration = ptagis_sites %>%
  mutate(
    #site_code = ifelse(site_code == 'LGRLDR', 'GRA', site_code),
    node = case_when(
      site_code %in% c('GRA', 'LGRLDR', 'LGR') ~ 'GRA',
      site_code %in% c('GRJ', 'GRX', 'GRS') ~ 'GRS',   
      site_code %in% c('LGS', 'GOJ', 'GOA') ~ 'GOA',
      site_code %in% c('LMN', 'LMJ', 'LMA') ~ 'LMA',
      site_code %in% c('IHR', 'ICH' ,'IHA') ~ 'IHR',
      site_code %in% c('MCN', 'MCJ', 'MCX', 'MC1', 'MC2') ~ 'MCN',
      site_code %in% c('JDJ', 'JDA', 'JO1', 'JO2') ~ 'JDA',
      site_code %in% c('TDA', 'TD1', 'TD2') ~ 'TDA',
      site_code %in% c('B2A', 'BONAFF', 'BO1', 'BO2', 'BON', 'BVJ', 'B1J', 'BVX', 'B2J', 'BO4', 'BWL', 'BO3', 'BCC') ~ 'BON',
      as.numeric(str_extract(rkm, '\\d+')) > 539 ~ 'UpperColumbia', #Upper Columbia      
      as.numeric(str_extract(rkm, '\\d+')) == 539 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'Yakima', #Yakima
      as.numeric(str_extract(rkm, '\\d+')) == 509 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'WallaWalla', #Walla Walla
      as.numeric(str_extract(rkm, '\\d+')) == 465 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'Umatilla', #Umatilla
      as.numeric(str_extract(rkm, '\\d+')) == 351 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'JohnDay', #JohnDay
      as.numeric(str_extract(rkm, '\\d+')) == 328 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'Deschutes', #Deschutes
      as.numeric(str_extract(rkm, '\\d+')) == 273 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'HoodRiver', #HoodRiver
      as.numeric(str_extract(rkm, '\\d+')) == 251 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'WindRiver', #WindRiver
      as.numeric(str_extract(rkm, '\\d+')) == 261 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'LittleWhite', #LittleWhite
      as.numeric(str_extract(rkm, '\\d+')) == 271 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'WhiteSalmon', #WhiteSalmon
      as.numeric(str_extract(rkm, '\\d+')) == 290 & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 0 ~ 'Klickitat', #Klickitat
      TRUE ~ node),
    
    node = ifelse(site_code %in% c('DWL', 'DWOR'),
                  'DWL',
                  node),
    node = ifelse(node == 'CROTRP',
                  'CRT',
                  node),
         node = ifelse(node %in% c('REDTRP','REDR'),
                       'RRT',
                       node),
         node = ifelse(site_code == 'AFC',
                       ifelse(grepl('MAINSTEM', antenna_group),
                              'AFC_D', #mainstem AFC becomes _D
                              'AFC_U'),# south fork and north fork become _U
                       node),
         node = ifelse(node == 'TUCH',
                       'TFH_U',
                       node),
         node = ifelse(site_code == 'MCCA',
                       'SALSFW',
                       node),
         node = ifelse(site_code == 'CARMEC',
                       'CRC_U',
                       node),
         node = ifelse(site_code == 'BIG2C',
                       'TAY_U',
                       node),
         node = ifelse(site_code == 'WIMPYC',
                       'WPC_U',
                       node),
         node = ifelse(site_code == 'IML' & config_id == 130 & antenna_id == '09',
                       'IML_U',
                       node),
         node = str_replace(node, '^BTC', 'BTL'),
         node = ifelse(site_code %in% c('YANKFK', 'CEY'),
                       'YFK_U',
                       node),
         # node = ifelse(site_code == 'SAWT', # included in parent-child
         #               'STL',
         #               node),
         node = ifelse(site_code == 'LOOKGC',
                       'LOOH',
                       node),
         node = ifelse(site_code == 'RPDTRP',
                       'RAPH',
                       node),
         node = ifelse(site_code == 'CHARLC',
                       'CCA_U',
                       node),
         node = ifelse(site_code == '18M',
                       str_replace(node, '18M', 'HEC'),
                       node),
        node = ifelse(site_code == 'BEARVC',
                      'BRC',
                      node),
        node = ifelse(site_code == 'POTREF',
                      'EFPW',
                      node)
         )

# Append Population Names
load('../../DFRM Projects/River_Mapping/data/polygons/SR_pops.rda')

SR_st_pop <- sf::st_as_sf(sth_pop) %>%
  select(st_ESU = ESU_DPS, st_MPG = MPG, st_POP_NAME = POP_NAME, st_TRT_POPID = TRT_POPID, st_GSI_Group = GSI_Group)

SR_ch_pop <- sf::st_as_sf(spsm_pop) %>%
  select(ch_ESU = ESU_DPS, ch_MPG = MPG, ch_POP_NAME = POP_NAME, ch_TRT_POPID = TRT_POPID, ch_GSI_Group = GSI_Group)

configuration <- configuration %>%
  filter(!is.na(latitude)) %>%
  sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)

configuration <- sf::st_join(configuration,
                             SR_st_pop)

configuration <- sf::st_join(configuration,
                             SR_ch_pop)

configuration <- configuration %>%
  mutate(across(c(contains('st_'), contains('ch_')), ~ifelse(grepl('^GRA$|^GRS$|^GOA$|^NPTH$|^DWL$', node), NA, .))) %>%
  mutate(across(c(ch_POP_NAME, ch_TRT_POPID), ~ifelse(grepl('^SW1$|^SW2$', node), 'SEUMA/SEMEA/SEMOO', .))) %>%
  mutate(st_POP_NAME = ifelse(grepl('^SC1$|^SC2$', node), 'South Fork Clearwater River', st_POP_NAME)) %>%
  mutate(ch_POP_NAME = ifelse(grepl('^SC1$|^SC2$', node), 'Upper South Fork Clearwater', ch_POP_NAME)) %>%
  mutate(st_TRT_POPID = ifelse(grepl('^SC1$|^SC2$', node), 'CRSFC-s',st_TRT_POPID)) %>%
  mutate(ch_TRT_POPID = ifelse(grepl('^SC1$|^SC2$', node), 'SCUMA', ch_TRT_POPID))


#View(configuration %>% sf::st_set_geometry(NULL))

# Create parent/child table----
config <- configuration %>% sf::st_set_geometry(NULL) # for use with parent child table
root_site <- 'BON' #'GRA' #'noGRS' #'test' #GRA' #'BON'
parent_child <- read_csv(paste0('./data/ConfigurationFiles/parent_child_',root_site,'.csv'))

parent_child <- parent_child[parent_child$child != 'CLWH',] # remove Clearwater hatchery for now...causing problems with SC1 and SC2
#parent_child <- parent_child[parent_child$child != 'POTRWF',] # 


# add rkm
parent_child <- parent_child %>%
  left_join(config %>%
              select(parent = site_code,
                     parent_rkm = rkm_total) %>%
              filter(!is.na(parent_rkm)) %>%
              distinct()
  ) %>%
  left_join(config %>%
              select(child = site_code,
                     child_rkm = rkm_total) %>%
              filter(!is.na(child_rkm)) %>%
              distinct()
  ) %>%
  select(parent, child, parent_rkm, child_rkm)


parent_child <- parent_child %>%
  mutate(child_rkm = case_when(
    child == 'UpperColumbia' ~ 540,
    child == 'Yakima' ~ 539,
    child == 'WallaWalla' ~ 509,
    child == 'Umatilla' ~ 465,
    child == 'JohnDay' ~ 351,
    child == 'Deschutes' ~ 328,
    child == 'HoodRiver' ~ 273,
    child == 'WindRiver' ~ 251,
    child == 'LittleWhite' ~ 261,
    child == 'WhiteSalmon' ~ 271,
    child == 'Klickitat' ~ 290,
    child == 'OXBO' ~ 850,
    TRUE ~ child_rkm
  ))

#plotNodes(parent_child)

parent_child <- left_join(parent_child,
                          config %>%
                            #sf::st_set_geometry(NULL) %>%
                            select(site_code, MPG = st_MPG, POP_NAME = st_POP_NAME, TRT_POPID = st_TRT_POPID) %>%
                            distinct(),
                          by = c('child' = 'site_code')) %>%
  arrange(MPG, POP_NAME, parent, child)

node_attributes <- tibble(label = union(parent_child$child, parent_child$parent)) %>%
  left_join(config %>%
             #sf::st_set_geometry(NULL) %>%
             select(label = site_code, MPG = st_MPG, POP_NAME = st_POP_NAME,
                    TRT_POPID = st_TRT_POPID, site_type, site_type_name, rkm_total) %>%
            distinct(),
            by = 'label') %>%
  mutate(detection_type = case_when(
    rkm_total > 695 & label != 'TPJ' ~ 'Spawner/Kelt/Repeat Spawner',
    label == 'OXBO' ~ 'Spawner/Kelt/Repeat Spawner',
    label == 'GRA' ~ 'Release',
    TRUE ~ 'Kelt/Repeat Spawner'
    ),
    group = case_when(
      grepl('Spawner/Kelt/Repeat Spawner', detection_type) ~ MPG,
      TRUE ~ 'Below LWG'
    ))

node_graph <- buildNetwork_tbl(parent_child, node_attributes)

plasma_pal <- c(viridis::plasma(n = 6, begin = .5), 'grey90')

site_network <- ggraph(node_graph, layout = 'tree') +
  geom_edge_bend() +
  #geom_node_point(aes(shape = detection_type, colour = group), size = 2) +
  geom_node_label(aes(label = label, fill = group), size = 1) +
  scale_fill_manual(values = plasma_pal,
                    breaks = c('Clearwater River', 'Hells Canyon', 'Grande Ronde River',
                               'Imnaha River', 'Salmon River', 'Lower Snake')) +
  guides(
    fill = guide_legend(
      title = '', #Spawner Locations',
      override.aes = aes(label = ''),
      nrow = 1
    )
  ) +
  theme_void() +
  theme(legend.position = 'bottom')

site_network

ggsave(paste0('./Figures/site_network_',root_site,'.png'), site_network, width = 14, height = 8.5)


# build network graphs for nodes
pc_nodes <- addParentChildNodes(parent_child, config) # only works with rk_edits to PITcleanr

node_attributes <- tibble(label = union(pc_nodes$child, pc_nodes$parent)) %>%
  left_join(config %>%
             # sf::st_set_geometry(NULL) %>%
              select(label = node, MPG = st_MPG, POP_NAME = st_POP_NAME,
                     TRT_POPID = st_TRT_POPID,
                     rkm_total) %>%
              distinct(),
            by = 'label') %>%
  group_by(label) %>%
  slice(which.max(rkm_total)) %>%
  mutate(detection_type = case_when(
    rkm_total > 695 & !grepl('TPJ', label) ~ 'Spawner/Kelt/Repeat Spawner',
    label == 'OXBO' ~ 'Spawner/Kelt/Repeat Spawner',
    label == 'GRA' ~ 'Release',
    TRUE ~ 'Kelt/Repeat Spawner'
  ),
  group = case_when(
    grepl('Spawner/Kelt/Repeat Spawner', detection_type) ~ MPG,
    TRUE ~ 'Below LWG'
  )) %>%
  select(label, group) %>%
  distinct()

node_graph <- buildNetwork_tbl(pc_nodes, node_attributes)

plasma_pal <- c(viridis::plasma(n = 6, begin = .5), 'grey90')

node_network <- ggraph(node_graph, layout = 'tree') +
  geom_edge_bend() +
  geom_node_label(aes(label = label, fill = group), size = 1) +
  scale_fill_manual(values = plasma_pal,
                    breaks = c('Clearwater River', 'Hells Canyon', 'Grande Ronde River',
                               'Imnaha River', 'Salmon River', 'Lower Snake')) +
  guides(
    fill = guide_legend(
      title = '', #Spawner Locations',
      override.aes = aes(label = ''),
      nrow = 1
    )
  ) +
  theme_void() +
  theme(legend.position = 'bottom')

node_network

ggsave(paste0('./Figures/node_network_',root_site,'.png'), node_network, width = 14, height = 8.5)

# Save file.
save(configuration, parent_child, pc_nodes, file = paste0('./data/ConfigurationFiles/site_config_',root_site,'.rda'))

# Map sites and nodes.
library(sf)
library(maps)
library(ggmap)
library(ggspatial)
library(cowplot)

source('./R/theme_map.R')
load('../../DFRM Projects/River_Mapping/data/polygons/SR_pops.rda')
load('../../DFRM Projects/River_Mapping/data/flowlines/large_rivers.rda')

states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
pnw <- states %>% filter(ID %in% c('idaho', 'oregon', 'washington')) %>%
  st_transform(crs = 4326)

sth_pop <- sth_pop %>%
  filter(TRT_POPID != 'CRNFC-s')# %>%
  #filter(TRT_POPID != 'SNTUC-s')

bb <- st_bbox(sth_pop)

pnw_map <- ggplot() +
  geom_sf(data = pnw, inherit.aes = FALSE) +
  geom_sf(data = sth_pop, fill = 'grey30', inherit.aes = TRUE) +
  #geom_sf(data = pnw_rivers, colour = 'cyan', inherit.aes = FALSE) +
  theme_map()


map_nodes <- tibble(site_code = union(parent_child$parent, parent_child$child)) %>%
  inner_join(configuration %>% 
               select(starts_with('site_')) %>%
               distinct())

snake_map <- ggplot() +
  geom_sf(data = pnw, fill = NA, inherit.aes = FALSE) +
  geom_sf(data = sth_pop, aes(fill = MPG), inherit.aes = TRUE) +
  geom_sf(data = pnw_rivers, colour = 'cyan', inherit.aes = FALSE) +
  geom_sf(data = snake_rivers, colour = 'cyan', inherit.aes = FALSE) +
  geom_sf(data = map_nodes, aes(geometry = geometry), size = 2) + 
  # ggrepel::geom_text_repel(data = map_nodes, aes(geometry = geometry, label = site_code),
  #                          size = 2,
  #                          stat = 'sf_coordinates') +
  scale_fill_manual(values = plasma_pal,
                    breaks = c('Clearwater River', 'Hells Canyon', 'Grande Ronde River',
                               'Imnaha River', 'Salmon River', 'Lower Snake')) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  guides(fill = guide_legend(
    title = '',
    nrow = 1)) +
  coord_sf(xlim = c(bb$xmin, bb$xmax), ylim = c(bb$ymin, bb$ymax)) +
  theme_map() +
  theme(legend.position = 'bottom')

snake_map

full_map <- ggdraw() +
  draw_plot(snake_map) +
  draw_plot(pnw_map, x = .75, y = .8, width = .2, height = .2)

full_map

ggsave('./figures/site_map.png', full_map, width = 8, height = 7)
