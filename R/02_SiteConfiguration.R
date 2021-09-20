# Author: Ryan N. Kinzer
# Purpose: Download and configure IPTDS infrastructure for tag observation
# processing and the DABOM model.
# Created: 9/27/19
# Last Modified: 
# Notes: 
# The output and saved file from this script is used for processing tag
# observations, with TRT population and GSI grouping designations
# and for visualizing infrastructure and mapping tag detections.
#-----------------------------------------------------------------

# load needed libraries
library(tidyverse)
#library(sf)
#library(magrittr)
#library(readxl)
#library(WriteXLS)
library(PITcleanr)


# download meta data for all ptagis INT and MRR sites
ptagis_sites = queryPtagisMeta()

# clean things up a bit
clean_sites = ptagis_sites %>%
  mutate(node = NA) %>%
  rename(site_type_name = site_type,
         site_type = type,
         config_id = configuration_sequence,
         antenna_group = antenna_group_name) %>%
  select(site_code,
         config_id,
         antenna_id,
         node,
         start_date,
         end_date,
         site_type,
         site_name,
         antenna_group,
         site_description,
         site_type_name,
         rkm,
         rkm_total,
         latitude,
         longitude) 

node_config <- clean_sites %>%
  mutate(node = ifelse(grepl('DOWNSTREAM', antenna_group, ignore.case = T) |
                         grepl('DNSTREAM', antenna_group, ignore.case = T) |
                         grepl('LOWER', antenna_group, ignore.case = T) |
                         grepl('BOTTOM', antenna_group, ignore.case = T),
                       paste0(site_code, '_D'),
                       node),
         node = ifelse(grepl('MIDDLE', antenna_group, ignore.case = T) |
                         grepl('MIDDLE', antenna_group, ignore.case = T),
                       paste0(site_code, '_M'),
                       node),
         node = ifelse(grepl('UPSTREAM', antenna_group, ignore.case = T) |
                         grepl('UPPER', antenna_group, ignore.case = T) |
                         grepl('TOP', antenna_group, ignore.case = T),
                       paste0(site_code, '_U'),
                       node) ,
         # node = ifelse(grepl('^LGR', site_code),
         #               'GRA',
         #               node),
         node = ifelse(is.na(node), site_code, node))

# # for any site that has some nodes with "D", "U", but some configurations with a single node, make that node "D"
node_config = node_config %>%
  group_by(site_code) %>%
  mutate(node_site = sum(node == site_code), # within a site_code check to see if any configurations have node equal to site_code; i.e., a node without an assigned _D, _U or _M
         node_site_D = sum(node == paste0(site_code, "_D"))) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(node = if_else(node_site > 0 & node_site_D > 0 & !(grepl("_D$", node) | grepl("_M$", node) | grepl("_U$", node)),
                        paste0(site_code, '_D'),
                        node)) %>%
  ungroup() %>%
  select(-node_site,
         -node_site_D)

# need to create meta data to batch change sites to downriver locations of
# interest based in RKM
# Ice harbor: 522.016, LoMo: 522.067, L.Goose: 522.113, LGR: 522.173

site_meta <- node_config %>%
  select(contains('site'), contains('rkm'), latitude, longitude) %>%
  distinct() %>%
  mutate(detection_group = case_when(
    as.numeric(str_extract(rkm, '\\d+')) == 234 ~ 'BON',
    as.numeric(str_extract(rkm, '\\d+')) < 234 ~ 'Lower_C',
    as.numeric(str_extract(rkm, '\\d+')) < 522 ~ 'Mid_C',
    as.numeric(str_extract(rkm, '\\d+')) > 522 ~ 'Upper_C',
    as.numeric(str_extract(rkm, '\\d+')) == 522 ~ 'Snake',
    TRUE ~ NA_character_)) %>%
  mutate(detection_group = case_when(
    detection_group == 'Snake' & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) > 173 ~ 'Above_LGD',    
    detection_group == 'Snake' & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) == 113 ~ 'LGO',
    detection_group == 'Snake' & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) == 67 ~ 'LMO',    
    detection_group == 'Snake' & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) == 16 ~ 'IHR',
  TRUE ~ detection_group))

# Create configuration file by joining metadata objects
configuration <- node_config %>%
  left_join(site_meta %>%
              select(site_code, site_type, detection_group)) %>%
  mutate(node = ifelse(detection_group == 'Above_LGD' | detection_group == 'Snake' | detection_group == 'GRS', node, detection_group))

# customize some nodes because of name changes across the years and combining of sites into single nodes
configuration = configuration %>%
  mutate(
    node = case_when(
      site_code %in% c('GRA', 'LGRLDR') ~ 'GRA',
      as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) == 173 ~ 'GRS', 
      TRUE ~ node),
    node = ifelse(site_code %in% c('CLWH', 'DWL', 'DWOR', 'NPTH', 'CLWR', 'CLWTRP'),
                  'Clearwater',
                  node),
    node = ifelse(site_code %in% c('CAPEHC', 'SULFUC', 'ELKC', 'CAMASC', 'LOONC', 'SALMF2'),
                  'MFSalmon',
                  node),
    node = ifelse(site_code %in% c('SLATEC', 'CHAMWF', 'SALR1', 'TOWERC'),
                  'Salmon',
                  node),
    node = ifelse(site_code %in% c('GRNTRP', 'OXBO', 'SNAKE4', 'SNAKE3'),
                  'Snake',
                  node),
    node = ifelse(node == 'CROTRP',
                  'CRT',
                  node),
         node = ifelse(node %in% c('REDTRP',
                                     'REDR'),
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
                       'STR',
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
         node = ifelse(site_code == 'SAWT',
                       'STL',
                       node),
         node = ifelse(site_code == 'LOOKGC',
                       'LOOH',
                       node),
         node = ifelse(site_code == 'RPDTRP',
                       'RAPH',
                       node),
         node = ifelse(site_code == 'CHARLC',
                       'CCA_U',
                       node),
         node = ifelse(site_code == 'HEC',
                       str_replace(node, 'HEC', '18M'),
                       node),
        node = ifelse(site_code == 'BEARVC',
                      'BRC',
                      node)
         ) %>%
  distinct()

tmp <- configuration %>% select(detection_group, site_code, rkm, node) %>% distinct()# summarise(n_Nodes = n_distinct(node)) %>% arrange(detection_group)

# Create parent/child table----
parent_child <- read_csv('./data/ConfigurationFiles/parent_child.csv')

plotNodes(parent_child,
          label_size = 2,
          layout = "tree")

#ggsave('./Figures/site_network.png', tmp, width = 14, height = 8.5)

tmp_pc <- tibble(site_code = union(parent_child$child, parent_child$parent))

tmp_config <- configuration %>%
  select(site_code, node) %>%
  distinct() %>%
  arrange(site_code, node) %>%
  group_by(site_code) %>%
  mutate(n_nodes = n_distinct(node)) %>%
  mutate(node_num = paste("node", 1:n(), sep = "_")) %>%
  ungroup()

node_long <- left_join(tmp_pc, tmp_config) %>%
  distinct() %>%
  group_by(site_code) %>%
  mutate(n_nodes = n_distinct(node)) %>%
  mutate(node_num = paste("node", 1:n(), sep = "_"))
  
if(sum(is.na(node_long$node)) > 0) {
  node_long %>%
    filter(is.na(node)) %>%
    mutate(message = paste(site_code, "has a node that is NA.\n")) %>%
    pull(message) %>%
    warning()
}

node_wide = node_long %>%
  tidyr::pivot_wider(names_from = "node_num",
                     values_from = "node")

pc_nodes = parent_child %>%
  left_join(node_wide %>%
              rename(n_parent_nodes = n_nodes),
            by = c("parent" = "site_code")) %>%
  left_join(node_wide %>%
              rename(n_child_nodes = n_nodes),
            by = c("child" = "site_code")) %>%
  group_by(parent, child) %>%
  tidyr::nest(node_info = -any_of(names(parent_child))) %>%
  ungroup() %>%
  mutate(pc = map(node_info,
                  .f = function(x) {
                    
                    if(x$n_parent_nodes == 1) {
                      pc_new = x %>%
                        select(parent = node_1.x,
                               child = node_1.y)
                    }

                    if(x$n_parent_nodes == 2) {
                      pc_new = x %>%
                        select(parent = node_1.x,
                               child = node_2.x) %>%
                        bind_rows(x %>%
                                    select(parent = node_2.x,
                                           child = node_1.y))
                    }

                    if(x$n_parent_nodes == 3) {
                      pc_new = x %>%
                        select(parent = node_1.x,
                               child = node_2.x) %>%
                        bind_rows(x %>%
                                    select(parent = node_2.x,
                                           child = node_3.x)) %>%
                        bind_rows(x %>%
                                    select(parent = node_3.x,
                                           child = node_1.y))
                    }

                    if(x$n_child_nodes == 2) {
                      pc_new = pc_new %>%
                        bind_rows(x %>%
                                    select(parent = node_1.y,
                                           child = node_2.y))
                    }

                    if(x$n_child_nodes == 3) {
                      pc_new = pc_new %>%
                        bind_rows(x %>%
                                    select(parent = node_1.y,
                                           child = node_2.y)) %>%
                        bind_rows(x %>%
                                    select(parent = node_2.y,
                                           child = node_3.y))
                    }
                    
                    return(pc_new)
                  })) %>%
  select(-parent, -child) %>%
  tidyr::unnest(cols = pc) %>%
  select(parent, child) %>%
  distinct()

# create down river nodes----
new_nodes <- tribble(~"parent", ~"child",
                    # 'GRA', 'Clearwater',
                     'GRA', 'Snake',
                     'GRA', 'Salmon',
                     'GRA', 'MFSalmon',
                     'GRS', 'LGO',
                     'LGO', 'LMO',
                     'LGO', 'LYFE',
                     'LMO', 'IHR',
                     'IHR', 'Upper_C',
                     'IHR', 'Mid_C',
                     'Mid_C', 'BON',
                     'BON', 'Lower_C',
                     'LGO', 'LTR_D')

parent_child_node <- bind_rows(new_nodes, pc_nodes)

# need to switch order of GRS nodes.
parent_child_node$child[parent_child_node$parent == 'GRA' & parent_child_node$child == 'GRS_D'] <- 'GRS_U'
parent_child_node$parent[parent_child_node$parent == 'GRS_U'] <- 'GRS_D'
parent_child_node$parent[parent_child_node$child == 'GRS_M'] <- 'GRS_U'
parent_child_node$child[parent_child_node$parent == 'GRS_M'] <- 'GRS_D'

plotNodes(parent_child_node,
          layout = "tree")

# save plot
pc_node_plot <- plotNodes(parent_child_node,
                          layout = "tree",
                          point_size = 1,
                          label_size = 1)

ggsave('./Figures/node_network.png', pc_node_plot, width = 14, height = 8.5)


# Append Population Names
load('../../DFRM Projects/River_Mapping/data/polygons/Snake_POP_metadata.rda')

tmp <- configuration %>% 
  filter(!is.na(latitude)) %>%
  sf::st_as_sf(coords = c('longitude', 'latitude'), crs = 4326)

tmp <- sf::st_join(tmp,
                   SR_st_pop %>%
                     select(MPG, POP_NAME))

configuration <- tmp %>% st_set_geometry(NULL)

# Save file.
save(configuration, parent_child_node, file = './data/ConfigurationFiles/site_config.rda')
