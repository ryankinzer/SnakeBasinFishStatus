# Purpose: Download PIT-tag observations at Snake Basin sites from DART.
# Author: Ryan N. Kinzer
# Created: 08/03/20
# Last Modified:
# Notes: 

# Load Packages----
library(tidyverse)
library(lubridate)
library(PITcleanr)

# Script Variables----
spp <- 'Steelhead'
spawn_yr <- 2017
#yr_range <- as.list(2019)
#names(yr_range) <- 2019

# build configuration table (requires internet connection)
org_config = buildConfig()

# customize some nodes based on DABOM framework
configuration = org_config %>%
  mutate(Node = ifelse(SiteID %in% c('VC2', 'VC1', 'LTR', 'MTR', 'UTR'),
                       SiteID,
                       Node),
         Node = ifelse(SiteID == 'SC2',
                       'SC2B0',
                       Node),
         Node = ifelse(SiteID %in% c('CROTRP',
                                     'CRT',
                                     'REDTRP',
                                     'REDR',
                                     'RRT'),
                       'SC2A0',
                       Node),
         Node = ifelse(Node == 'ACB',
                       'ACBB0',
                       Node),
         Node = ifelse(Node == 'CCA',
                       'CCAB0',
                       Node),
         Node = ifelse(SiteID == 'AFC',
                       ifelse(grepl('MAINSTEM', AntennaGroup),
                              'AFCB0',
                              'AFCA0'),
                       Node),
         Node = ifelse(SiteID == 'HBC',
                       'HYCA0',
                       Node),
         Node = ifelse(SiteID %in% c('TUCH', 'TFH'),
                       'TUCH',
                       Node),
         Node = ifelse(SiteID == 'MCCA',
                       'STR',
                       Node),
         Node = ifelse(SiteID == 'CARMEC',
                       'CRCA0',
                       Node),
         Node = ifelse(SiteID == 'BIG2C',
                       'TAYA0',
                       Node),
         Node = ifelse(SiteID == 'WIMPYC',
                       'WPCA0',
                       Node),
         Node = ifelse(SiteID == 'IML' & ConfigID == 130 & AntennaID == '09',
                       'IMLA0',
                       Node),
         Node = str_replace(Node, '^BTC', 'BTL'),
         Node = ifelse(SiteID %in% c('YANKFK', 'CEY'),
                       'YFKA0',
                       Node),
         Node = ifelse(SiteID == 'SAWT',
                       'STL',
                       Node),
         Node = ifelse(SiteID == 'LOOH',
                       'LOOKGC',
                       Node),
         Node = ifelse(SiteID == 'RPDTRP',
                       'RAPH',
                       Node),
         Node = ifelse(SiteID == 'CHARLC',
                       'CCAB0',
                       Node),
         Node = ifelse(Node == 'KEN',
                       'KENB0',
                       Node),
         Node = ifelse(Node == 'HYC',
                       'HYCB0',
                       Node),
         Node = ifelse(Node == 'YFK',
                       'YFKB0',
                       Node),
         Node = ifelse(Node == 'LLR',
                       'LLRB0',
                       Node),
         Node = ifelse(Node == 'LRW',
                       'LRWB0',
                       Node),
         Node = ifelse(SiteID == '18M',
                       str_replace(Node, '18M', 'HEC'),
                       Node)) %>%
  distinct()


# Node network for DABOM
site_df = writeLGRNodeNetwork()

# remove some sites that have been combined with others (see the modifications to the configuration file)
site_df = site_df %>%
  filter(!SiteID %in% c('TFH',
                        'MCCA',
                        'WIMPYC',
                        'YANKFK', 'CEY',
                        'SAWT',
                        'LOOH',
                        'CARMEC',
                        'BIG2C',
                        'RPDTRP'))

dart_obs <- processDART_LGR(species = spp, spawnYear = spawn_yr,
                            configuration = configuration,
                            truncate = TRUE)

proc_ch <- dart_obs[[1]]

mark_data <- dart_obs[[2]] %>%
  mutate(spawn_yr = spawn_yr,
         rel_year = year(rel_date),
         brood_year = rel_year -1) %>%
  select(spawn_yr, TagID = tag_id, file_id, contains('mark_'), contains('rel_'), contains('t_'), rrkm, length, contains('trans_')) %>%
  distinct()

gra_H_obs <- mark_data %>%
  filter(t_rear_type == 'H') %>%
  filter(!grepl('LGR', mark_site)) %>%
  filter(!grepl('TRP', mark_site)) %>%
  group_by(spawn_yr, t_species, t_run, t_rear_type, mark_site, rel_site, rrkm, rel_year, brood_year) %>%
  tally()

write_csv(gra_H_obs, path = paste0('./data/PBTData/spawn_yr_',spawn_yr,'_H_gra_obs.csv'))

iptds_H_obs <- proc_ch %>%
  left_join(mark_data) %>%
  filter(t_rear_type == 'H') %>%
  filter(!grepl('LGR', mark_site)) %>%
  filter(!grepl('TRP', mark_site)) %>%
  filter(Group != 'NA') %>%
  group_by(spawn_yr, t_species, t_run, t_rear_type, Group, mark_site, rel_site, rrkm, rel_year, brood_year) %>%
  summarise(obs = n_distinct(TagID))

tmp <- full_join(gra_H_obs, iptds_H_obs)

write_csv(tmp, path = paste0('./data/PBTData/spawn_yr_',spawn_yr,'_H_obs.csv'))
