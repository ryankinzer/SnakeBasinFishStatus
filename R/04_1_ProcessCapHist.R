# Author: Ryan N. Kinzer
# Purpose: Process LGR tag histories for DABOM using PITcleanr v2.0
# Created: 6/28/2021
# Last Modified: 
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(PITcleanr)
source('./R/identifyFishType.R')

#-----------------------------------------------------------------
# set up folder structure
PITcleanrFolder = 'data/PITcleanr_v2'
if(!dir.exists(PITcleanrFolder)) {
  dir.create(PITcleanrFolder)
}

# load configuration and site_df data
load('./data/ConfigurationFiles/site_config.rda')

#-----------------------------------------------------------------
# read observations from PTAGIS and process with PITcleanr
# set species
spp <- 'Steelhead'
yr <- 2020

## Compress Observations----
ptagis_file <- paste0('./data/CompleteTagHistories/LGR_',spp,'_',yr,'.csv')

comp_obs <- compress(ptagis_file, configuration = configuration, # would like to keep site_code
                     ignore_event_vs_release = TRUE,
                     units = 'days')

# determine trap date, and remove detections prior to that

if(spp == 'Chinook'){
  start_date = lubridate::ymd(paste0(yr,'0301'))
  end_date = lubridate::ymd(paste0(yr,'0817'))
} else {
  start_date = lubridate::ymd(paste0(yr-1,'0701'))
  end_date = lubridate::ymd(paste0(yr,'0630'))
}

obs_direct = comp_obs %>%
  # get the first detection of each tag at Lower Granite Dam
  left_join(comp_obs %>%
              filter(node == "GRA",
                     event_type_name %in% c("Mark", "Recapture")) %>%  #trap observations
              filter(min_det >= start_date) %>% # eliminates observations prior to spawn year
              group_by(tag_code) %>%
              filter(min_det == min(min_det)) %>%
              summarise(start_date = min_det,
                        .groups = "drop"),
            by = "tag_code") %>%
  # filter any detections before the "start_date"
  filter(min_det >= start_date) %>%
  filter(!(node == 'GRA' & event_type_name == 'Observation')) %>% # removes additional detection at GRA upstream of trap
  group_by(tag_code) %>%
  # re-calculate the "slots" for each tag_code
  mutate(slot = slot - min(slot) + 1) %>%
  # add direction using "addDirection()
  addDirection(parent_child = parent_child_node) %>%
  ungroup()

# test out obs type
tmp <- obs_direct %>%
  group_by(tag_code) %>%
  nest() %>%
  mutate(
    #passage_type = map(data,
    #                .f = ~passageType(., fallback_days)),
    obs_type = map(data,
                   .f = ~obsType(., spawn_year = yr)),
  ) %>%
  select(-data) %>%
  unnest(obs_type) %>%
  ungroup() %>%
  mutate(index = 1:n()) %>%
  select(index, tag_code, node, min_det, obs_type, everything()) %>%
  arrange(tag_code, min_det)

write_csv(tmp, file = './data/PITcleanr_v2/test_obs.csv')

# previous methods.
all_obs <- identifyFishType(obs_direct, fallback_days = , pop_days = , lgd_days = 5) # maybe this should return a separate dataset and only run for steelhead?

all_obs <- all_obs %>%
  select(-node_order, -path, -direction)

# will not work if addDirection is already included in compress_obs - removed
# in line above

prepped_all_obs = filterDetections(compress_obs = all_obs,
                               parent_child = parent_child_node,
                               max_obs_date = NULL)

no_kelt_obs <- all_obs %>%
  filter(!kelt_obs)
  
prepped_obs = filterDetections(compress_obs = no_kelt_obs,
                              parent_child = parent_child_node,
                              max_obs_date = NULL)

list('all_obs' = prepped_all_obs,
     'prepped_obs' = prepped_obs) %>%
writexl::write_xlsx(
  path = paste0(PITcleanrFolder,'/PreppedObs_',spp,'_',yr,'.xlsx')
)
