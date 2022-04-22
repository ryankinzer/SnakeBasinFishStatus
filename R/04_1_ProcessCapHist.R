# Author: Ryan N. Kinzer
# Purpose: Process LGR tag histories for DABOM using PITcleanr v2.0
# Created: 6/28/2021
# Last Modified: 
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(PITcleanr)
source('../SR_Steelhead/R/identifyFishType.R')

#-----------------------------------------------------------------
# set up folder structure
PITcleanrFolder = 'data/PITcleanr_v2'
if(!dir.exists(PITcleanrFolder)) {
  dir.create(PITcleanrFolder)
}

#-----------------------------------------------------------------
# read observations from PTAGIS and process with PITcleanr
# set species
spp <- 'Chinook'
yr <- 2021
root_site <- 'BON'  #'noGRS'#'test'# 'BON'# or 'GRA'

# load configuration and site_df data
load(paste0('./data/ConfigurationFiles/site_config_',root_site,'.rda'))
#load('./data/ConfigurationFiles/site_config.rda')


## Compress Observations----
ptagis_file <- paste0('./data/CompleteTagHistories/LGR_',spp,'_',yr,'.csv')
    
observations = readCTH(ptagis_file)
    
qc_tags <- qcTagHistory(observations)
orphans <- qc_tags$orphan_tags

orphan_obs <- observations %>%
  filter(tag_code %in% orphans) %>%
  group_by(mark_site_code_value, event_site_code_value) %>%
  summarise(n = n_distinct(tag_code))
    
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
    
# DABOM is capable of fitting a model with both H and W
# if we don't want H then they need to be filtered out.
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
      filter(min_det >= start_date) %>%
      group_by(tag_code) %>%
      mutate(slot = 1:n()) %>%
      filter(!(node == 'GRA' & slot == 2 & event_type_name == 'Observation')) %>% # removes the subsequent detection in GRA ladder directly after release at trap
      group_by(tag_code) %>%
      # re-calculate the "slots" for each tag_code
      mutate(slot = 1:n()) %>%
      addDirection(parent_child = pc_nodes) %>%
      ungroup() %>%
      mutate(id = 1:n()) %>%
      select(id, everything())
          
  if(spp == 'Steelhead'){
    
    #TODO: determine lifestage, filter only spawner detections, bind with non-spawner obs so we have a complete dataset
  
      obs_direct <- steelhead_lifestage(obs_direct, spawn_year = yr, max_spawn_month = 3) %>% # confident of kelt in May
        select(id, tag_code, life_stage, everything()) 
         
      spawner_obs <- obs_direct %>%
          filter(life_stage == 'spawner') %>%
          filterDetections(parent_child = NULL,
                           max_obs_date = NULL) # could use this for to exclude kelt/repeat spawners, set to max date if left null

      dabom_obs <- obs_direct %>%
        left_join(spawner_obs %>%
                    select(id, contains('keep_obs'))) %>%
        select(id, tag_code, auto_keep_obs, user_keep_obs,
               site_code, node, direction, everything())
      # 
      # modObs <- filtered_obs %>%
      #   filter(tag_code == '3DD.003BD99F31') %>%
      #   filter(life_stage == 'spawner') %>%
      #   filter(!(direction %in% c('downstream', 'no movement'))) %>%
      #   group_by(tag_code) %>%
      #   #distinct(node, .keep_all=TRUE) %>%
      #   slice(1:which.max(node_order)) %>%
      #   mutate(dabom_obs = TRUE) %>%
      #   ungroup() %>%
      #   select(id, dabom_obs)
      # 
      # dabom_obs = filtered_obs %>%
      #   left_join(modObs,
      #             by = 'id') %>%
      #   mutate(dabom_obs = ifelse(is.na(dabom_obs), FALSE, dabom_obs)) %>%
      #   mutate(dabom_obs = case_when(
      #     user_keep_obs == FALSE ~ FALSE,
      #     is.na(user_keep_obs) ~ NA,
      #     TRUE ~ dabom_obs
      #   )) %>%
        
      # dabom_obs <- spawner_obs %>%  
      #   select(id, tag_code, life_stage, auto_keep_obs, user_keep_obs, dabom_obs, site_code, node, direction, everything())
      
  } else {
    # for Chinook
    dabom_obs <- filterDetections(compress_obs = obs_direct,   # should find final spawn location and pathway, does each obs/node fall in the pathway?  if not flag.
                                  parent_child = NULL,
                                  max_obs_date = NULL) %>%
      mutate(life_stage = 'spawner') %>%
      select(id, tag_code, life_stage, auto_keep_obs, user_keep_obs,
             site_code, node, direction, everything())
  }
    
    
# corrects calls for fish that reascended and were seen at GRA, GRS and finally GRA....we don't want these fish assigned to GRS
tmp <- dabom_obs %>%
  filter(life_stage == 'spawner') %>%
  filter(node %in% c('GRA', 'GRS')) %>%
  group_by(tag_code) %>%
  slice(which.max(min_det)) %>%
  select(tag_code, max_LGR = node)

dabom_obs <- dabom_obs %>%
  left_join(tmp) %>%
  mutate(auto_keep_obs = ifelse(is.na(auto_keep_obs), NA,
                                ifelse(node == 'GRS' & max_LGR == 'GRA', FALSE, auto_keep_obs)),
         user_keep_obs = ifelse(is.na(user_keep_obs), NA,
                                ifelse(node == 'GRS' & max_LGR == 'GRA', FALSE, user_keep_obs)))


# if its only GRS and GRA we should flag, possible re-ascension that should not be assigned to GRS: record 36
# kelt calls shouldn't be at 'no movement" directions or the same node: record 40-42
# record 98 has spawner call at GRS that occurs late april .... could be avoided with different month
# check record 124 - GRS obs after May is called a spawner, probably b/c of upstream movement at GRS?

dabom_obs %>%
  filter(life_stage == 'spawner') %>%
  filter(is.na(user_keep_obs)) %>%
  summarise(n = n_distinct(tag_code))

# if(spp == 'Steelhead'){
#   list('dabom_obs' = dabom_obs,
#        'all_obs' = all_obs) %>%
#     writexl::write_xlsx(
#       path = paste0(PITcleanrFolder,'/PreppedObs_',spp,'_',yr,'.xlsx')
#     )
#} else {
  list('dabom_obs' = dabom_obs) %>%
    writexl::write_xlsx(
      path = paste0(PITcleanrFolder,'/PreppedObs_',spp,'_',yr,'.xlsx')
    )
#}
