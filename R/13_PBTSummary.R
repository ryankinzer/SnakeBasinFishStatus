# Purpose: Summarize hatchery PBT and pit-tag detection data for each node.
# Author: Ryan N. Kinzer
# Created: 8/19/19
# Last Modified: 7/9/19
# Notes: 

# Load Packages----
library(tidyverse)
library(lubridate)
library(PITcleanr)

# Script Variables----
spp <- 'Steelhead'
yr_range <- as.list(2010:2019)
names(yr_range) <- 2010:2019

# Combine all observation data---------------------------------------
PBTobs <- map_df(.x = yr_range,
                 .id = 'spawn_year',
                 .f = function(x){
                   
                   read_csv(paste0('./data/PBTData/HatcheryTagObservations/PBT_Complete Tag History_',
                                   x,
                                   '.csv'))
                   
                 }) # %>%
#   group_by(spawn_year) %>%
#   nest()

PBTobs %>%
  filter() %>%
  select(`Mark Site Code Value`, `Release Site Code Value`, `Release Site RKM Value`, `Release Site RKM Total`) %>%
  distinct() %>%
write_csv('./data/PBTData/PTAGIS_release_sites.csv')
# Load Configuration Data

load(paste0('./data/DABOMready/LGR_Steelhead_',yr_range[9],'.rda'))
valid_paths <- proc_list$ValidPaths
node_order <- proc_list$NodeOrder
rm(proc_list)

names(PBTobs) <- gsub("\\.", " ",names(PBTobs))


# Get Tag List-----------------------------------------------
LGR_hat_tags <- PBTobs %>% 
  # mutate(MarkFile = map(data,
  #                       function(df){
  # df %>%                      
  #filter(!str_detect(`Event Conditional Comments Name`, 'Tagged as Adult')) %>%
  filter(`Event Type Name` == 'Mark') %>%
  mutate(AdiposeFin = case_when(
    str_detect(`Event Conditional Comments Name`, 'Adipose Fin Clip') ~ 'Clipped',
    TRUE ~ 'Intact'),
    CodedWireTag = case_when(
      str_detect(`Event Conditional Comments Name`, 'Coded Wire Tag') ~ TRUE,
      TRUE ~ FALSE)) %>%
  mutate(spawn_year = as.integer(spawn_year),
         Release_Date = mdy(`Release Date MMDDYYYY`),
         release_month = month(Release_Date),
         release_year = year(Release_Date),
         brood_year = release_year - 1,
         ocean_age = spawn_year - release_year,
         age = spawn_year - brood_year) %>%
left_join(
  PBTobs %>% filter(`Event Site Code Value` == 'GRA') %>%
  mutate(TrapDate = mdy_hms(`Event Date Time Value`)) %>%
  group_by(`Tag Code`) %>%
  slice(which.min(TrapDate)) %>%
    select(`Tag Code`, TrapDate)) %>%
select(TagID = `Tag Code`, TrapDate, brood_year, release_year, release_month,
       Release_Date, spawn_year, age, ocean_age,
      contains('Mark'), `Release Site Code Value`,
       AdiposeFin, CodedWireTag, contains('Conditional'))

# Summarize Tag Counts
PBTmarks <- LGR_hat_tags %>%
  group_by(spawn_year, `Mark Site Code Value`, `Release Site Code Value`,
           AdiposeFin, brood_year, ocean_age, age) %>%
  summarise(graTags = n()) %>%
  left_join(configuration %>%
              select(SiteID, SiteName), by = c('Release Site Code Value' = 'SiteID')) %>%
  select(spawn_year, ReleaseSite = SiteName, `Release Site Code Value`, everything())


PCdat <- LGR_hat_tags %>%
  group_by(`Mark Site Code Value`, `Release Site Code Value`,
           AdiposeFin) %>%
  summarise(graTags = n(),
            `10%` = quantile(TrapDate, .1),
            `50%` = quantile(TrapDate, .5),
            `90%` = quantile(TrapDate, .9))

library(WriteXLS)
list(PBTmarks,
     PBTgrp) %>%
  WriteXLS('./data/PBTData/Steelhead_GRA_arrival_time.xlsx',
           SheetNames = list('year', 'release_grp'))



# Process Observations-----------------------------------------------

# Assign Nodes ----
valid_obs <- assignNodes(valid_tag_df = LGR_hat_tags,
                         observation = PBTobs,
                         configuration = configuration,
                         parent_child_df = parent_child,
                         truncate = T)
# Process obs.
proc_obs <- writeCapHistOutput(valid_obs,
                               valid_paths,
                               node_order,
                               save_file = FALSE)

proc_obs <- proc_obs %>%
  mutate(UserProcStatus = AutoProcStatus)

# Assign spawn locs.
PBTspawnLoc <- LGR_hat_tags %>%
  left_join(estimateSpawnLoc(proc_obs),
            by = 'TagID') %>%
  mutate(Group = as.character(Group),
         Group = ifelse(is.na(Group), 'GRA', Group))


PBTsumm <- full_join(PBTmarks,
                     PBTspawnLoc %>%
                       group_by(spawn_year, `Mark Site Code Value`, `Release Site Code Value`, 
                                AdiposeFin, AssignSpawnNode, release_year) %>%
                       summarise(nTags = n()))

# Load PBT abundance estimates----
PBTest <- read_excel('./data/PBTData/Steelhead_PBT_Abundance.xlsx',
                     sheet = 'Data',
                     skip = 1) %>%
  filter(method == 'PBT') %>%
  mutate(spawn_year = paste0('20',str_sub(`Run Year`, nchar(`Run Year`) -1, nchar(`Run Year`))),
         release_grp = )




tmp <- configuration %>% filter(grepl('North Fork Clearwater', SiteName))
