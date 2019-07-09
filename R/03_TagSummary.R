
#------------------------------------------------------------------------------
# Read in processed observation data after final biologist call is made.
# Then summarize for final spawn location and IDFG genetic work.
#------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(PITcleanr)
library(readxl)

spp <- 'Steelhead'
yr <- 2018
#timestp <- '20190305' # for file path

#------------------------------------------------------------------------------
# Load configuration, parent_child and processed datasets from PITcleanr 
#------------------------------------------------------------------------------
load(paste0('data/PreppedData/LGR_', spp, '_', yr,'.rda'))

proc_ch = read_excel(paste0('./data/CleanedProcHist/LGR_',spp,'_EDITTED_',yr,'.xlsx')) %>%
  mutate(TrapDate = ymd(TrapDate),
         ObsDate = ymd_hms(ObsDate),
         lastObsDate = ymd_hms(lastObsDate)) %>%
  filter(ModelObs=='TRUE')

call_diff <- proc_ch %>%
  mutate(call_diff = case_when(
    AutoProcStatus != UserProcStatus ~ 1,
    TRUE ~ 0)) %>%
  group_by(TagID) %>%
  mutate(error = ifelse(sum(call_diff)>0,'error','ok'))%>%
  filter(error == 'error')

n_distinct(call_diff$TagID)

#------------------------------------------------------------------------------
# Get a quick summary of detections
#------------------------------------------------------------------------------
# number of tags observed at each node
n_tags <- proc_ch %>%
  group_by(Group, Node) %>%
  distinct(TagID, .keep_all = TRUE) %>%
  summarise(n = n()) %>%
  arrange(Group, n) 
# nodes with zero tags in proc_ch file
zero_tags <- proc_list$NodeOrder %>%
  anti_join(proc_ch, by = 'Node') %>%
  arrange(Group)
# nodes with zero tags that have observations in proc_list$ValidObs
miss_obs <- inner_join(zero_tags, proc_list$ValidObs, by = 'Node') %>%
  group_by(Group, SiteID, Node)# %>%
#  summarise(n = n_distinct(TagID))

parent_child = read.csv(paste0('./data/ConfigurationFiles/parent_child_',spp,'_',yr,'.csv'))
site_df <- read.csv(paste0('./data/ConfigurationFiles/site_df_',spp,'_',yr,'.csv'))

valid_paths <- getValidPaths(parent_child, 'GRA')
node_order <- createNodeOrder(valid_paths, proc_list$my_config, site_df, step_num = 3)
eff <- estNodeEff(proc_ch, node_order)
