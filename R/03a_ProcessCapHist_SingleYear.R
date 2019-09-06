#------------------------------------------------------------------------------
# Script runs PITcleanr to gathering and clean DABOM data.
#
# Author: Ryan Kinzer and Rick Orme
#------------------------------------------------------------------------------
# Load Packages and Functions
#------------------------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(magrittr)
library(PITcleanr)
library(readxl)
library(WriteXLS)
library(lubridate)


#---- set up folder structure
PITcleanrFolder = 'data/PITcleanr'
if(!dir.exists(PITcleanrFolder)) {
  dir.create(PITcleanrFolder)
}

#------------------------------------------------------------------------------
# Identify which species and spawn year data is needed. 
#------------------------------------------------------------------------------

#timestp <- gsub('[^0-9]','', Sys.Date())

spp = 'Steelhead'  # either 'Chinook' or 'Steelhead'
yr = 2019        # tagging operations started in spawn year 2009

#------------------------------------------------------------------------------
# Download all PTAGIS interrogation and MRR sites, antennas and configuations,
#  then it appends Node names.
#------------------------------------------------------------------------------

org_config = buildConfig()

#------------------------------------------------------------------------------
# Change and customize the configuration file.
#------------------------------------------------------------------------------

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

#------------------------------------------------------------------------------
# Create Lower Granite site network and available paths
#------------------------------------------------------------------------------

site_df = writeLGRNodeNetwork()

# remove some sites that have been combined with others (see the modifications to the configuration file)

# May no longer be needed
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

# where is trap data?
trap_path = 'data/TrappingDbase/tblLGDMasterCombineExportJodyW.csv'

# start date is July 1 of the previous year for steelhead, or March 1 of current year for Chinook
startDate = if_else(spp == 'Steelhead',
                    paste0(yr-1, '0701'),
                    if_else(spp == 'Chinook',
                            paste0(yr, '0301'),
                            NULL))


parent_child = createParentChildDf(site_df,
                                  configuration,
                                   startDate = startDate)

#get observations from PTAGIS
observations = read_csv(paste0('data/CompleteTagHistories/LGR_', spp, '_', yr, '.csv'))

#observations$`Antenna ID` <- str_pad(observations$`Antenna ID`, 2, pad = '0')

#------------------------------------------------------------------------------
# Combines node names with raw PTAGIS observations and runs cleaning/processing
# algorithms
#------------------------------------------------------------------------------

proc_list = processCapHist_LGD(species = spp,
                               spawnYear = yr,
                               configuration = configuration,
                               trap_path = trap_path,
                               parent_child = parent_child,
                               observations = observations,
                               site_df = site_df,
                               truncate = T,
                               step_num = 3,
                               save_file = T,
                               file_name = paste0(PITcleanrFolder, '/LGR_', spp, '_', yr, '.xlsx'))

# save list object as data frame for saving
proc_ch <- proc_list$ProcCapHist %>%
  as_tibble()

#------------------------------------------------------------------------------
# The code below is all wrapped in the processCapHist_LGD function above.  I 
# have listed out each step below to assist in debugging, and to improve the
# understanding of the black-box.
#------------------------------------------------------------------------------
# cat("Constructing valid pathways\n")
# valid_paths = getValidPaths(parent_child)
# 
# if (class(valid_paths) == "character") {
#   print(paste("The following nodes returned an error:", 
#               paste(valid_paths, collapse = ", ")))
#   return(NULL)
# }
# 
# cat("Filtering trap database\n")
# trap_df = filterLGRtrapDB(trap_path, species, spawnYear)
# if (filter_by_PBT) {
#   trap_df = trap_df %>% filter(!grepl("H$", SRR))
# }
# 
# cat("Getting valid tags\n")
# valid_tag_df = trap_df %>% group_by(TagID = LGDNumPIT) %>% 
#   summarise(TrapDate = min(CollectionDate, na.rm = T))
# 
# cat("Assigning nodes\n")
# valid_obs = assignNodes(valid_tag_df, observations, configuration, 
#                         parent_child, truncate)
# 
# cat("Creating node order\n")
# node_order = createNodeOrder(valid_paths = valid_paths, configuration = configuration, 
#                              site_df = site_df, step_num = step_num)
# 
# cat("Processing assigned nodes\n")
# save_df = writeCapHistOutput(valid_obs, valid_paths, node_order, 
#                              last_obs_date, save_file, file_name)

#------------------------------------------------------------------------------
# Save Data for biologist review
#------------------------------------------------------------------------------

write.csv(proc_ch, 
     file = paste0('./data/PreppedData/ProcCapHist_NULL_', spp, '_', yr,'.csv'),
     row.names = FALSE)

save(spp, yr, startDate, site_df, configuration, parent_child, proc_list, 
     file = paste0('./data/PreppedData/LGR_', spp, '_', yr,'.rda'))

# After biologist review
spp <- 'Steelhead'
yr <- 2019

DABOMdataFolder = 'data/DABOMready'
if(!dir.exists(DABOMdataFolder)) {
  dir.create(DABOMdataFolder)
}

load(paste0('./data/PreppedData/LGR_', spp, '_', yr,'.rda'))

# for Chinook, remove some observations from basins we don't expect Chinook to go to
if(spp == 'Chinook') {
  proc_list$ProcCapHist %<>%
    filter(!Group %in% c('Asotin', 'Lapwai', 'Potlatch', 'JosephCreek', 'CowCreek', 'CarmenCreek', 'Almota', 'Alpowa', 'Penawawa'))
  
}

# replace default PITcleanr output with files from Rick Orme
proc_list$ProcCapHist = read_excel(paste0('data/CleanedProcHist/LGR_', spp, '_EDITTED_', yr, '.xlsx')) %>%
  mutate_at(vars(AutoProcStatus, UserProcStatus, ModelObs, ValidPath),
            list(as.logical)) %>%
  # filter(UserProcStatus) %>%
  select(one_of(names(proc_list$ProcCapHist))) %>%
  mutate_at(vars(TrapDate),
            list(ymd)) %>%
  # mutate_at(vars(TrapDate),
  #           list(as.POSIXct),
  #           tz = "America/Los_Angeles") %>%
  mutate_at(vars(ObsDate, lastObsDate),
            list(ymd_hms)) %>%
  mutate_at(vars(Group),
            list(factor),
            levels = levels(proc_list$ProcCapHist$Group)) %>%
  arrange(TagID, ObsDate)


# save some stuff
save(spp, yr, startDate, site_df, configuration, parent_child, proc_list,
     file = paste0(DABOMdataFolder, '/LGR_', spp, '_', yr, '.rda'))
