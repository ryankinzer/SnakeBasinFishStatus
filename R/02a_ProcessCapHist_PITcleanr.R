# Author: Kevin See
# Purpose: Clean up LGR tag histories for DABOM
# Created: 5/1/2019
# Last Modified: 7/9/19
# Notes: Rick Orme provided feedback on PITcleanr output on 7/1/19

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(magrittr)
library(PITcleanr)
library(readxl)
library(WriteXLS)
library(lubridate)

#-----------------------------------------------------------------
# set species
spp = 'Steelhead'

# set up folder structure
PITcleanrFolder = 'data/PITcleanr'
if(!dir.exists(PITcleanrFolder)) {
  dir.create(PITcleanrFolder)
}

DABOMdataFolder = 'data/DABOMready'
if(!dir.exists(DABOMdataFolder)) {
  dir.create(DABOMdataFolder)
}


#-----------------------------------------------------------------
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


# where is trap data?
trap_path = 'data/tblLGDMasterCombineExportJodyW.csv'

#-----------------------------------------------------------------
# read observations from PTAGIS and process with PITcleanr
for(yr in 2010:2018) {
  
  cat(paste('Starting year', yr, '\n'))
  
  # start date is July 1 of the previous year for steelhead, or March 1 of current year for Chinook
  startDate = if_else(spp == 'Steelhead',
                      paste0(yr-1, '0701'),
                      if_else(spp == 'Chinook',
                              paste0(yr, '0301'),
                              NULL))
  
  # build parent-child table
  parent_child = createParentChildDf(site_df,
                                     configuration,
                                     startDate = startDate)
  
  # get raw observations from PTAGIS
  # These come from running a saved query on the list of tags to be used
  observations = read_csv(paste0('data/CompleteTagHistories/LGR_', spp, '_', yr, '.csv'))
    
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
  
  rm(startDate, parent_child, observations, proc_list)
  
}

