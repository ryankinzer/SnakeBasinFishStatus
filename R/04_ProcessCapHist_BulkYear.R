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
# set up folder structure
PITcleanrFolder = 'data/PITcleanr'
if(!dir.exists(PITcleanrFolder)) {
  dir.create(PITcleanrFolder)
}

DABOMdataFolder = 'data/DABOMready'
if(!dir.exists(DABOMdataFolder)) {
  dir.create(DABOMdataFolder)
}

# load configuration and site_df data
load('./data/ConfigurationFiles/site_config.rda')

#-----------------------------------------------------------------
# read observations from PTAGIS and process with PITcleanr
# set species
spp = 'Chinook'
# have the processed capture history files been cleaned by the 
# biologist and are they ready for DABOM
bio_cleaned = TRUE
# where is trap data?
trap_path = 'data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'

for(yr in 2010:2019) {
  
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
  
  if(bio_cleaned){
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
      
  }
  
  rm(startDate, parent_child, observations, proc_list)
  
}

