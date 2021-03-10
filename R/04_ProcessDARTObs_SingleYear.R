# Author: Ryan N. Kinzer
# Purpose: Process LGR tag histories for DABOM using DART query
# Created: 2/12/2021
# Last Modified: 
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
#library(magrittr)
library(PITcleanr)
#library(readxl)
#library(WriteXLS)
#library(lubridate)

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
spp <- 'Steelhead'
yr <- 2020

startDate = if_else(spp == 'Steelhead',
                    paste0(yr-1, '0701'),
                    if_else(spp == 'Chinook',
                            paste0(yr, '0301'),
                            NULL))

dart_obs <- processDART_LGR(species = spp,
                            spawnYear = yr,
                            configuration = configuration,
                            truncate = T)

obs <- dart_obs$dart_obs
proc_ch <- dart_obs$proc_ch

# Run proc_ch with PTAGIS observation file.

trap_path = 'data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'

parent_child <- createParentChildDf(site_df,
                                    configuration,
                                    startDate = startDate)
 
# valid_paths <- getValidPaths(parent_child)
# 
# node_order = createNodeOrder(valid_paths,
#                              configuration = configuration,
#                              site_df,
#                              step_num = 3)

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

proc_ch <- proc_list$ProcCapHist
