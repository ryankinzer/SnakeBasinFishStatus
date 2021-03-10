# Author: Kevin See
# Purpose: Generate valid tag lists for LGR
# Created: 5/1/2019
# Last Modified: 7/9/19
# Notes: 

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(lubridate)
library(PITcleanr)

#-----------------------------------------------------------------
# set species
spp = 'Steelhead'

# where is trap data?
trap_path = 'data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'

# set up folder structure
tagsFolder = 'data/ValidTagLists'
if(!dir.exists(tagsFolder)) {
  dir.create(tagsFolder)
}

yr_range <- 2020
#-----------------------------------------------------------------
# pull out valid tags from trap database, and save them in a file to upload to PTAGIS
validTags = tibble(Year = yr_range) %>%
  mutate(tags = map(Year,
                    .f = function(x) {
                      filterLGRtrapDB(trap_path = trap_path,
                                      species = spp,
                                      spawnYear = x,
                                      saveValidTagList = T,
                                      validTagFileNm = paste0(tagsFolder, '/LGR_', spp, '_', x, '.txt')) %>%
                        select(TagID = LGDNumPIT) %>%
                        distinct()
                    }))

validTags
