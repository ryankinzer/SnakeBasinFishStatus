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
# set species and year
spp = 'Steelhead'
yr <- 2021

# where is trap data?
trap_path = 'data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'

# set up folder structure
tagsFolder = 'data/ValidTagLists'
if(!dir.exists(tagsFolder)) {
  dir.create(tagsFolder)
}


#-----------------------------------------------------------------
# pull out valid tags from trap database, and save them in a file to upload to PTAGIS

trap_df = read_csv(trap_path)

sppCode = ifelse(spp == 'Chinook', 1,
                 ifelse(spp == 'Steelhead', 3, NA))

# keep only correct species, spawnyear and adults (returning fish),
# as well as fish determined to be valid, with ad intact adipose fins and non-missing PIT tags
valid_df = trap_df %>%
  filter(grepl(paste0('^', sppCode), SRR),      # keep only the desired species
         SpawnYear == paste0('SY', yr),  # keep only the desired spawn year
         LGDLifeStage == 'RF',                  # keep only adults (returning fish)
         LGDValid == 1,                         # keep only records marked valid
         LGDMarkAD == 'AI',                     # keep only adipose-intact records
         !is.na(LGDNumPIT))                     # remove any records with missing PIT tag code


validTagFileNm = paste0(tagsFolder, '/LGR_', spp, '_', yr, '.txt')

tag_list <- valid_df %>%
    select(LGDNumPIT) %>%
    write.table(validTagFileNm,
                quote = F,
                row.names = F,
                col.names = F,
                sep = '\t')

# QA/QC Complete Tag History----
ptagis_file <- paste0('./data/CompleteTagHistories/LGR_',spp,'_',yr,'.csv')

qc_detections <- qcTagHistory(ptagis_file)
qc_detections
