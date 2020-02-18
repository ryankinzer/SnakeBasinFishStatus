# Purpose: Download hatchery PIT-detections at Lower Granite dam for all years.
# Author: Ryan N. Kinzer
# Created: 8/19/19
# Last Modified:
# Notes: 

# Load Packages----
library(tidyverse)
library(lubridate)
library(PITcleanr)

# Script Variables----
spp <- 'Steelhead'
yr_range <- as.list(2010:2019)
names(yr_range) <- 2010:2019

# Download Data ---------------------------------------------------
  map(.x = yr_range,
      .f = function(.x){ 
        
        adultLGR <- queryPITtagAdult(
                            site = 'GRA',
                            species = spp,
                            spawn_yr = .x)
                   
                   hatTags <- adultLGR %>%  
                     mutate(Release_Year = year(ymd(Release_Date))) %>%
                     filter(str_detect(SpRRT, 'H')) %>% # remove wild fish
                     filter(Return_Year != Release_Year) %>%# remove newly tagged adults returning
                     distinct(Tag_ID) %>%
                     select(Tag_ID)
                   
                   hatTags %>%
                     write.table(paste0('./data/PBTData/TagList/LGR_Hatchery_',spp,'_TagList_',.x,'.csv'),
                                 quote = F,
                                 row.names = F,
                                 col.names = F,
                                 sep = '\t')
                   }
)

  





