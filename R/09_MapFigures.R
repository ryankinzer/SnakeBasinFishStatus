# Author: Ryan N. Kinzer
# Purpose: Produce report map based figures for results.
# Created: 9/30/19
# Last Modified: 
# Notes: 

#-----------------------------------------------------------------

# load needed libraries
library(tidyverse)
library(sf)

# load function and assign group variables to all sites; TRT POPs and GSI
source('./R/assign_POP_GSI.R')

# set species
spp = 'Steelhead'
# set year
yr = 2018

# mapping data
pop_ls <- assign_POP_GSI(species = spp, configuration, site_df)
grp_df <- pop_ls[[1]]
map_df <- pop_ls[[2]]

#site_labs <- distinct(grp_df, SiteID)

# left_join(map_df, modSexDf, by = c('TRT_POPID' = 'TRT')) %>%
#  ggplot() +
#    geom_sf(aes(fill = propF)) +
#    scale_fill_continuous(type = 'viridis') +
# #   geom_sf(data = grp_df) +
#    geom_sf_text(aes(label = round(propF,2)), size = 3)

