# Purpose: Helper function to assign MPG and TRTs to Sites and Nodes.
# Author: Ryan Kinzer
# Created: 7/12/19
# Modified:

site_trt_designations <- function(species, configuration){

  load('./data/ConfigurationFiles/DABOM_map_data.rda')

  config <- configuration %>%
    filter(SiteID %in% site_df$SiteID) %>%
    select(SiteID, SiteType, SiteName, RKM, Latitude, Longitude, Node) %>%
    distinct() %>%
    sf::st_as_sf(coords = c('Longitude', 'Latitude'),
                 crs = 4326)
  
  if(species == 'Steelhead'){
    df <- config %>%
      sf::st_join(SR_st_pop %>%
            select(ESU_DPS, MPG, POP_NAME, TRT = TRT_POPID)) %>%
      as_tibble() %>%
      select(-geometry)
  }
  
  if(species == 'Chinook'){
    df <- config %>%
      sf::st_join(SR_ch_pop %>%
                    select(ESU_DPS, MPG, POP_NAME, TRT = TRT_POPID)) %>%
      as_tibble() %>%
      select(-geometry)
                
  }
  
  return(df)

}