# Testing script to learn and debug new PITcleanR package.
devtools::install_github("BiomarkABS/PITcleanr", ref = "develop", build_vignettes = T)

library(tidyverse)
library(PITcleanr)

browseVignettes("PITcleanr")

# example ptagis file doesn't work
# ptagis_file <- read.csv(system.file("extdata", 
#                                     "LGR_Chinook_2014.csv", 
#                                     package = "PITcleanr",
#                                     mustWork = TRUE))
# comp_obs <- compress(ptagis_file)


# Set species and year
spp <- 'Steelhead'
yr = 2018

# load ptagis complete tag history file and prep for functions
#ptagis_file <- read_csv(paste0('data/CompleteTagHistories/LGR_', spp, '_', yr, '.csv'))
# Need: Tag, Mark Species, Mark Rear Type, Event Type, Event Site Type, Event Site Code, Event Date Time, Antenna, Antenna Group Configuration, Event Release Site Code, Event Release Date Time

# Error with compress because file names need to be lower snake case 
#names(ptagis_file) <- gsub(' ','_', tolower(names(ptagis_file)))

# Error with compress because file names need to be lower snake case 
# ptagis_file <- ptagis_file %>%
#   mutate(#across(.cols = c('event_date_mmddyyyy', 'event_release_date_mmddyyyy'), .fns = lubridate::parse_date_time, '%m/%d/%Y %I:%M:%S %p'),
#     across(.cols = c('event_date_time_value', 'event_release_date_time_value'), .fns = lubridate::parse_date_time, '%m/%d/%Y %I:%M:%S %p'),
#     mark_species_name = spp,
#     event_site_type_description = 'INT')


# Or read the file and format directly
# ptagis_file <- readCTH(paste0('data/CompleteTagHistories/LGR_', spp, '_', yr, '.csv')) %>%
#   mutate(mark_species_name = spp)




# get configuration file
array_configuration <- buildConfig()
#ptagis_meta <- queryPtagisMeta()

site_meta <- array_configuration %>%
  select(contains('site'), contains('rkm'), latitude, longitude) %>%
  distinct() %>%
  mutate(detection_group = case_when(
    as.numeric(str_extract(rkm, '\\d+')) <= 234 ~ 'Lower_C',
    as.numeric(str_extract(rkm, '\\d+')) < 522 ~ 'Mid_C',
    as.numeric(str_extract(rkm, '\\d+')) > 522 ~ 'Upper_C',
    as.numeric(str_extract(rkm, '\\d+')) == 522 ~ 'Above_GRA')) %>%
  mutate(detection_group = case_when(
    detection_group == 'Above_GRA' & as.numeric(sub('.','',str_extract(rkm, '\\.\\d+'))) < 173 ~ 'Below_GRA',
      TRUE ~ detection_group))

configuration <- array_configuration %>%
  left_join(site_meta %>%
              select(site_code, detection_group),
            by = 'site_code') %>%
  mutate(node = ifelse(detection_group == 'Above_GRA', node, detection_group))

# customize some nodes based on DABOM framework
configuration = array_configuration %>%
  mutate(node = ifelse(site_code %in% c('VC2', 'VC1'),
                       site_code,
                       node),
         node = ifelse(site_code == 'SC2',
                       'SC2B0',
                       node),
         node = ifelse(site_code %in% c('CROTRP',
                                        'CRT',
                                        'REDTRP',
                                        'REDR',
                                        'RRT'),
                       'SC2A0',
                       node),
         node = ifelse(node == 'ACB',
                       'ACBB0',
                       node),
         node = ifelse(node == 'CCA',
                       'CCAB0',
                       node),
         node = ifelse(site_code == 'AFC',
                       ifelse(grepl('MAINSTEM', antenna_group),
                              'AFCB0',
                              'AFCA0'),
                       node),
         node = ifelse(site_code == 'HBC',
                       'HYCA0',
                       node),
         node = ifelse(site_code == 'MCCA',
                       'STR',
                       node),
         node = ifelse(site_code == 'CARMEC',
                       'CRCA0',
                       node),
         node = ifelse(site_code == 'BIG2C',
                       'TAYA0',
                       node),
         node = ifelse(site_code == 'WIMPYC',
                       'WPCA0',
                       node),
         node = ifelse(site_code == 'IML' & config_id == 130 & antenna_id == '09',
                       'IMLA0',
                       node),
         node = str_replace(node, '^BTC', 'BTL'),
         node = ifelse(site_code %in% c('YANKFK', 'CEY'),
                       'YFKA0',
                       node),
         node = ifelse(site_code == 'SAWT',
                       'STL',
                       node),
         node = ifelse(site_code == 'LOOH',
                       'LOOKGC',
                       node),
         node = ifelse(site_code == 'RPDTRP',
                       'RAPH',
                       node),
         node = ifelse(site_code == 'CHARLC',
                       'CCAB0',
                       node),
         node = ifelse(node == 'KEN',
                       'KENB0',
                       node),
         node = ifelse(node == 'HYC',
                       'HYCB0',
                       node),
         node = ifelse(node == 'YFK',
                       'YFKB0',
                       node),
         node = ifelse(node == 'LLR',
                       'LLRB0',
                       node),
         node = ifelse(node == 'LRW',
                       'LRWB0',
                       node),
         node = ifelse(site_code == '18M',
                       str_replace(node, '18M', 'HEC'),
                       node),
         node = ifelse(node == 'LAP',
                       'LAPB0',
                       node),
         node = ifelse(site_code == 'KRS',
                       'KRS',
                       node)) %>%
  distinct()

# get parent-child file
# parent_child_file = system.file("extdata",
#                                 "parent_child_GRA.csv",
#                                 package = "PITcleanr",
#                                 mustWork = TRUE)

parent_child <- read_csv('./data/ConfigurationFiles/parent_child.csv')
plotNodes(parent_child,
          layout = "tree")

parent_child_node = addParentChildNodes(parent_child,
                                        configuration)
plotNodes(parent_child_node,
          layout = "tree")


# Down load data from DART - includes complete tag history after fish was observed at LGR 
comp_dart <- compressDART(
  species = spp,
  loc = 'GRA',
  spawn_year = yr,
  configuration = configuration
)

comp_obs <- comp_dart$compress_obs






# check the file for orphan and disown tags
qc_detections <- qcTagHistory(ptagis_file)

comp_obs <- compress(ptagis_file, configuration = configuration, units = 'days')


# determine trap date, and remove detections prior to that
obs_direct = comp_obs %>%
  # get the first detection of each tag at Lower Granite Dam
  left_join(comp_obs %>%
              filter(node == "GRA",
                     event_type_name %in% c("Mark", "Recapture")) %>%
              group_by(tag_code) %>%
              filter(min_det == min(min_det)) %>%
              summarise(start_date = min_det,
                        .groups = "drop"),
            by = "tag_code") %>%
  # filter any detections before the "start_date"
  filter(min_det >= start_date) %>%
  group_by(tag_code) %>%
  # re-calculate the "slots" for each tag_code
  mutate(slot = slot - min(slot) + 1) %>%
  ungroup() %>%
  # add direction using "addDirection()
  addDirection(parent_child = parent_child_node)


if(spp == 'Steelhead'){
  min_date = paste0(yr,'0701')
} else {
  min_date = paste0(yr,'0301')
}

sites_sf <- extractSites(ptagis_file,
                         as_sf = TRUE,
                         min_date = min_date,
                         configuration = configuration)
# focus on sites within Snake River Basin
sites_sf = sites_sf %>%
  # all sites in the Snake have a river kilometer that starts with 754
  filter(grepl("522.", rkm))

load('../../DFRM Projects/River_Mapping/data/flowlines/large_rivers.rda')

nhd_list = queryFlowlines(sites_sf = sites_sf,
                          root_site_code = "GRA",
                          min_strm_order = 2,
                          dwnstrm_sites = T,
                          dwn_min_stream_order_diff = 2)

ggplot() +
  geom_sf(data = snake_rivers, inherit.aes = FALSE) +
  geom_sf(data = sites_sf, inherit.aes = FALSE) +
  geom_sf_label(data = sites_sf, aes(label = site_code), size = 2)




# compress function now works if we run the above 2 code sections
# keeping the fish info and batch release tagging info would be good
#comp_obs <- compress(ptagis_file, units = 'days')

#ptagis_meta <- queryPtagisMeta() # is any information lost with a select call?