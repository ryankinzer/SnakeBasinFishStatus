# Testing script to learn and debug new PITcleanR package.

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
spp <- 'Chinook'
yr = 2021

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


# load configuration and site_df data
load('./data/ConfigurationFiles/site_config.rda')

# Down load data from DART - includes complete tag history after fish was observed at LGR 
comp_dart <- compressDART(
  species = spp,
  loc = 'GRA',
  spawn_year = yr,
  configuration = configuration
)

comp_obs <- comp_dart$compress_obs
dart_obs <- comp_dart$dart_obs

lgrldr <-  dart_obs %>%
  filter(mark_site == 'LGRLDR',
         mark_rear_type_name == 'W') %>%
  distinct(tag_code) %>% pull(tag_code)

lgrldr_obs <- comp_obs %>%
  filter(tag_code %in% lgrldr)

tmp <- lgrldr_obs %>% group_by(node) %>% summarise(n = n_distinct(tag_code))


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