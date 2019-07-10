#------------------------------------------------------------------------------
# Script runs PITcleanr to gathering and clean DABOM data.
#
# Author: Ryan Kinzer and Rick Orme
#------------------------------------------------------------------------------
# Load Packages and Functions
#------------------------------------------------------------------------------
library(tidyverse)
library(PITcleanr)
source('./R/loadLGTrappingDBase.R')

#------------------------------------------------------------------------------
# Identify which species and spawn year data is needed. 
#------------------------------------------------------------------------------

#timestp <- gsub('[^0-9]','', Sys.Date())

spp = 'Steelhead'  # either 'Chinook' or 'Steelhead'
yr = 2016        # tagging operations started in spawn year 2009

#------------------------------------------------------------------------------
# Download all PTAGIS interrogation and MRR sites, antennas and configuations,
#  then it appends Node names.
#------------------------------------------------------------------------------

org_config = buildConfig()

#------------------------------------------------------------------------------
# Change and customize the configuration file.
#------------------------------------------------------------------------------

my_config = org_config %>%
  mutate(Node = ifelse(SiteID %in% c('VC2', 'VC1', 'LTR', 'MTR', 'UTR'),
                       SiteID,
                       Node),
         Node = ifelse(SiteID %in% c('CROTRP',
                                     'CRT',
                                     'REDTRP',
                                     'REDR',
                                     'RRT'),
                       'above_SC2',
                       Node),
         Node = ifelse(SiteID == 'ACB', 'ACB', Node),
         Node = ifelse(SiteID == 'AFC',
                       ifelse(grepl('MAINSTEM', AntennaGroup),
                              'AFCB0',
                              'AFCA0'),
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
         Node = ifelse(Node == 'LLR',
                       'LLRB0',
                       Node),
         Node = ifelse(Node == 'LRW',
                       'LRWB0',
                       Node),
         Node = ifelse(SiteID == '18M',
                       paste0('X', Node),
                       Node)) %>%
  distinct()

#------------------------------------------------------------------------------
# Save customized configuration file as .csv with time stamp.
#------------------------------------------------------------------------------
my_config <- my_config %>% filter(Node != 'IML')

write.csv(my_config, file = paste0('./data/ConfigurationFiles/my_config_',spp,'_',yr,'.csv'),
          row.names = FALSE)

#------------------------------------------------------------------------------
# Create Lower Granite site network and available paths
#------------------------------------------------------------------------------

site_df = writeLGRNodeNetwork()

# remove some sites that have been combined with others (see the modifications to the configuration file)

# May no longer be needed
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

#--------------------------
# Double check SiteId names
#--------------------------
anti_join(site_df, my_config, by = 'SiteID') # zero observations is good!
#not_in_myconfig <- anti_join(my_config, site_df, by = 'SiteID') # lots of observations
#write.csv(not_in_myconfig, file = './not_in_myconfig.csv')

# Save Output
write.csv(site_df, file = paste0('./data/ConfigurationFiles/site_df_',spp,'_',yr,'.csv'),
          row.names = FALSE)

#------------------------------------------------------------------------------
# Create parent-child node network based on site_df and customized
# configuration.
#
# Term: Node Detection Landscape 
#
# Joins DABOM site hierarchy with config table to get node names, then it 
# only keeps node available and needed for DABOM model
#
# dates are used to identify which antenna configuration in PTAGIS to use
#
# could be maintained externally
#------------------------------------------------------------------------------

parent_child = createParentChildDf(site_df,
                                   my_config,
                                   startDate = ifelse(spp == 'Chinook',
                                                      paste0(yr, '0301'),
                                                      paste0(yr-1, '0701')))

write.csv(parent_child, 
     file = paste0('./data/ConfigurationFiles/parent_child_',spp,'_',yr,'.csv'),
     row.names = FALSE)

#------------------------------------------------------------------------------
## Get valid tags from Lower Granite trapping database and relavent fields,
# and saves a .txt file of tag codes for upload into PTAGIS complete tag 
# history query
#------------------------------------------------------------------------------

# Export .csv file from access database
#trapDB_filepath <- './data/TrappingDBase/LGTrappingExportJodyW.accdb'
#con <- loadLGTrappingDBase(trapDB_filepath)

#trap_df <- DBI::dbReadTable(con, 'tblLGDMasterCombineExportJodyW')
#write_csv(trap_df, path = 'data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv')

trap_filepath = './data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'

valid_df = filterLGRtrapDB(trap_path = trap_filepath,
                           species = spp,
                           spawnYear = yr,
                           saveValidTagList = T,
                           validTagFileNm = paste0('data/ValidTagLists/LGR_', spp, '_', yr, '.txt'))

#------------------------------------------------------------------------------
# Read in the result of the PTAGIS complete tag history query
#------------------------------------------------------------------------------
observations = read_csv(paste0('data/CompleteTagHistories/LGR_', spp, '_', yr, '.csv'))

observations$`Antenna ID` <- str_pad(observations$`Antenna ID`, 2, pad = '0')

#------------------------------------------------------------------------------
# Combines node names with raw PTAGIS observations and runs cleaning/processing
# algorithms
#------------------------------------------------------------------------------

proc_list = processCapHist_LGD(species = spp,
                               spawnYear = yr,
                               configuration = my_config,
                               parent_child = parent_child,
                               trap_path = trap_filepath,
                               filter_by_PBT = T,
                               observations = observations,
                               truncate = T,
                               site_df = site_df,
                               step_num = 3,
                               save_file = F,
                               file_name = paste0('data/PreppedData/ProcCapHist_', spp, '_', yr, '.csv'))

# save list object as data frame for saving
proc_ch <- proc_list$ProcCapHist %>%
  as_tibble()

#------------------------------------------------------------------------------
# The code below is all wrapped in the processCapHist_LGD function above.  I 
# have listed out each step below to assist in debugging, and to improve the
# understanding of the black-box.
#------------------------------------------------------------------------------
# cat("Constructing valid pathways\n")
# valid_paths = getValidPaths(parent_child)
# 
# if (class(valid_paths) == "character") {
#   print(paste("The following nodes returned an error:", 
#               paste(valid_paths, collapse = ", ")))
#   return(NULL)
# }
# 
# cat("Filtering trap database\n")
# trap_df = filterLGRtrapDB(trap_path, species, spawnYear)
# if (filter_by_PBT) {
#   trap_df = trap_df %>% filter(!grepl("H$", SRR))
# }
# 
# cat("Getting valid tags\n")
# valid_tag_df = trap_df %>% group_by(TagID = LGDNumPIT) %>% 
#   summarise(TrapDate = min(CollectionDate, na.rm = T))
# 
# cat("Assigning nodes\n")
# valid_obs = assignNodes(valid_tag_df, observations, configuration, 
#                         parent_child, truncate)
# 
# cat("Creating node order\n")
# node_order = createNodeOrder(valid_paths = valid_paths, configuration = configuration, 
#                              site_df = site_df, step_num = step_num)
# 
# cat("Processing assigned nodes\n")
# save_df = writeCapHistOutput(valid_obs, valid_paths, node_order, 
#                              last_obs_date, save_file, file_name)

#------------------------------------------------------------------------------
# Save Data
#------------------------------------------------------------------------------
write.csv(proc_ch, 
     file = paste0('./data/PreppedData/ProcCapHist_NULL_', spp, '_', yr,'.csv'),
     row.names = FALSE)

# add items to proc_list
proc_list$parent_child <- parent_child
proc_list$site_df <- site_df

# save entire list to feed into DABOM package
save(proc_list,
     file = paste0('data/PreppedData/LGR_', spp, '_', yr,'.rda'))

