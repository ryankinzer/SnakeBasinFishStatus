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
yr = 2018        # tagging operations started in spawn year 2009

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
#
# retrieves all specified species and spawn year,
#                        returning adults (LGDLifeStage = RF)
#                        LGDValid = 1,
#                         LGDMarkAD = AI
#                       fish with tag codes.
#
#  AND ONLY INCLUDES SRR = 15X (trap call not PTAGIS mark call) 
#       if species = 'Chinook';
#    spring/summer chinook catch at trap MUST BE CALLED SRR = '15X'
#------------------------------------------------------------------------------

trapDB_filepath <- './data/LGTrappingExportJodyW.accdb'
con <- loadLGTrappingDBase(trapDB_filepath)

trap_df <- DBI::dbReadTable(con, 'tblLGDMasterCombineExportJodyW')
write_csv(trap_df, path = 'data/tblLGDMasterCombineExportJodyW.csv')

trap_filepath = './data/tblLGDMasterCombineExportJodyW.csv'

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
#
# returns a list with 5 objects
# 1. ValidPaths - path of fish from GRA to each Node; created with parent child table and 
#                   the getValidPaths() functions
# 2. NodeOrder - metadata for each Node - possible DABOM model information and parameter settingss
#                     runs from createNodeOrder(), writeLGRNodeNetwork(), createParentChildDf control
#                     much of the DABOM groups and model landscape
# 3. ValidTrapData - LowerGraniteDB for the correct species and spawnyear and no 'H' in SRR
# 4. ValidObs - duplicates in ValidTrapData removed and trap date is min(CollectionDate),
#                 LGR trap data is combined with observation data and Nodes are assigned,
#                 observations are also truncated to only valid nodes and obs dates and 
#                 the extra detections at a site (min obs date within an uninterupted sequence by node.)
# 5. ProcCapHist - from writeCapHistOutput(), has processed observations using two seperate
#               cleaning algorithms, direction of travel and migration direction
#
# Fix/add max date to ProcCapHist object which comes from writeCapHistOutput().
#    need Obsdate to be minObsdate at the node and need to include a 
#    maxObsdate at the node for each detection sequence.
#
# Fix ModelObs call - which comes from the writeSpwnPaths or assignNodes
#     need to flag all fish seen in two branches
#
# Fix maxUpDate field from writeSpwnPaths() could be incorrect if a fish moves downstream
# and then goes partially back up, maybe just flag observations for reveiw
#
# Fix truncate = F, currently doesn't run.
#
# Fix save_file portion if we think its needed. Change in function to save as .csv
#
# Fix filter_by_PBT to include/exclude wild, hnc
#
# step_num = points to main branch column in site_df - NEEDS TO BE 3
# filter_by_PGT = removes all 'XXH' from SRR field in Lower Granite DB
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
                               file_name = paste0('data/ProcCapHist_', spp, '_', yr, '.csv'))

# save list object as data frame for saving
proc_ch <- proc_list$ProcCapHist %>%
  as_tibble()

#------------------------------------------------------------------------------
# Save Data
#------------------------------------------------------------------------------
write.csv(proc_ch, 
     file = paste0('./data/PreppedData/ProcCapHist_NULL_', spp, '_', yr,'.csv'),
     row.names = FALSE)

# save entire list to feed into DABOM package
save(proc_list,
     file = paste0('data/PreppedData/LGR_', spp, '_', yr,'.rda'))

