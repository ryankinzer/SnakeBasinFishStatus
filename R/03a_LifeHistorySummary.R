# Author: Kevin See
# Purpose: Summarise sex and age structure from PITcleanr
# Created: 7/1/19
# Last Modified: 7/9/19
# Notes: 

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(sf)
library(magrittr)
library(readxl)
library(WriteXLS)
library(PITcleanr)

#-----------------------------------------------------------------
# what sites are in which GSI boundaries
gsiSites = read_excel('data/GSI/SITE-NODE meta data.xlsx',
                      sheet = 'RO_org_config',
                      range = 'O5:T5749') %>%
  select(SiteID = `DABOM site`,
         Node = `DABOM node`,
         sthdGSI = steelGSI,
         sthdTRT = steelTRT,
         chnkGSI = `CHNK GSI`,
         chnkTRT = ChinookTRT) %>%
  distinct() %>%
  # fix one Chinook GSI
  mutate(chnkGSI = if_else(SiteID == 'PCM',
                           'HELLSC',
                           chnkGSI)) %>%
  # break out a couple sites into A and B nodes
  mutate(Node = recode(Node,
                       'SC2' = 'SC2B0',
                       'AFC' = 'AFCB0'))
gsiSites %<>%
  bind_rows(gsiSites %>%
              filter(Node %in% c('SC2B0', 'AFCB0')) %>%
              mutate(Node = str_replace(Node, 'B0$', 'A0')))


gsiSites %<>%
  gather(var, value, matches('GSI'), matches('TRT')) %>%
  mutate(Species = if_else(grepl('sthd', var),
                           'Steelhead',
                           if_else(grepl('chnk', var),
                                   'Chinook',
                                   as.character(NA)))) %>%
  mutate(var = str_remove(var, 'sthd'),
         var = str_remove(var, 'chnk')) %>%
  distinct() %>%
  spread(var, value)

#-----------------------------------------------------------------
# set CRS code
myCRS = 5070

# spatial file of steelhead population boundaries
sthdPops = st_read('data/ConfigurationFiles/SnakeRiver_Sthd_Pops.shp') %>%
  st_transform(crs = myCRS)

# pull out MPGs
sthdMPG = sthdPops %>%
  group_by(MPG) %>%
  summarise(nPopID = n_distinct(NWR_POPID)) %>%
  ungroup() %>%
  mutate_at(vars(MPG),
            list(fct_drop))

# pull out TRT populations
sthdTRT = sthdPops %>%
  filter(ESU_DPS == 'Snake River Basin Steelhead DPS') %>%
  group_by(NWR_POPID, NWR_NAME) %>%
  summarise_at(vars(SHAPE_AREA),
               list(sum)) %>%
  ungroup() %>%
  mutate_at(vars(NWR_POPID, NWR_NAME),
            list(fct_drop))

rm(sthdPops)

# set species
spp = 'Steelhead'
# set year
# yr = 2014

#-----------------------------------------------------------------
# take tag summaries from PITcleanr, remove duplicate tags and summarise by sex, age and brood year
for(yr in 2010:2018) {
  cat(paste('Working on', spp, 'in', yr, '\n'))
  
  load(paste0('data/DABOMready/LGR_', spp, '_', yr, '.rda'))
  
  # spatial file of detection sites
  detectSites = configuration %>%
    filter(SiteID %in% site_df$SiteID) %>%
    select(SiteID, SiteType, SiteName, RKM, Latitude, Longitude, Node) %>%
    distinct() %>%
    st_as_sf(coords = c('Longitude', 'Latitude'),
             crs = 4326) %>%
    st_transform(crs = myCRS) %>%
    st_join(sthdMPG %>%
              select(MPG)) %>%
    st_join(sthdTRT %>%
              select(-SHAPE_AREA)) %>%
    # add which main branch of DABOM site is on
    mutate(siteNode = str_remove(Node, 'A0$'),
           siteNode = str_remove(siteNode, 'B0$')) %>%
    left_join(site_df %>%
                select(siteNode = SiteID, modBranch = Step3))
  
  # summary of tags
  tagSumm = summariseTagData(capHist_proc = proc_list$ProcCapHist %>%
                               filter(UserProcStatus),
                             trap_data = proc_list$ValidTrapData) %>%
    # how many records for each tag?
    group_by(TagID) %>%
    mutate(nRec = n()) %>%
    ungroup() %>%
    arrange(TagID, CollectionDate)
  
  # deal with duplicated tags
  dupTags = tagSumm %>%
    filter(nRec > 1)
  n_distinct(dupTags$TagID)
  
  dupTagsKeep1 = dupTags %>%
    group_by(TagID) %>%
    filter(!grepl('\\?', BioScaleFinalAge),
           BioScaleFinalAge != 'N:A',
           !is.na(BioScaleFinalAge)) %>%
    slice(1) %>%
    ungroup()
  
  dupTags %>%
    anti_join(dupTagsKeep1 %>%
                select(TagID)) %>%
    group_by(TagID) %>%
    slice(1) %>%
    ungroup() %>%
    bind_rows(dupTagsKeep1) -> dupTagsKeep
  
  tagSumm %>%
    filter(nRec == 1) %>%
    bind_rows(dupTagsKeep) -> tagSumm
  
  tagSumm %<>%
    left_join(gsiSites %>%
                filter(Species == spp) %>%
                select(Node, obsGSI = GSI, TRT) %>%
                distinct(),
              by = c('AssignSpawnNode' = 'Node')) %>%
    left_join(detectSites %>%
                as_tibble() %>%
                select(Node, siteNode, modBranch, MPG, NWR_POPID) %>%
                distinct(),
              by = c('AssignSpawnNode' = 'Node')) %>%
    mutate_at(vars(MPG, NWR_POPID),
              list(as.character)) %>%
    mutate_at(vars(obsGSI, TRT, modBranch, MPG, NWR_POPID),
              list(~ if_else(AssignSpawnNode == 'GRA',
                             as.character(NA),
                             .)))
  
  # xtabs(~ TRT + NWR_POPID, tagSumm, drop.unused.levels = T)
  
  #-----------------------------------------------------------------
  # join the spawn node to detectSites to bring in MPG, and summarise prop. female by model branch and MPG
  modSexDf = tagSumm %>%
    group_by(MPG, TRT, modBranch, GenSex) %>%
    summarise(nTags = n_distinct(TagID)) %>%
    ungroup() %>%
    filter(GenSex %in% c('F', 'M')) %>%
    spread(GenSex, nTags,
           fill = 0) %>%
    mutate(nSexed = F + M) %>%
    mutate(propF = F / (F + M),
           propF_se = sqrt((propF * (1 - propF)) / (F + M))) %>%
    select(MPG, TRT, modBranch, nSexed, everything())
  
  #-----------------------------------------------------------------
  # join the spawn node to detectSites to bring in MPG, and summarise age props by model branch and MPG
  # modAgeDf = tagSumm %>%
  #   group_by(MPG, TRT, modBranch, BioScaleFinalAge) %>%
  #   summarise(nTags = n_distinct(TagID)) %>%
  #   ungroup() %>%
  #   filter(!grepl('\\?', BioScaleFinalAge),
  #          BioScaleFinalAge != 'N:A') %>%
  #   spread(BioScaleFinalAge, nTags,
  #          fill = 0) %>%
  #   mutate(nAged = select(., -(MPG:modBranch)) %>% rowSums) %>%
  #   select(MPG, TRT, modBranch, nAged, everything())
  
  # calculates total age of each fish
  modAgeDf = tagSumm %>%
    filter(!grepl('\\?', BioScaleFinalAge),
           BioScaleFinalAge != 'N:A') %>%
    mutate(fwAge = str_sub(BioScaleFinalAge, 1, 1),
           fwAge = as.integer(fwAge)) %>%
    mutate(swAge_tmp = str_remove(BioScaleFinalAge,
                                  '.*:'),
           swAge_tmp = str_replace_all(swAge_tmp,
                                       'S', '1'),
           swAge_tmp = str_replace_all(swAge_tmp,
                                       's', '1'),
           swAge_tmp = str_split(swAge_tmp, '')) %>%
    mutate(swAge = map_int(swAge_tmp,
                           .f = function(x) sum(as.integer(x)))) %>%
    select(-swAge_tmp) %>%
    mutate(totalAge = fwAge + swAge) %>%
    group_by(MPG, TRT, modBranch, totalAge) %>%
    summarise(nTags = n_distinct(TagID)) %>%
    ungroup() %>%
    mutate(totalAge = paste0('age', totalAge)) %>%
    spread(totalAge, nTags,
           fill = 0) %>%
    mutate(nAged = select(., -(MPG:modBranch)) %>% rowSums) %>%
    select(MPG, TRT, modBranch, nAged, everything())
  
  # calculates brood year of each fish
  modBrdYrDf = tagSumm %>%
    filter(!grepl('\\?', BioScaleFinalAge),
           BioScaleFinalAge != 'N:A') %>%
    mutate(fwAge = str_sub(BioScaleFinalAge, 1, 1),
           fwAge = as.integer(fwAge)) %>%
    mutate(swAge_tmp = str_remove(BioScaleFinalAge,
                                  '.*:'),
           swAge_tmp = str_replace_all(swAge_tmp,
                                       'S', '1'),
           swAge_tmp = str_split(swAge_tmp, '')) %>%
    mutate(swAge = map_int(swAge_tmp,
                           .f = function(x) sum(as.integer(x)))) %>%
    select(-swAge_tmp) %>%
    mutate(totalAge = fwAge + swAge) %>%
    mutate(BrdYr = as.integer(str_extract(SpawnYear, '[:digit:]+')) - totalAge) %>%
    group_by(MPG, TRT, modBranch, BrdYr) %>%
    summarise(nTags = n_distinct(TagID)) %>%
    ungroup() %>%
    spread(BrdYr, nTags,
           fill = 0) %>%
    mutate(nAged = select(., -(MPG:modBranch)) %>% rowSums) %>%
    select(MPG, TRT, modBranch, nAged, everything())
  
  
  
  list(TagSummary = tagSumm,
       SexRatio = modSexDf,
       AgeFreq = modAgeDf,
       BroodYear = modBrdYrDf) %>%
    WriteXLS(paste0('data/output/tagSummaries/LGR_', spp, '_', yr, '.xlsx'),
             AdjWidth = T,
             BoldHeaderRow = T,
             AutoFilter = F,
             FreezeRow = 1)
  
  rm(tagSumm, dupTags, dupTagsKeep, dupTagsKeep1, modAgeDf, modSexDf, modBrdYrDf, proc_list, configuration, startDate, parent_child, site_df, detectSites)
  
}


