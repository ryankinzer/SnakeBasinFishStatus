# Author: Kevin See
# Purpose: Summarise sex and age structure from PITcleanr
# Created: 7/1/19
# Last Modified: 7/9/19
# Notes: 

# load needed libraries
library(tidyverse)
library(sf)
library(magrittr)
library(readxl)
library(WriteXLS)
library(PITcleanr)

# set up folder structure
LifeHistory = 'data/LifeHistoryData' # for processed files
if(!dir.exists(LifeHistory)) {
  dir.create(LifeHistory)
}

# load configuration and site_df data
load('./data/ConfigurationFiles/site_config.rda')

# load function and assign group variables to all sites; TRT POPs and GSI
source('./R/assign_POP_GSI.R')

# set species
spp = 'Steelhead'

pop_ls <- assign_POP_GSI(species = spp, configuration, site_df)
grp_df <- pop_ls[[1]]
map_df <- pop_ls[[2]]

#-----------------------------------------------------------------
# take tag summaries from PITcleanr, remove duplicate tags and summarise by sex, age and brood year
for(yr in 2010:2018) {
  cat(paste('Working on', spp, 'in', yr, '\n'))
  
  load(paste0('data/DABOMready/LGR_', spp, '_', yr, '.rda'))
  
  # summary of tags
  tagSumm = summariseTagData(capHist_proc = proc_list$ProcCapHist %>%
                               filter(UserProcStatus),
                             trap_data = proc_list$ValidTrapData) %>%
    # how many records for each tag?
    group_by(TagID) %>%
    mutate(Species = spp, nRec = n()) %>%
    ungroup() %>%
    mutate_if(is.factor, as.character) %>%
    arrange(TagID, CollectionDate)
  
  # deal with duplicated tags - are these fallback and reascension tags?
  # we could deal with this be changing our join in summariseTagData()
  dupTags = tagSumm %>%
    filter(nRec > 1)
  
  #n_distinct(dupTags$TagID)
  
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
    left_join(grp_df, by = c('AssignSpawnNode' = 'Node'))

  
  # Calculate age (freshwater, saltwater, total & broody year)
  tagSumm %<>%
    mutate(BioScaleFinalAge = if_else(grepl('\\?', BioScaleFinalAge) | BioScaleFinalAge == 'N:A',
                                    as.character(NA),
                                    BioScaleFinalAge)) %>%
    mutate(fwAge = str_sub(BioScaleFinalAge, 1, 1),
           fwAge = as.integer(fwAge)) %>%
    mutate(swAge_tmp = str_remove(BioScaleFinalAge,
                                  '.*:')) %>%
    # for minijacks, replace MJ with 0
    rowwise() %>%
    mutate(swAge_tmp = if_else(Species == 'Chinook',
                               str_replace_all(swAge_tmp,
                                               'MJ', '0'),
                               swAge_tmp)) %>%
    ungroup() %>%
    # for repeat spawners, replace the S with a 1 (for one year in the ocean between spawning)
    mutate(swAge_tmp = str_replace_all(swAge_tmp,
                                       'S', '1'),
           swAge_tmp = str_replace_all(swAge_tmp,
                                       's', '1'),
           swAge_tmp = str_split(swAge_tmp, '')) %>%
    mutate(swAge = map_int(swAge_tmp,
                           .f = function(x) sum(as.integer(x)))) %>%
    select(-swAge_tmp) %>%
    mutate(totalAge = fwAge + swAge + 1) %>%
    # assign brood year
    mutate(BrdYr = as.integer(str_extract(SpawnYear, '[:digit:]+')) - totalAge)
  
  
  # xtabs(~ TRT + NWR_POPID, tagSumm, drop.unused.levels = T)
  
  #-----------------------------------------------------------------
  # join the spawn node to detectSites to bring in MPG, and summarise prop. female by model branch and MPG
  modSexDf = tagSumm %>%
    group_by(MPG, TRT, Group, GenSex) %>%
    summarise(nTags = n_distinct(TagID)) %>%
    ungroup() %>%
    filter(GenSex %in% c('F', 'M')) %>%
    spread(GenSex, nTags,
           fill = 0) %>%
    mutate(nSexed = F + M) %>%
    mutate(propF = F / (F + M),
           propF_se = sqrt((propF * (1 - propF)) / (F + M))) %>%
    select(MPG, TRT, Group, nSexed, everything()) %>%
    mutate(Group = if_else(is.na(Group), 'NotObserved', Group))
  
  #-----------------------------------------------------------------
  # summarise age props (and brood year) by model branch and MPG
  
  modAgeDf = tagSumm %>%
    filter(!is.na(totalAge)) %>%
    group_by(MPG, TRT, Group, totalAge) %>%
    summarise(nTags = n_distinct(TagID)) %>%
    ungroup() %>%
    mutate(totalAge = paste0('age', totalAge)) %>%
    spread(totalAge, nTags,
           fill = 0) %>%
    mutate(nAged = select(., -(MPG:Group)) %>% rowSums) %>%
    select(MPG, TRT, Group, nAged, everything()) %>%
    mutate(Group = if_else(is.na(Group), 'NotObserved', Group))
  
  modBrdYrDf = tagSumm %>%
    filter(!is.na(BrdYr)) %>%
    group_by(MPG, TRT, Group, BrdYr) %>%
    summarise(nTags = n_distinct(TagID)) %>%
    ungroup() %>%
    spread(BrdYr, nTags,
           fill = 0) %>%
    mutate(nAged = select(., -(MPG:Group)) %>% rowSums) %>%
    select(MPG, TRT, Group, nAged, everything()) %>%
    mutate(Group = if_else(is.na(Group), 'NotObserved', Group))
  
  list(TagSummary = tagSumm,
       SexRatio = modSexDf,
       AgeFreq = modAgeDf,
       BroodYear = modBrdYrDf) %>%
    WriteXLS(paste0(LifeHistory,'/LGR_', spp, '_', yr, '.xlsx'),
             AdjWidth = F,
             BoldHeaderRow = T,
             AutoFilter = F,
             FreezeRow = 1)
  
#  rm(tagSumm, dupTags, dupTagsKeep, dupTagsKeep1, modAgeDf, modSexDf, modBrdYrDf, proc_list, configuration, startDate, parent_child, site_df)
  
}


