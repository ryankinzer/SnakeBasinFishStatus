
#------------------------------------------------------------------------------
# Read in processed observation data after final biologist call is made.
# Then summarize for final spawn location, and sex and age summaries. GSI
# comparisons can also be made with PIT-tag observations.
#------------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(PITcleanr)
library(readxl)

spp <- 'Steelhead'
yr <- 2018
#timestp <- '20190305' # for file path

#------------------------------------------------------------------------------
# Load configuration, parent_child and processed datasets from PITcleanr 
#------------------------------------------------------------------------------
load(paste0('data/PreppedData/LGR_', spp, '_', yr,'.rda'))

proc_ch = read_excel(paste0('./data/CleanedProcHist/LGR_',spp,'_EDITTED_',yr,'.xlsx')) %>%
  mutate(TrapDate = ymd(TrapDate),
         ObsDate = ymd_hms(ObsDate),
         lastObsDate = ymd_hms(lastObsDate),
         AutoProcStatus = as.logical(AutoProcStatus),
         UserProcStatus = as.logical(UserProcStatus),
         ModelObs = as.logical(ModelObs),
         ValidPath = as.logical(ValidPath)) %>%
  filter(ModelObs)


#------------------------------------------------------------------------------
# Check auto processed calls against user/biologist call.
#------------------------------------------------------------------------------
call_diff <- proc_ch %>%
  mutate(call_diff = case_when(
    AutoProcStatus != UserProcStatus ~ 1,
    TRUE ~ 0)) %>%
  group_by(TagID) %>%
  mutate(error = ifelse(sum(call_diff)>0,'error','ok'))%>%
  filter(error == 'error')

n_distinct(call_diff$TagID)

#------------------------------------------------------------------------------
# Get a quick summary of detections
#------------------------------------------------------------------------------
# number of tags observed at each node
n_tags <- proc_ch %>%
  group_by(Group, Node) %>%
  distinct(TagID, .keep_all = TRUE) %>%
  summarise(n = n()) %>%
  arrange(Group, n) 
# nodes with zero tags in proc_ch file
zero_tags <- proc_list$NodeOrder %>%
  anti_join(proc_ch, by = 'Node') %>%
  arrange(Group)
# nodes with zero tags that have observations in proc_list$ValidObs
miss_obs <- inner_join(zero_tags, proc_list$ValidObs, by = 'Node') %>%
  group_by(Group, SiteID, Node) %>%
  summarise(n = n_distinct(TagID))

#------------------------------------------------------------------------------
# Node efficiencies
#------------------------------------------------------------------------------
eff <- estNodeEff(proc_ch, proc_list$NodeOrder)


#------------------------------------------------------------------------------
# Life History Summary - for IDFG and summarizing sex, age and GSI representation
#------------------------------------------------------------------------------

lifehistory_summ = summariseTagData(capHist_proc = proc_ch,
                                    trap_data = proc_list$ValidTrapData) %>%
  mutate(equal_robots = ifelse(PtagisEventLastSpawnSite == AssignSpawnSite, TRUE,FALSE),
         equal_robots = ifelse(AssignSpawnNode == 'GRA', TRUE, equal_robots),
         equal_robots = ifelse(is.na(equal_robots), FALSE, equal_robots))

#------------------------------------------------------------------------------
# what is the agreement with IDFG
#------------------------------------------------------------------------------

# proportion of final spawn sites that disagree
length(which(!lifehistory_summ$equal_robots))/length(lifehistory_summ$equal_robots)

# number of final spawn sites that disagree
length(which(!lifehistory_summ$equal_robots))

# save file for IDFG Genetics Lab
write.csv(lifehistory_summ, 
          file = paste0('./data/LifeHistoryData/lifehistory_summ_', spp, '_', yr,'.csv'))

not_equal <- lifehistory_summ %>%
  filter(!equal_robots) %>%
  select(TagID:TagPath, contains('Ptagis'))

#------------------------------------------------------------------------------
# Save Data
#------------------------------------------------------------------------------

proc_list[["life_hist"]] <- lifehistory_summ
proc_list[["sex"]]
proc_list[["age"]]
# save entire list to feed into DABOM package
save(proc_list,
     file = paste0('data/PreppedData/LGR_', spp, '_', yr,'_',timestp,'.rda'))





#------------------------------------------------------------------------------
# Sex Ratios
#------------------------------------------------------------------------------