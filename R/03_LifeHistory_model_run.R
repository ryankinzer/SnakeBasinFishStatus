# Summarize life history data from cleaned capture history files output from
# PITcleanR.


# Load Packages
library(tidyverse)
library(readlx)

filepath <- './data/CleanedProcHist'
# files <- list.files(filepath)
# n_files <- length(files)

# Set species and Spawn Year
spp = 'Steelhead'
yr = 2010

# Load configuration files and valid tag lists.
load(paste0('data/PreppedData/LGR_', spp, '_', yr,'_', timestp,'.rda'))
# read in data
proc_ch = read_excel(paste0(filepath,'/','LGR_',spp,'_EDITTED_',yr,'.xlsx'))

#------------------------------------------------------------------------------
# Biological Summaries for IDFG and Life History
# assigns spawn location, last observation date and filters tag obs for only those
# tags used in the DABOM model.
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

# save file for John Powell
write.csv(lifehistory_summ, 
          file = paste0('./data/LifeHistoryData/lifehistory_summ_', spp, '_', yr,'_',timestp,'.csv'))

not_equal <- lifehistory_summ %>%
  filter(!equal_robots) %>%
  select(TagID:TagPath, contains('Ptagis'))


#------------------------------------------------------------------------------
# Save Data
#------------------------------------------------------------------------------

proc_list[["life_hist"]] <- lifehistory_summ

# save entire list to feed into DABOM package
save(proc_list,
     file = paste0('data/PreppedData/LGR_', spp, '_', yr,'_',timestp,'.rda'))



