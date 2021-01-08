#------------------------------------------------------------------------------
# Script gathers data and runs the STADEM model.  A loop is built in to cycle
# through Chinook and steelhead and all years begining in 2010.
#
# Author: Ryan Kinzer
#------------------------------------------------------------------------------

# Download STADEM Package
if(!require(STADEM)){
  devtools::install_github('KevinSee/STADEM', ref = 'master', force = TRUE,
                           upgrade = "always")
}



library(tidyverse)
library(stringr)
library(lubridate)
library(magrittr)
library(forcats)
library(STADEM)

# set up folder structure
stademFolder = 'STADEM_results'
if(!dir.exists(stademFolder)) {
  dir.create(stademFolder)
}

modelFolder = 'ModelFiles'
if(!dir.exists(modelFolder)) {
  dir.create(modelFolder)
}

## Requires a copy of the Lower Granite Dam Trap Database and the odbc driver
if(.Platform$OS.type != 'unix') {
  source('./R/loadLGTrappingDBase.R')
  trap_filepath <- './data/TrappingDBase/LGTrappingExportJodyW.accdb'
  con <- loadLGTrappingDBase(trap_filepath)
  trap_dbase <- DBI::dbReadTable(con, 'tblLGDMasterCombineExportJodyW')
}

# this works for Mac, or any system really
if(.Platform$OS.type == 'unix') {
  # trap_filepath <- './data/TrappingDBase/LGTrappingExportJodyW.accdb'
  trap_filepath <- './data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'
  trap_dbase = readLGRtrapDB(trap_filepath)
}

# save .csv of dbase for later use
write_csv(trap_dbase, file = './data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv')

# set species and spawn year
species = c('Chinook', 'Steelhead')  # either Chinook or Steelhead

# running only a single year to keep previous estimates the same
#year = 2010:2018        # tagging operations started at Lower Granite with spawn year 2009.
year = 2020

# Loop through species and years
for(i in 1:length(species)){
  for(j in 1:length(year)){
    
    spp <- species[i]
    yr <- year[j]
    
    if(spp == 'Chinook'){
      incl_jacks = TRUE
    } else{
      incl_jacks = FALSE
    }
    
    if(spp == 'Chinook'){
      start_date = paste0(yr,'0301')
      end_date = paste0(yr,'0817')
    } else {
      start_date = paste0(yr-1,'0701')
      end_date = paste0(yr,'0630')
    }
    
    
    stadem_list = compileGRAdata(spp = spp,
                                 start_date = start_date,
                                 end_date = end_date,
                                 strata_beg = 'Mon',
                                 trap_dbase = trap_dbase,
                                 incl_jacks = incl_jacks)   
    
    # creates JAGs data list.
    jags_data_list = prepJAGS(stadem_list[['weeklyData']])
    
    # JAGs needs to access a .txt file of the model code
    model_file_nm = './ModelFiles/STADEM_LGR_model.txt'
    
    # what distribution to use for window counts?
    win_model = c('pois', 'neg_bin', 'neg_bin2', 'quasi_pois', 'log_space')[2]
    
    # run model - runSTADEMmodel sets the params to save inside the fnc, and it
    # does not save all the params available in the model!
    stadem_mod = runSTADEMmodel(file_name = model_file_nm,
                                mcmc_chainLength = 40000,
                                mcmc_burn = 10000,
                                mcmc_thin = 30,
                                mcmc_chains = 4,
                                jags_data = jags_data_list,
                                seed = 5,
                                weekly_params = T,
                                win_model = win_model)
    
    #ests <- stadem_mod$summary
    
    # save results
    save(stadem_mod, stadem_list,
         file = paste0(stademFolder,'/LGR_STADEM_', spp, '_', yr, '.rda'))
    
  } # close j loop
} # close i loop
