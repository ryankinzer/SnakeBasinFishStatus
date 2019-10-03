#------------------------------------------------------------------------------
# Script loads data processed through PIT-cleanr and checked by the
# biologist then runs the DABOM package and model.  The script then uses 
# built in DABOM functions to extract and save the results
#
# Author: Ryan Kinzer
#------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(PITcleanr)
library(DABOM)

# set up folder structure
dabomFolder = 'DABOM_results'
if(!dir.exists(dabomFolder)) {
  dir.create(dabomFolder)
}

#------------------------------------------------------------------------------
# Set species and year of interest
#------------------------------------------------------------------------------
spp <- 'Steelhead'
yr <- 2018 

timestp <- gsub('[^0-9]','', Sys.Date())

#------------------------------------------------------------------------------
# Load required DABOM data 
#------------------------------------------------------------------------------
load(paste0('data/DABOMready/LGR_',spp,'_',yr,'.rda'))

proc_ch <- proc_list$ProcCapHist %>%
  filter(UserProcStatus)

#------------------------------------------------------------------------------
# Remove all nodes within steelhead branches that we don't care about
# estimating for Chinook
#------------------------------------------------------------------------------

if(spp == 'Chinook') {
  # remove some Chinook detections from the data; INCREASES THE BLACK BOX!!!
  sthd_only_nodes = proc_list$NodeOrder %>%
    filter(Group %in% c('Asotin', 'Lapwai', 'Potlatch', 'JosephCreek', 'CowCreek', 'CarmenCreek', 'Almota', 'Alpowa', 'Penawawa')) %>%
    select(Node) %>%
    as.matrix() %>%
    as.character()
  
   proc_ch = proc_ch %>%
     filter(!Node %in% sthd_only_nodes)
   
 }

#------------------------------------------------------------------------------
# Because HLM was giving us problems...
# Switch Potlatch detections - move POTREF and POTRWF to HLMA0 with det = 1.0
#------------------------------------------------------------------------------
proc_ch <- proc_ch %>%
  mutate(Node = ifelse(Node %in% c('POTREF', 'POTRWF'), 'HLMA0', Node))

#------------------------------------------------------------------------------
# In 2010 and 2011, the Asotin Creek weir was upstream of ACB. For those years, 
# take all detections at ASOTIC and move them to ACBA0
#------------------------------------------------------------------------------
if(spp == 'Steelhead' & yr %in% c(2010, 2011)) {
  proc_ch <- proc_ch %>%
    mutate(Node = if_else(Node == 'ASOTIC',
                          'ACBA0',
                          Node))
}

#------------------------------------------------------------------------------
# Create default LGR branch occupancy JAGs model code.
#------------------------------------------------------------------------------

# should initial movement probabilities be time-varying?
time_varying = TRUE

# file path to the default and initial model 

basic_modNm = 'ModelFiles/LGR_DABOM.txt'

writeDABOM_LGD(file_name = basic_modNm,
               time_varying = time_varying)

#------------------------------------------------------------------------------
# Alter default model code for species and year of 
# interest; sets prior for some detection node efficiencies at 0 or 100%
# based on actual tag detection data; 0% if no tags were seen
#------------------------------------------------------------------------------

# filepath for specific JAGS model code for species and year
mod_path = paste0('ModelFiles/LGR_DABOM_', spp, '_', yr, '.txt')

#writes species and year specific jags code
fixNoFishNodes(basic_modNm,
               mod_path,
               proc_ch,
               proc_list$NodeOrder)

#------------------------------------------------------------------------------
# Create capture history matrices for each main branch to be used in 
# the JAGS data list
#------------------------------------------------------------------------------
dabom_list = createDABOMcapHist(proc_ch,
                                proc_list$NodeOrder,
                                split_matrices = T)
#------------------------------------------------------------------------------
# Only Used to Debug
#------------------------------------------------------------------------------
# full_dabom = createDABOMcapHist(proc_ch,
#                                proc_list$NodeOrder,
#                                split_matrices = F)
#------------------------------------------------------------------------------

# Creates a function to spit out initial values for MCMC chains
init_fnc = setInitialValues_LGD(dabom_list)

#Create all the input data for the JAGS model
jags_data = createJAGSinputs_LGD(dabom_list)


# CHANGE TIME VARY DATE FOR SPECIES!!!!!!

if(time_varying) {
  jags_data = c(jags_data,
                addTimeVaryData(proc_ch,
                                node_order = proc_list$NodeOrder,
                                start_date = paste0(yr-1,'0701'), 
                                end_date = paste0(yr,'0630')))
}

#------------------------------------------------------------------------------
# Tell JAGS which parameters in the model that it should save.
# the fnc is hard coded and needs to be updated if there are changes!
#------------------------------------------------------------------------------

jags_params = setSavedParams_LGD(time_varying = time_varying)

#------------------------------------------------------------------------------
# Run the model

# Recommended MCMC parameters are:
#   
#   * `n.chains`: 4
# * `n.iter`: 5,000
# * `n.burnin`: 2,500
# * `n.thin`: 10
# 4*(5000+2500) = 30000
# 1 iteration takes about .18 minutes
# about 3.75 days!!!!!

library(jagsUI)
set.seed(12)
dabom_mod <- jags.basic(data = jags_data,
                        inits = init_fnc,
                        parameters.to.save = jags_params,
                        model.file = mod_path,
                        n.chains = 4,
                        n.iter = 5000,
                        n.burnin = 2500,
                        n.thin = 10,
                        DIC = T)



#--------------------------------------------------------------------------------
# Save the results
#--------------------------------------------------------------------------------
proc_list[["proc_ch"]] <- proc_ch

save(dabom_mod, dabom_list, proc_list,
     file = paste0(dabomFolder,'/LGR_DABOM_Bio_', spp, '_', yr,'_',timestp,'.rda'))

detach('package:jagsUI')