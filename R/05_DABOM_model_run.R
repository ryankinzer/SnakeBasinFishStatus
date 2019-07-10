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

#------------------------------------------------------------------------------
# Set species and year of interest
#------------------------------------------------------------------------------
spp <- 'Chinook'
yr <- 2018 

timestp <- gsub('[^0-9]','', Sys.Date())

#------------------------------------------------------------------------------
# Load configuration, parent_child and processed datasets from PITcleanr 
#------------------------------------------------------------------------------
#load(paste0('data/PreppedData/LGR_', spp, '_', yr,'_20190125.rda'))
load(paste0('data/PreppedData/LGR_Chinook_2018_20190304.rda'))

#proc_ch <- proc_list$ProcCapHist %>%
#  mutate(UserProcStatus = AutoProcStatus)
#------------------------------------------------------------------------------
# Load biologist corrected capture history file with the final TRUE/FALSE calls
# in the UserProcStatus column, and then filter out all FALSE observations.
#------------------------------------------------------------------------------
#proc_ch <- read_delim(paste0('data/CleanedData/ProcCapHist_EDITTED_', spp, '_', yr,'_20190125.txt'), header = TRUE, sep = '\t') %>%
proc_ch <- read_delim('data/CleanedData/ProcCapHist_EDITTED_Steelhead_2018_20190125_RO_2-29-19.txt', delim = '\t') %>%
    mutate(TrapDate = mdy(TrapDate),
         ObsDate = mdy_hms(ObsDate),
         lastObsDate = mdy_hms(lastObsDate),
         AutoProcStatus = ifelse(AutoProcStatus == 1, TRUE, FALSE),
         UserProcStatus = ifelse(UserProcStatus == 1, TRUE, FALSE),
         ModelObs = ifelse(ModelObs == 1, TRUE, FALSE)) %>%
  filter(ModelObs)

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
# Create default LGR branch occupancy JAGs model code.
#------------------------------------------------------------------------------

# should initial movement probabilities be time-varying?
time_varying = TRUE

# file path to the default and initial model 
#  (works for both Chinook and steelhead and all years)

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
# Used to Debug
#------------------------------------------------------------------------------
full_dabom = createDABOMcapHist(proc_ch,
                               proc_list$NodeOrder,
                               split_matrices = F)
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
                                spawn_yr = yr,
                                spp = spp,
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
     file = paste0('DABOM_results/LGR_DABOM_Bio_', spp, '_', yr,'_',timestp,'.rda'))

detach('package:jagsUI')
#------------------------------------------------------------------------------
# Summarize the results - everything in this section should be moved 
# to another script but we can leave it here for now!!!
#------------------------------------------------------------------------------

load(file = paste0('DABOM_results/LGR_DABOM_Bio_', spp, '_', yr, '_20190304.rda'))

# Detection Probabilities Directly from `DABOM`, we can extract estimates of
# the detection probability of the observation nodes. These are average
# probabiliities across the entire season of the run, and how these nodes match
# up to actual arrays and antennas is defined in the configuration file.

detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                   capHist_proc = proc_list$proc_ch)
head(detect_summ)

# Detections Probs From PITcleanr
config_filepath <- './data/ConfigurationFiles/my_config_20190304.csv'
my_config <- read_csv(config_filepath)

sitedf_filepath <- './data/ConfigurationFiles/site_df_20190304.csv'
site_df <- read_csv(sitedf_filepath, col_types = cols(.default = 'c'))

parentchild_filepath <- './data/ConfigurationFiles/parent_child_20190304.csv'
parent_child <- read_csv(parentchild_filepath)

valid_paths <- getValidPaths(parent_child, 'GRA')
node_order <- createNodeOrder(valid_paths, my_config, site_df, step_num = 3)

eff <- estNodeEff(proc_list$proc_ch, node_order, method = 'Chapman')

eff_df <- left_join(detect_summ, eff, by = 'Node')


eff_df %>%
  filter(Node != 'GRA') %>%
ggplot(aes(x = median, y = detEff, size = tagsAboveNode, size = n_tags)) +
  geom_point()+
  geom_text(aes(label = Node),size = 2, hjust = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme_bw()

eff_df %>%
  filter(Node != 'GRA') %>%
ggplot(aes(x= n_tags, y = estTagsAtNode, size = n_tags)) +
  geom_point()+
  geom_text(aes(label = Node),size = 2, hjust = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  theme_bw()

# Estimate Escapement with **STADEM** and **DABOM**

load(paste0('STADEM_results/LGR_STADEM_', spp, '_', yr, '.rda'))

# compiles posteriors of escapement
trib_summ = calcTribEscape_LGD(dabom_mod,
                               stadem_mod,
                               stadem_param_nm = 'X.new.wild',
                               bootstrap_samp = 5, #2000
                               node_order = proc_list$NodeOrder,
                               summ_results = T,
                               pt_est_nm = 'median',
                               cred_int_prob = 0.95)

report_summ = calcRepGrpEscape_LGD(dabom_mod = dabom_mod,
                                   stadem_mod = stadem_mod,
                                   node_order = proc_list$NodeOrder,
                                   pt_est_nm = 'median',
                                   spp = spp)

#------------------------------------------------------------------------------
## Save Results to an excel file - CURRENTLY NOT WORKING!!!
#------------------------------------------------------------------------------
list('Report Groups' = report_summ,
     'All Escapement' = trib_summ,
     'Detection' = detect_summ) %>%
  WriteXLS(ExcelFileName = paste0('DABOM_results/Escapement_', spp, '_', yr, '.xlsx'),
           AdjWidth = T,
           AutoFilter = T,
           BoldHeaderRow = T,
           FreezeRow = 1)