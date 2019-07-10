#------------------------------------------------------------------------------
# Script loads DABOM and STADEM model results and combines posteriors to estimate
# abundance at each DABOM site.  Population abundance posteriors is then combined
# with sex and age proporiton abundances to estimate the abundance of each
# life history group.
#
# Author: Ryan Kinzer
#------------------------------------------------------------------------------
# Load packages
#------------------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(PITcleanr)
library(DABOM)
library(WriteXLS)
source('./R/definePopulations.R')
source('./R/summarisePosterior.R')

# set up folder structure
AbundanceFolder = 'Abundance_results' # for processed files
if(!dir.exists(AbundanceFolder)) {
  dir.create(AbundanceFolder)
}

dabom_files <- list.files('./DABOM_results')
stadem_files <- list.files('./STADEM_results')

#------------------------------------------------------------------------------
# set species, spawn year and time stamp
spp <- 'Steelhead'  # either Chinook or Steelhead
yr <- 2018        # tagging operations started at Lower Granite with spawn year 2009.
timestp <- gsub('[^0-9]','', Sys.Date())

#for(yr in 2016:2018){
  
  load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',yr,'.rda')) 
  
# steps are necessary b/c each DABOM filename contains a time stamp
  spp_files <- dabom_files[str_detect(dabom_files, spp)]
  tmp_file <- spp_files[str_detect(spp_files, as.character(paste0("_",yr,"_")))]
  load(file = paste0('./DABOM_results/',tmp_file))
  
  load(paste0('./Sex_results/Population_SexProp_',spp,'_',yr,'.rda'))
  #load(paste0('./Age_results/Population_AgeProp_',spp,'_',yr,'.rda'))

#------------------------------------------------------------------------------
## Gather Detection Probabilities
#------------------------------------------------------------------------------
detect_summ = summariseDetectProbs(dabom_mod = dabom_mod,
                                   capHist_proc = proc_list$proc_ch) %>%
  mutate(spawn_yr = yr,
         species = spp,
         cv = sd/median) %>%
  select(spawn_yr, species, Node, n_tags, estimate = median, sd, cv, lowerCI, upperCI)

#------------------------------------------------------------------------------
## Gather All Transition Probabilities for STADEM
#------------------------------------------------------------------------------
trans_summ <- summariseTransProbs_LGD(dabom_mod,
                                      cred_int_prob = 0.95) %>%
  mutate(spawn_yr = yr,
         species = spp) %>%
  select(spawn_yr, species, everything())

#------------------------------------------------------------------------------
## Gather Weekly Transition Probabilities 
#------------------------------------------------------------------------------
wk_trans_summ <- compileWeekTransProbs(dabom_mod) %>%
  mutate(spawn_yr = yr,
         species = spp) %>%
  group_by(spawn_yr, species, week, branch) %>%
  summarise(estimate = median(prob),
            sd = sd(prob),
            cv = sd/estimate,
            lowerCI = quantile(prob, probs = .025),
            upperCI = quantile(prob, probs = .975)) 

#------------------------------------------------------------------------------
## Gather Tributary Estimates 
#------------------------------------------------------------------------------
trib_summ = calcTribEscape_LGD(dabom_mod,
                               stadem_mod,
                               stadem_param_nm = 'X.new.wild',
                               bootstrap_samp = 2000, #2000
                               node_order = proc_list$NodeOrder,
                               summ_results = T,
                               pt_est_nm = 'median',
                               cred_int_prob = 0.95) %>%
  mutate(spawn_yr = yr,
         species = spp) %>%
  select(spawn_yr, species, everything())

#------------------------------------------------------------------------------
## Gather Population Group Estimates -- NEED TO WORK ON!!!!!
#------------------------------------------------------------------------------
pop_df <- definePopulations(spp)

#------------------------------------------------------------------------------
## Gather Posterior Estimates -- NEED TO WORK ON!!!!!
#------------------------------------------------------------------------------

# abundance at main branches, upstream sites and black-boxes
abundance_post = calcTribEscape_LGD(dabom_mod,
                                 stadem_mod,
                                 node_order = proc_list$NodeOrder,
                                 summ_results = F) %>%
  mutate(spawn_yr = yr, species = spp) %>%
  select(spawn_yr, species, everything())

# Create posteriors of population abundance by combining branches
N_post = pop_df %>%
  left_join(abundance_post, by = 'area') %>%
  group_by(TRT, iter) %>%
  summarise_at(vars(escape),
               funs(sum),
               na.rm = T) %>%
  ungroup() %>%
  mutate(spawn_yr = yr, species = spp) %>%
  rename(N = escape) %>%
  select(spawn_yr, species, everything())

sex_post <- sex_mod$sims.list$p %>%
  as_tibble() %>%
  gather(key = 'popNum', value = p) %>%
  mutate(popNum = as.integer(gsub('V','',popNum))) %>%
  group_by(popNum) %>%
  mutate(iter = 1:n()) %>%
  left_join(modSexDf %>%
              filter(!is.na(TRT)) %>%
              mutate(popNum = sex_jagsData$popNum),
            by = 'popNum')

N_pop <- N_post %>%
  group_by(TRT) %>%
  sample_n(size = 1000, replace = T) %>%
  mutate(iter = 1:n())

p_pop <- sex_post %>%
  group_by(TRT) %>%
  sample_n(1000, replace = T) %>%
  mutate(iter = 1:n())


pop_post <- left_join(N_pop, p_pop, by = c('TRT', 'iter')) %>%
  mutate(N_females = N * p)


N_pop_summ <- summarisePosterior(pop_post, N, TRT)
p_f <- summarisePosterior(pop_post, p, TRT)
N_f_summ <-summarisePosterior(pop_post, N_females, TRT)
  





# ggplot(pop_post, aes(x = N_females)) +
#   geom_density() +
#   facet_wrap(~TRT, scales = 'free')











save(detect_summ,trans_summ, wk_trans_summ, trib_summ, report_summ,
     file = paste0(AbundanceFolder, '/LGR_PIT_estimates_',timestp,'.rda'))

#------------------------------------------------------------------------------
# Load model run/estimates
load('./DABOM_estimates/LGR_PIT_estimates_20180314.rda')

library(WriteXLS)
#testPerl()

WriteXLS(x = c('detect_summ', 'trib_summ', 'report_summ'),
         ExcelFileName = paste0('./DABOM_estimates/LGR_PIT_estimates_',timestp,'.xlsx'),
         SheetNames = c('Detection_Eff','Site_ests','Population_ests'))

#}
