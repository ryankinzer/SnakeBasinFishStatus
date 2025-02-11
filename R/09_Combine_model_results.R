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
source('./R/definePopulations2.R')
source('./R/summarisePosterior.R')
source('./R/assign_POP_GSI.R')

# load configuration file
load('./data/ConfigurationFiles/site_config.rda')

# set up folder structure
AbundanceFolder = 'Abundance_results' # for processed files
if(!dir.exists(AbundanceFolder)) {
  dir.create(AbundanceFolder)
}

#dabom_files <- list.files('./DABOM_results')
#stadem_files <- list.files('./STADEM_results')

#------------------------------------------------------------------------------
# set species, spawn year and time stamp
species <- 'Steelhead' #c('Steelhead', 'Chinook')  # either Chinook or Steelhead
yr <- 2020        # tagging operations started at Lower Granite with spawn year 2009.
#timestp <- gsub('[^0-9]','', Sys.Date())

for(spp in species){
  
  if(spp == 'Steelhead'){
    year_range <- c(2010:2019)
  } else {
    year_range <- c(2010:2019)
  }

for(yr in year_range){
  
  load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',yr,'.rda')) 
  load(paste0('./DABOM_results/LGR_DABOM_',spp,'_',yr,'.rda'))
  # steps are necessary b/c each DABOM filename contains a time stamp
  #spp_files <- dabom_files[str_detect(dabom_files, spp)]
  #tmp_file <- spp_files[str_detect(spp_files, as.character(paste0("_",yr,"_")))]
  
  
  load(paste0('./Sex_results/Population_SexProp_',spp,'_',yr,'.rda'))
  load(paste0('./Age_results/Population_AgeProp_',spp,'_',yr,'.rda'))
  
  #------------------------------------------------------------------------------
  ## summarise tag information
  #------------------------------------------------------------------------------
  # get which population each site is assigned to
  site_pop = assign_POP_GSI(species = spp, 
                            configuration = configuration,
                            site_df = site_df)[[1]] %>%
    as_tibble() %>%
    select(SiteID:SiteName, Node, MPG:GSI_Group)
  
  
  tag_summ = summariseTagData(capHist_proc = proc_list$proc_ch,
                              trap_data = proc_list$ValidTrapData) %>%
    left_join(site_pop %>%
                select(Node, MPG, POP_NAME, TRT),
              by = c('AssignSpawnNode' = 'Node'))
  
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
  # how many samples from each posterior?
  nSamps = 1000
  
  # abundance at main branches, upstream sites and black-boxes
  abundance_post = calcTribEscape_LGD(dabom_mod,
                                      stadem_mod,
                                      node_order = proc_list$NodeOrder,
                                      summ_results = F) %>%
    mutate(spawn_yr = yr, species = spp) %>%
    select(spawn_yr, species, everything())
  
  # Create posteriors of population abundance by combining branches
  N_post = pop_df %>%
    left_join(abundance_post, by = 'site') %>%
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
  
  poss_ages <- modAgeDf %>%
    filter(!is.na(TRT)) %>%
    select(starts_with('age'))

  poss_ages <- as_tibble(na.omit(as.numeric(
    gsub('[[:alpha:]]','', names(poss_ages[,!colSums(poss_ages) == 0])))))  %>%
    mutate(age_fct = as.integer(as.factor(value))) %>%
    rename(age = value)

  age_post <- age_mod$sims.list$pi %>%
    as.data.frame.table() %>%
    as_tibble() %>%
    rename(iter = Var1,
           popNum = Var2,
           age_fct = Var3,
           ageProp = Freq) %>%
    mutate_at(vars(iter, popNum, age_fct),
              list(as.integer)) %>%
    left_join(poss_ages, by = 'age_fct') %>%
    select(-age_fct) %>%
    #mutate(age = age + 1) %>%
    left_join(modAgeDf %>%
                filter(!is.na(TRT)) %>%
                mutate(popNum = age_jagsData$popNum) %>%
                select(popNum, MPG:TRT) %>%
                distinct(),
              by = 'popNum')
  
  N_pop <- N_post %>%
    group_by(TRT) %>%
    sample_n(size = nSamps, replace = T) %>%
    mutate(iter = 1:n())
  
  p_pop <- sex_post %>%
    group_by(TRT) %>%
    sample_n(nSamps, replace = T) %>%
    mutate(iter = 1:n()) %>%
    select(iter, TRT, p)
  
  age_pop <- age_post %>%
    mutate(age = paste0('age', age)) %>%
    spread(age, ageProp) %>%
    group_by(TRT) %>%
    sample_n(nSamps, replace = T) %>%
    mutate(iter = 1:n()) %>%
    select(iter, TRT, starts_with('age'))
  
  pop_post <- N_pop %>%
    left_join(p_pop, 
              by = c('TRT', 'iter')) %>%
    left_join(age_pop,
              by = c('TRT', 'iter')) %>%
    mutate(N_females = N * p) %>%
    mutate_at(vars(starts_with('age')),
              list(~ . * N))
  
  # summarize total escapement by population
  N_pop_summ <- summarisePosterior(pop_post, N, TRT) %>%
    mutate(spawn_yr = yr, 
           species = spp) %>%
    left_join(tag_summ %>%
                group_by(TRT) %>%
                summarise(n_tags = n_distinct(TagID))) %>%
    select(spawn_yr, species, TRT, n_tags, everything())
  # summarize proportion female
  p_f <- summarisePosterior(pop_post, p, TRT, round = F) %>%
    mutate(spawn_yr = yr, 
           species = spp) %>%
    select(spawn_yr, species, everything())
  # summarize number of females
  N_f_summ <- summarisePosterior(pop_post, N_females, TRT) %>%
    mutate(spawn_yr = yr, 
           species = spp) %>%
    select(spawn_yr, species, everything())
  
  # summarize the age proportions
  
  age_prop_summ <- age_post %>%
    group_by(age) %>%
    nest() %>%
    mutate(sumPost = map(data,
                         ~summarisePosterior(data = .x, ageProp, TRT, round = F)
    )
    ) %>%
    select(age, sumPost) %>%
    unnest(sumPost) %>%
    mutate(spawn_yr = yr,
           species = spp) %>%
    select(spawn_yr, species, age, everything())
  
  # age_prop_summ = summarisePosterior(age_pop, age2, TRT, round = F) %>%
  #   mutate(age = 2) %>%
  #   bind_rows(summarisePosterior(age_pop, age3, TRT, round = F) %>%
  #               mutate(age = 3)) %>%
  #   bind_rows(summarisePosterior(age_pop, age4, TRT, round = F) %>%
  #               mutate(age = 4)) %>%
  #   bind_rows(summarisePosterior(age_pop, age5, TRT, round = F) %>%
  #               mutate(age = 5)) %>%
  #   bind_rows(summarisePosterior(age_pop, age6, TRT, round = F) %>%
  #               mutate(age = 6)) %>%
  #   mutate(spawn_yr = yr, 
  #          species = spp) %>%
  #   select(spawn_yr, species, age, everything())
  
  # summarize the escapement by age
  age_summ <- pop_post %>%
    select(-N_females) %>%
    gather(age, N_age, starts_with('age')) %>%
    mutate(age = gsub('age','',age)) %>%
    group_by(age) %>%
    nest() %>%
    mutate(sumPost = map(data,
                         ~summarisePosterior(data = .x, N_age, TRT, round = F)
    )
    ) %>%
    select(age, sumPost) %>%
    unnest(sumPost) %>%
    mutate(spawn_yr = yr,
           species = spp) %>%
    select(spawn_yr, species, TRT, age, everything())
  # age_summ <- summarisePosterior(pop_post, age2, TRT) %>%
  #   mutate(age = 2) %>%
  #   bind_rows(summarisePosterior(pop_post, age3, TRT) %>%
  #               mutate(age = 3)) %>%
  #   bind_rows(summarisePosterior(pop_post, age4, TRT) %>%
  #               mutate(age = 4)) %>%
  #   bind_rows(summarisePosterior(pop_post, age5, TRT) %>%
  #               mutate(age = 5)) %>%
  #   bind_rows(summarisePosterior(pop_post, age6, TRT) %>%
  #               mutate(age = 6)) %>%
  #   mutate(spawn_yr = yr, 
  #          species = spp) %>%
  #   select(spawn_yr, species, age, everything())

  # save posteriors
  save(N_pop, p_pop, age_pop,
       file = paste0(AbundanceFolder, '/LGR_Posteriors_', spp, '_', yr, '.rda'))
  
  # save summaries
  save(detect_summ, trans_summ, wk_trans_summ, trib_summ, N_pop_summ, p_f, N_f_summ, age_prop_summ, age_summ,
       file = paste0(AbundanceFolder, '/LGR_Summary_', spp, '_', yr, '.rda'))

  }
}

#------------------------------------------------------------------------------
# combine summaries by year
#------------------------------------------------------------------------------
valid_est <- read_csv('./data/ConfigurationFiles/valid_trt_estimates.csv')

for(spp in species){
  
  if(spp == 'Steelhead'){
    year_range <- c(2010:2020)
    prod_range <- 2010:2015
  } else {
    year_range <- c(2010:2019)
    prod_range <- 2010:2014
  }
  
  allSumm = as.list(year_range) %>%
    rlang::set_names() %>%
    map(.f = function(x) {
      load(paste0(AbundanceFolder, '/LGR_Summary_', spp, '_', x[1], '.rda'))
      list(detect_summ = detect_summ,
           trib_summ = trib_summ,
           N_pop_summ = N_pop_summ,
           p_f = p_f, 
           N_f_summ = N_f_summ, 
           age_prop_summ = age_prop_summ, 
           age_summ = age_summ)
    })
  
  
  detect_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'detect_summ')
  
  trib_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'trib_summ')
  
  N_pop_all = inner_join(valid_est,
                        allSumm %>%
                          map_df(.id = NULL,
                          .f = 'N_pop_summ'))
  
  p_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'p_f')
  
  N_f_all = inner_join(valid_est,
                      allSumm %>%
                        map_df(.id = NULL,
                        .f = 'N_f_summ'))
  
  age_prop_all = allSumm %>%
    map_df(.id = NULL,
           .f = 'age_prop_summ')
  
  age_all = inner_join(valid_est,
                      allSumm %>%
                        map_df(.id = NULL,
                        .f = 'age_summ') %>%
                        mutate(brood_yr = spawn_yr - as.numeric(age)) %>%
                        select(spawn_yr, species, brood_yr, everything()))
  
  #------------------------------------------------------------------------------
  # Build brood year tables
  # Still need to figure out how to best filter out years where we don't have complete brood year returns
  
  # combine posteriors by year
  allBrYr = as.list(year_range) %>%
    rlang::set_names() %>%
    map_df(.f = function(x) {
      load(paste0(AbundanceFolder, '/LGR_Posteriors_', spp, '_', x[1], '.rda'))
      N_pop %>%
        left_join(age_pop %>%
                    gather(age, ageProp, starts_with('age'))) %>%
        filter(!is.na(age)) %>%
        mutate(brood_yr = spawn_yr - as.integer(str_remove(age, 'age')),
               Nage = N * ageProp) %>%
        #group_by(spawn_yr, TRT) %>%
        #select(spawn_yr, species, TRT, iter, brood_yr, N, Nage) %>%
        ungroup()
    })
  
  allBrYr %>%
    group_by(species, TRT, spawn_yr) %>%
    summarise(nBrYrs = n_distinct(brood_yr),
              minBrYr = min(brood_yr),
              maxBrYr = max(brood_yr))
  
  allBrYr %>%
    group_by(species, TRT, brood_yr) %>%
    summarise(nSpYrs = n_distinct(spawn_yr),
              minSpYr = min(spawn_yr),
              maxSpYr = max(spawn_yr))
  
  spawnRec_post = allBrYr %>%
    group_by(iter, species, TRT, brood_yr) %>%
    summarise_at(vars(R = Nage),
                 list(sum)) %>%
    left_join(allBrYr %>%
                select(spawn_yr:iter, N) %>%
                distinct() %>%
                rename(brood_yr = spawn_yr,
                       S = N)) %>%
    select(iter:brood_yr, S, R) %>%
    mutate(lambda = R / S) %>%
    ungroup()
  
  spawnRec_summ = spawnRec_post %>%
    gather(var, value, S:lambda) %>%
    group_by(species, TRT, brood_yr, var) %>%
    summarise_at(vars(value),
                 list(mean = mean, 
                      median = median,
                      # mode = estMode,
                      SE = sd),
                 na.rm = T) %>%
    select(-median, -SE) %>%
    spread(var, mean) %>%
    select(species, brood_yr, TRT, S, R, lambda) %>%
    filter(!is.na(lambda))
  
  #-----------------------------------------------------------------------------
  # RK test
  brood_table <- allBrYr %>%
    select(spawn_yr:iter, N) %>%
    distinct() %>%
    group_by(spawn_yr, species, TRT) %>%
    summarise(S = median(N, na.rm=TRUE)) %>%
    full_join(allBrYr %>%
                group_by(brood_yr, species, TRT, age) %>%
                summarise(Nage = median(Nage, na.rm=TRUE)) %>%
                spread(age, Nage), 
              by = c('spawn_yr' = 'brood_yr', 'species', 'TRT')) %>%
    rename(brood_yr = spawn_yr) %>%
    arrange(species, TRT, brood_yr) %>%
    ungroup() %>%
    mutate(nRyrs = rowSums(!is.na(select(., contains('age'))))) %>%
    mutate(R = rowSums(select(.,contains('age')),na.rm = TRUE),
           lambda = ifelse(!is.na(S) & R != 0,R/S,""))
  
  prod_df = as.list(prod_range) %>%
    rlang::set_names() %>%
    map_df(.f = function(x) {
      
      tmp <- spawnRec_post %>% filter(species == spp,
                                      brood_yr == x[1])
      
      
      Sdf <- summarisePosterior(tmp, S, TRT, round = T) %>%
        mutate(brood_yr = x[1],
               species = spp,
               variable = 'Spawners')
      Rdf <- summarisePosterior(tmp, R, TRT, round = T) %>%
        mutate(brood_yr = x[1],
               species = spp,
               variable = 'Recruits')
      Ldf <- summarisePosterior(tmp, lambda, TRT, round = F) %>%
        mutate(brood_yr = x[1],
               species = spp,
               variable = 'lambda')
      
      bind_rows(Sdf, bind_rows(Rdf, Ldf)) %>% select(brood_yr, species, TRT, variable, everything()) %>% arrange(species, TRT, brood_yr)
      
    })
  
  # STADEM Estimates
  # Lower Granite Estimates
  
  allLGR <- as.list(year_range) %>%
    rlang::set_names() %>%
    map_df(.id = 'spawn_year',
           .f = function(x){
             
          load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',x[1],'.rda'))
             
          stadem_mod$summary %>% as_tibble(rownames = 'param') %>%
               filter(grepl('X.tot.new.',param)) %>%
               mutate(species = spp,
                      origin = case_when(grepl('all', param) ~ 'Total',
                                         grepl('wild',param) ~ 'Natural',
                                         grepl('hatch', param) ~ 'Hatchery Clipped',
                                         grepl('hnc', param) ~ 'Hatchery No-Clipped')) %>%
               select(species,
                      origin,
                      estimate = `50%`,
                      lowerCI = `2.5%`,
                      upperCI = `97.5%`,
                      mean,
                      sd)
             
           })
  
  allLGRtags <- as.list(year_range) %>%
    rlang::set_names() %>%
    map_df(.id = 'spawn_year',
           .f = function(x){
             load(paste0('./DABOM_results/LGR_DABOM_',spp,'_',x[1],'.rda'))
             
             tag_dat <- proc_list$ValidTrapData %>%
               summarise(valid_tags = n_distinct(LGDNumPIT)) %>% 
               mutate(origin = 'Natural')
           })
  
  stadem_df <- left_join(allLGR, allLGRtags, by = c('spawn_year', 'origin')) %>%
    select(spawn_year, species, origin, valid_tags, everything())
  
  
  # Get unique tags per site......should be moved to DABOM tribCalcEstimate fnc.
  site_tags <- as.list(year_range) %>%
    rlang::set_names() %>%
    map_df(.id = 'spawn_year',
           .f = function(x){
             load(paste0('./DABOM_results/LGR_DABOM_',spp,'_',x[1],'.rda'))
             
             proc_list$proc_ch %>%
               filter(UserProcStatus) %>%
               group_by(SiteID) %>%
               summarise(n_tags = n_distinct(TagID)) %>%
               mutate(species = spp) %>%
               select(species, everything())
             
           })
  
  #Save all data as .xlsx
  list('LGR Esc' = stadem_df, 
       'Pop Total Esc' = N_pop_all,
       'Pop Female Esc' = N_f_all,
       'Pop Age Esc' = age_all,
       'Pop Brood Table' = brood_table,
       'Pop Stock Recruit' = prod_df,
       'Pop Female Props' = p_all,
       'Pop Age Props' = age_prop_all,
       'Site Esc' = trib_all,
       'Site Unique Tags' = site_tags,
       'Node Detect Eff' = detect_all) %>%
   writexl::write_xlsx(paste0(AbundanceFolder, '/LGR_AllSummaries_',spp,'.xlsx'))
  
}
