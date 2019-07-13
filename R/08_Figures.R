# Purpose: Create model diagnostic and result figures.
# Authors: Kevin See and Ryan N. Kinzer
# Created: 7/12/19
# Modified: 

library(tidyverse)
library(lubridate)
library(PITcleanr)
source('./R/site_trt_designations.R')

spp = 'Steelhead'
yr_range = 2016:2018

# Need GSI to PIT Comparisons!!!

# Detection ------------------------------------------------------------------------------
# make an observed vs predicted of detection probs.
allEffDf = as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x) {
           load(paste0('data/DABOMready/LGR_',spp,'_',x[1],'.rda'))
            
           eff <- estNodeEff(filter(proc_list$ProcCapHist,UserProcStatus), node_order = proc_list$NodeOrder)  
           
           trt_df <- site_trt_designations(spp, configuration)
           
           load(paste0('Abundance_results/LGR_Summary_', spp, '_', x[1], '.rda'))
           
           detect_summ %>%
             left_join(eff, by = 'Node') %>%
             left_join(proc_list$NodeOrder %>%
                         select(Node, Group), by = 'Node') %>%
             left_join(trt_df %>%
                         select(Node, MPG, TRT) %>% distinct(), by = 'Node')
         })


ObsPred_detect_p <- allEffDf %>%
  filter(!detEff == 1 & !estimate == 1, !detEff == 0 & !estimate == 0) %>%
  ggplot(aes(x = detEff, y = estimate, size = n_tags)) +
  geom_point()+
  geom_text(aes(label = Node),size = 2, hjust = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap(~spawn_year) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Predicted with PITcleanR',
       y = 'Predicted with DABOM',
       size = 'Unique Tags Observed')

ggsave('Figures/ObsVsPred_DetectProp.pdf',
       ObsPred_detect_p,
       width = 8,
       height = 8)

detect_p <- allEffDf %>%
  filter(estimate != 0 & estimate != 1) %>%
  filter(!is.na(Group)) %>%
  ggplot(aes(x = Node,
             y = estimate,
             color = Group,
             group = spawn_year)) +
  geom_linerange(aes(ymin = lowerCI,
                     ymax = upperCI),
                 position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_colour_viridis_d() +
  coord_flip() +
  facet_wrap(~spawn_year, scales = 'free', ncol = 3) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Node',
       y = 'Detection Probability',
       colour = 'Main Branch')

detect_p

ggsave('Figures/DetectProp.pdf',
       detect_p,
       width = 8,
       height = 10)

# Abundance----
# plot population abundance
allAbundDf = as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x) {
           
           load(paste0('data/DABOMready/LGR_',spp,'_',x[1],'.rda'))
           load(paste0('Abundance_results/LGR_Summary_', spp, '_', x[1], '.rda'))
           
           trt_df <- site_trt_designations(spp, configuration)
           
           N_pop_summ %>%
           left_join(trt_df %>%
                       select(MPG, TRT) %>% distinct(), by = 'TRT')
         })

pop_N <- allAbundDf %>%
  ggplot(aes(x = TRT,
             y = median,
             color = spawn_year,
             group = spawn_year)) +
  geom_linerange(aes(ymin = lowerCI,
                     ymax = upperCI),
                 position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_colour_viridis_d() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'TRT',
       y = 'Abundance',
       colour = 'Spawn Year')

pop_N

# Sex ------------------------------------------------------------------------------
# make an observed vs predicted female proportion plot, by population
allSexDf = as.list(2010:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x) {
           load(paste0('Sex_results/Population_SexProp_', spp, '_', x[1], '.rda'))
           
           sex_mod$summary %>%
             as_tibble(rownames = 'param') %>%
             filter(grepl('^p', param)) %>%
             mutate(popNum = str_extract(param, '[:digit:]+'),
                    popNum = as.integer(popNum)) %>%
             left_join(modSexDf %>%
                         filter(!is.na(TRT)) %>%
                         mutate(popNum = as.integer(as.factor(TRT))),
                       by = 'popNum') %>%
             bind_rows(sex_mod$summary %>%
                         as_tibble(rownames = 'param') %>%
                         filter(param == 'mu_ilogit')) %>%
             mutate_at(vars(MPG, TRT, modBranch),
                       list(fct_explicit_na),
                       na_level = 'Overall') %>%
             mutate(TRT = fct_reorder(TRT, mean),
                    TRT = fct_relevel(TRT, 'Overall',
                                      after = 0))
         })

sex_ObsVsPred_p = allSexDf %>%
  filter(!is.na(propF)) %>%
  ggplot(aes(x = propF,
             y = mean,
             color = nSexed)) +
  geom_abline(linetype = 2) +
  geom_point(aes(size = nSexed)) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  scale_color_viridis_c() +
  facet_wrap(~ TRT,
             scales = 'free') +
  labs(x = 'Observed',
       y = 'Predicted',
       title = 'Female Proportions',
       color = '# Sexed',
       size = '')

sex_ObsVsPred_p

ggsave('Figures/ObsVsPred_SexProp.pdf',
       sex_ObsVsPred_p,
       width = 8,
       height = 8)

sex_Pred <- allSexDf %>%
  #filter(spawn_year == '2018') %>%
  ggplot(aes(x = TRT,
             y = mean,
             color = spawn_year,
             group = spawn_year)) +
  geom_linerange(aes(ymin = `2.5%`,
                     ymax = `97.5%`),
                 position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0.5,
             linetype = 2) +
  ylim(0,1) +
  scale_color_viridis_d() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = 'bottom') +
  theme(strip.text.y = element_text(angle = .45)) +
  facet_grid(MPG~., scales = "free_y", space = 'free_y', labeller = label_wrap_gen(width = 10)) +
  labs(x = '',
       y = 'Female Proportion',
       title = paste(spp, 'Sex Proportions'),
       color = 'Spawn Year')

sex_Pred

ggsave('Figures/SexProp.pdf',
       sex_Pred,
       width = 8,
       height = 8)

# Age ------------------------------------------------------------------------------
# make an observed vs predicted age proportion plot, by population
allAgeDf = as.list(2010:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           load(paste0('Age_results/Population_AgeProp_', spp, '_', x[1], '.rda'))
           
           age_mod$summary %>%
             as_tibble(rownames = 'param') %>%
             filter(grepl('^pi', param)) %>%
             mutate(popNum = str_extract(param, '[:digit:]+'),
                    age = str_split(param, ',', simplify = T)[,2],
                    age = str_remove(age, '\\]')) %>%
             mutate_at(vars(popNum, age),
                       list(as.integer)) %>%
             mutate(age = age + 1) %>%
             left_join(modAgeDf %>%
                         filter(!is.na(TRT)) %>%
                         mutate(popNum = as.integer(as.factor(TRT))) %>%
                         mutate_at(vars(starts_with('age')),
                                   list(~ . / nAged)) %>%
                         gather(age, obsProp, starts_with('age')) %>%
                         mutate(age = str_remove(age, 'age'),
                                age = as.integer(age)))
         })

age_ObsVsPred_p = allAgeDf %>%
  filter(obsProp > 0) %>%
  ggplot(aes(x = obsProp,
             y = mean,
             color = as.factor(age))) +
  geom_abline(linetype = 2) +
  geom_point(aes(size = nAged)) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  scale_color_brewer(palette = 'Set1') +
  facet_wrap(~ TRT,
             scales = 'free') +
  labs(x = 'Observed',
       y = 'Predicted',
       title = 'Age Proportions',
       color = 'Age',
       size = '# Aged')

age_ObsVsPred_p

ggsave('Figures/ObsVsPred_AgeProp.pdf',
       age_ObsVsPred_p,
       width = 8,
       height = 8)

# plot showing estimated age proportions (mean) for each TRT population, by year
ageProp_p = allAgeDf %>%
  select(Year, MPG, TRT, age, mean) %>%
  distinct() %>%
  ggplot(aes(x = fct_rev(TRT),
             y = mean,
             fill = as.factor(age))) +
  geom_bar(stat = 'identity',
           position = position_stack(reverse = T)) +
  scale_fill_brewer(palette = 'Set1') +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 5),
        legend.position = 'bottom') +
  facet_wrap(~ Year) +
  labs(fill = 'Age',
       x = 'TRT',
       y = 'Proportion')

ggsave('Figures/AgePropEst.pdf',
       ageProp_p,
       width = 8,
       height = 8)


muVecDf = as.list(2010:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           load(paste0('Age_results/Population_AgeProp_', spp, '_', x[1], '.rda'))
           
           age_mod$summary %>%
             as_tibble(rownames = 'param') %>%
             filter(grepl('^avgPi', param)) %>%
             mutate(runType = str_extract(param, '[:digit:]'),
                    runType = if_else(runType == 1, 'A', 
                                      if_else(runType == 2, 
                                              'B', as.character(NA))),
                    age = str_sub(param, -2, -2),
                    age = as.integer(age),
                    age = age + 1) %>%
             mutate_at(vars(runType),
                       list(as.factor)) %>%
             mutate_at(vars(age),
                       list(~as.factor(as.character(.)))) %>%
             mutate(TRT = paste0(runType, '-run'),
                    MPG = paste0(runType, '-run'))
         })

muProp_p = muVecDf %>%
  select(Year, TRT, age, mean) %>%
  distinct() %>%
  ggplot(aes(x = fct_rev(TRT),
             y = mean,
             fill = as.factor(age))) +
  geom_bar(stat = 'identity',
           position = position_stack(reverse = T)) +
  scale_fill_brewer(palette = 'Set1') +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = 'right') +
  facet_wrap(~ Year) +
  labs(fill = 'Age',
       x = 'Run Type',
       y = 'Proportion')

ggsave('Figures/AvgAgePropEst_mu.pdf',
       muProp_p,
       width = 8,
       height = 8)
