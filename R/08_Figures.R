# Purpose: Create model diagnostic and result figures.
# Authors: Kevin See and Ryan N. Kinzer
# Created: 7/12/19
# Modified: 

spp = 'Steelhead'

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
