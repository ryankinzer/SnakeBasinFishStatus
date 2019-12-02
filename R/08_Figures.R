# Purpose: Create model diagnostic and result figures.
# Authors: Kevin See and Ryan N. Kinzer
# Created: 7/12/19
# Modified: 

library(tidyverse)
library(lubridate)
library(PITcleanr)
library(scales)
source('./R/assign_POP_GSI.R')

spp = 'Steelhead'

if(spp == 'Steelhead'){
  yr_range = c(2010:2019)
} else {
  yr_range = c(2010:2018)
}

# Site/Management Designations
load('./data/ConfigurationFiles/site_config.rda')
pop_ls <- assign_POP_GSI(spp, configuration, site_df)
grp_df <- pop_ls[[1]] %>% sf::st_drop_geometry()
map_df <- pop_ls[[2]] %>% sf::st_drop_geometry()

# Lower Granite Estimates

allLGR <- as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x){
           
           load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',x[1],'.rda'))
         
           stadem_mod$summary %>% as_tibble(rownames = 'param')
           
           })

allLGRdata <- as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x){
           load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',x[1],'.rda'))
           
           stadem_list$weeklyData
         }
  )

allLGRtags <- as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x){
           
           #load(paste0('./STADEM_results/LGR_STADEM_',spp,'_',x[1],'.rda'))
           load(paste0('./DABOM_results/LGR_DABOM_',spp,'_',x[1],'.rda'))
           
           tag_dat <- proc_list$ValidTrapData
           
           if(spp == 'Steelhead'){
              strata <- STADEM::weeklyStrata(lubridate::ymd(paste0(x[1]-1,'0701')),lubridate::ymd(paste0(x[1],'0630')))
           } else {
              strata <- STADEM::weeklyStrata(lubridate::ymd(paste0(x[1],'0301')),lubridate::ymd(paste0(x[1],'0817')))
           }
           
           week = vector('integer', nrow(tag_dat))
           for(i in 1:length(strata)) {
             week[which(tag_dat$CollectionDate %within% strata[[i]])] = i
           }
           
           tag_dat %>%
             mutate(week = week) %>%
             arrange(CollectionDate)
         })

# Plot unique passage
lgr_new_fig <- allLGR %>%
  filter(grepl('X.tot.new', param)) %>%
  mutate(grp = case_when(
    grepl('all', param) ~ 'Total Passage',
    grepl('hatch', param) ~ 'Hatchery Ad-Clipped',
    grepl('hnc', param) ~ 'Hatchery No-Clipped',
    grepl('wild', param) ~ 'Natural Origin')
    ) %>%
  filter(grp != 'Total Passage') %>%
  ggplot(aes(x = as.integer(spawn_year), y = `50%`, group = grp)) +
  #geom_area(aes(fill = grp))
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .2) +
  geom_line(aes(colour = grp)) +
  geom_point(aes(colour = grp)) +
  facet_wrap(~ grp) +
  scale_color_viridis_d() +
  scale_y_continuous(labels = comma)+
  scale_x_continuous(breaks = pretty_breaks()) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Unique ', spp, ' Abundance at Lower Granite Dam'),
       x = 'Spawn Year',
       y = 'Abundance',
       colour = 'Origin Group')
lgr_new_fig

ggsave(paste0('Figures/LGR_unique_abund',spp,'.png'),
       lgr_new_fig,
       width = 7,
       height = 5)

# STADEM with Window
lgr_win_fig <- allLGR %>%
  filter(grepl('X.all\\[', param)) %>%
  mutate(week = as.integer(str_extract(param, '\\d+'))) %>%
  left_join(allLGRdata, by = c('spawn_year', 'week' = 'week_num')) %>%
  mutate(adj_win = win_cnt/(day_tags/tot_tags)) %>%
  mutate(adj_win = ifelse(adj_win == Inf, NA, adj_win)) %>%
  mutate(trap_est = ifelse(trap_valid,trap_est,NA)) %>%
 # filter(spawn_year == 2019) %>%
  ggplot(aes(x = Start_Date, group = spawn_year)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .2) +
  geom_line(aes(y = trap_est, colour = 'Trap')) +
  geom_line(aes(y = adj_win, colour = 'WindowAdj')) +
  geom_line(aes(y = win_cnt, colour = 'WindowRaw')) +
  geom_point(aes(y = `50%`, colour = 'STADEM')) +
  geom_point(aes(y = trap_est, colour = 'Trap')) +
  geom_point(aes(y = adj_win, colour = 'WindowAdj')) +
  geom_point(aes(y = win_cnt, colour = 'WindowRaw')) +
  geom_line(aes(y = `50%`, colour = 'STADEM')) +
  scale_colour_manual(values = c(STADEM = 'black',Trap = 'red', WindowAdj = 'blue', WindowRaw = 'skyblue'), labels = c('STADEM','Trap', 'Window (adj.)', 'Window (raw)')) +
  facet_wrap(~ spawn_year, scales = 'free') +
theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Total ',spp, ' Abundance at Lower Granite Dam'),
       x = 'Date',
       y = 'Abundance',
       colour = 'Estimate')
lgr_win_fig

ggsave(paste0('Figures/LGR_win_abund',spp,'.png'),
       lgr_win_fig,
       width = 7,
       height = 5)

# Weekly Natural-origin
lgr_nat_fig <- allLGR %>%
  filter(grepl('X.new.wild', param)) %>%
  mutate(week = as.integer(str_extract(param, '\\d+'))) %>%
  left_join(allLGRtags %>%
              group_by(spawn_year, week) %>%
              summarise(dabom_tags = n()),
            by = c('spawn_year', 'week')) %>%
  left_join(allLGRdata, by = c('spawn_year', 'week' = 'week_num')) %>%
  mutate(p_tagged = dabom_tags/`50%`) %>%
  ggplot(aes(x = Start_Date)) +
  geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .2)+
    geom_line(aes(y = `50%`)) +
    geom_point(aes(y = `50%`, fill = trap_open), colour = 'black', shape = 21, size = 2) +
    geom_bar(aes(y = dabom_tags),fill = 'darkblue', stat='identity') +
    facet_wrap(~ spawn_year, scales = 'free_x') +
  scale_fill_manual(labels = c('Closed', 'Open'), values = c('white', 'darkblue')) +
theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Unique Natural-origin ', spp, ' Weekly Abundance at Lower Granite Dam'),
       subtitle = 'Vertical bars represent the number of natural-origin fish tagged and released for DABOM.',
       x = 'Date',
       y = 'Abundance',
       fill = 'Trap Operation')
lgr_nat_fig

ggsave(paste0('Figures/LGR_nat_week',spp,'.png'),
       lgr_nat_fig,
       width = 7,
       height = 5)


# Tag Rate
# Weekly Natural-origin
lgr_tag_fig <- allLGR %>%
  filter(grepl('X.new.wild', param)) %>%
  mutate(week = as.integer(str_extract(param, '\\d+'))) %>%
  left_join(allLGRtags %>%
              group_by(spawn_year, week) %>%
              summarise(dabom_tags = n()), by = c('spawn_year', 'week')) %>%
  left_join(allLGRdata, by = c('spawn_year', 'week' = 'week_num')) %>%
  mutate(p_tagged = dabom_tags/`50%`) %>%
  ggplot(aes(x = Start_Date)) +
  geom_line(aes(y = p_tagged)) +
  geom_point(aes(y = p_tagged)) +
  facet_wrap(~ spawn_year, scales = 'free') +
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0(spp, ' Natural-Origin Tagging Rate at Lower Granite Dam'),
       subtitle = 'Tagged proportions greater than 1.0 suggest tagging reascension fish.',
       x = 'Date',
       y = 'Proportion Tagged',
       colour = 'Origin Group')
lgr_tag_fig

ggsave(paste0('Figures/LGR_tag_rate',spp,'.png'),
       lgr_tag_fig,
       width = 7,
       height = 5)

# LGR Data for Comparison of GSI to Obs sites.
gsi_dat <- allLGRtags %>%
  left_join(grp_df, by = c('PtagisEventLastSpawnSite' = 'SiteID')) %>%
  select(spawn_year, MasterID, LGDNumPIT, CollectionDate, week,
         MPG, TRT, SiteID = PtagisEventLastSpawnSite, GenSex, BioScaleFinalAge, GenStock, GSI_Group) %>%
  mutate(date = format(as.Date(CollectionDate, '%y-%m-%d'), '%m-%d'),
         date = lubridate::ymd(paste0('2020-',date)),
         sy = as.integer(spawn_year),
         MPG = fct_relevel(MPG, c('Dry Clearwater', 'Wet Clearwater', 'Lower Snake', 'Upper Salmon River', 'Grande Ronde / Imnaha', 'Middle Fork Salmon River', 'South Fork Salmon River', NA))) 
  #filter(spawn_year == 2017)

# Reascencion and Night-Passage Proportions

lgr_passage <- allLGR %>%
  filter(grepl('X.night.wild', param) | grepl('X.reasc.wild', param)) %>%
  mutate(week = as.integer(str_extract(param, '\\d+')),
         grp = case_when(
           str_detect(param, 'night') ~ 'Night Passage',
           str_detect(param, 'reasc') ~ 'Reascension'
         )) %>%
  left_join(allLGRdata, by = c('spawn_year', 'week' = 'week_num')) %>%
  ggplot(aes(x = Start_Date, y = `50%`)) +
  geom_line(aes(colour = grp, group = grp)) +
  scale_colour_viridis_d(end = .75) +
  #scale_x_date(labels = format("%b-%d")) +
  facet_wrap(~ spawn_year, scales = 'free') +
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0(spp, ' Weekly Natural-Origin Night-time Passage and Reascension'),
       x = 'Date',
       y = 'Abundance',
       colour = '')
lgr_passage
  
ggsave(paste0('Figures/LGR_passage',spp,'.png'),
       lgr_passage,
       width = 7,
       height = 5)

#------------------------------------------------------------------------------
# Tributary and Population Estimates
#------------------------------------------------------------------------------
# Travel Time to arrays
allproc_ch <- as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x){
           load(paste0('data/DABOMready/LGR_',spp,'_',x[1],'.rda'))
           
           proc_ch <- proc_list$ProcCapHist %>%
             filter(UserProcStatus) %>%
             select(-UserComment)
           
           proc_ch %>% left_join(grp_df %>% select(-Node), by = 'SiteID') %>%
             left_join(distinct(configuration, SiteID, RKMTotal), by = 'SiteID')
           
         })

trt_ord <- unique(trt_arrival_df$TRT[order(trt_arrival_df$RKMTotal)])
mpg_ord <- c('Dry Clearwater', 'Wet Clearwater', 'Lower Snake', 'Upper Salmon River', 'Grande Ronde / Imnaha', 'Middle Fork Salmon River', 'South Fork Salmon River')

trt_arrival_df <- allproc_ch %>%
  filter(!is.na(TRT)) %>%
  group_by(TagID) %>%
  slice(which.min(ObsDate)) %>%
  mutate(year = case_when(
    spp == 'Steelhead' & lubridate::yday(TrapDate) >= 183 ~ 2019,
    TRUE ~ 2020),
    lgr_date = format(as.Date(TrapDate, '%y-%m-%d %hh:%mm:%ss'), '%m-%d'),
    lgr_date = lubridate::ymd(paste0(year,lgr_date)),
    date = format(as.Date(ObsDate, '%y-%m-%d %hh:%mm:%ss'), '%m-%d'),
    date = lubridate::ymd(paste0(year,date)),
    sy = as.integer(spawn_year),
    #MPG = factor(MPG, levels = mpg_ord),
    travel_time = difftime(ObsDate,TrapDate, units = 'days'))

# Arrival/run timing to LGR
lgr_run <- trt_arrival_df %>%
  ggplot(aes(x = lgr_date, group = TRT, colour = MPG)) +
  stat_ecdf( size = 1) +
  facet_grid(spawn_year ~ .) +
  scale_colour_viridis_d(end = 1, option = 'D') +
  scale_x_date(date_labels = format('%b-%d')) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Natural-Origin ',spp, ' Arrival Timing at Lower Granite'),
       subtitle = '.',
       x = 'Date',
       y = 'ECDF',
       colour = 'Major Population Group')

lgr_run

ggsave(paste0('Figures/LGR_run_timing',spp,'.png'),
       lgr_run,
       width = 7,
       height = 9)


travel_fig <- trt_arrival_df %>%
  ggplot(aes(x = travel_time, group = TRT, fill = MPG, colour = MPG)) +
  #geom_density() +
  stat_ecdf(size = 1) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  facet_grid(spawn_year ~ .) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Natural-Origin ',spp, ' Travel Time to Spawn Location'),
       x = 'Days from LGR Trap Date',
       y = 'ECDF',
       colour = '')

travel_fig

ggsave(paste0('Figures/travel_time_',spp,'.png'),
       travel_fig,
       width = 7,
       height = 9)

# 
# speed_fig <- trt_arrival_df %>%  
#   mutate(TRT = factor(TRT, levels = ord)) %>%
#     ggplot(aes(x = travel_time, y = RKMTotal, colour = TRT)) +
#   geom_jitter(alpha = .125, height = 0, width = 0)  +
#   geom_smooth(method = 'lm', se = FALSE) +
#   scale_colour_viridis_d(end = 1) +
#   #facet_wrap( ~ spawn_year, scales = 'free_x') +
#   theme_bw() +
#   theme(strip.text.y = element_text(angle = .45)) +
#   theme(legend.position = 'bottom') +
#   labs(title = paste0('Natural-Origin ', spp, ' Travel Speed to Spawn Location'),
#        x = 'Travel Time (Days)',
#        y = 'Total River Kilometer',
#        colour = 'Major Population Group')
# 
# speed_fig
# 
# ggsave(paste0('Figures/LGR_arrival_',spp,'.png'),
#        trt_arrival_fig,
#        width = 8,
#        height = 8)

# Detection -------------------------------------------------------------------
# make an observed vs predicted of detection probs.
allEffDf = as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x) {
           load(paste0('data/DABOMready/LGR_',spp,'_',x[1],'.rda'))
           load(paste0('Abundance_results/LGR_Summary_', spp, '_', x[1], '.rda')) 
           
           eff <- estNodeEff(filter(proc_list$ProcCapHist,UserProcStatus), node_order = proc_list$NodeOrder)
           
           trt_df <- assign_POP_GSI(spp, configuration, site_df)[[1]] %>% as_tibble()
           
           detect_summ %>%
             left_join(eff, by = 'Node') %>%
             left_join(proc_list$NodeOrder %>%
                         select(Node, Group), by = 'Node') %>%
             left_join(trt_df %>%
                         select(Node, MPG, TRT) %>% distinct() %>%
                         as_tibble(), by = 'Node')
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
  labs(title = paste0(spp, ' Detection Probabilities'),
       x = 'Predicted with PITcleanR',
       y = 'Predicted with DABOM',
       size = 'Unique Tags Observed')

ObsPred_detect_p

ggsave(paste0('Figures/ObsVsPred_DetectProp_',spp,'.png'),
       ObsPred_detect_p,
       width = 8,
       height = 8)

detect_p <- allEffDf %>%
  filter(estimate != 0 & estimate != 1) %>%
  filter(!is.na(Group)) %>%
  ggplot(aes(x = Node,
             y = estimate,
             color = spawn_year,
             group = spawn_year)) +
  geom_linerange(aes(ymin = lowerCI,
                     ymax = upperCI),
                 position = position_dodge(width = .5)) +
  geom_point(size = 2, position = position_dodge(width = .5)) +
  scale_colour_viridis_d() +
  coord_flip() +
  facet_wrap(~Group, scales = 'free', ncol = 4) +
  #facet_grid(MPG~., scales = "free_y", space = 'free_y', labeller = label_wrap_gen(width = 10)) #+
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Natural-Origin ',spp, ' Detection Probabilities at DABOM Nodes'),
       x = 'Node',
       y = 'Detection Probability',
       colour = 'Spawn Year')

detect_p

ggsave(paste0('Figures/DetectProp_',spp,'.png'),
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
           
           #trt_df <- site_trt_designations(spp, configuration, site_df)
           
           N_pop_summ %>%
           left_join(map_df %>%
                       select(MPG, TRT_POPID), by = c('TRT' = 'TRT_POPID')) %>%
             filter(median != 0)
         })

# Get pop names for time series ends
ts_names <- allAbundDf %>%
  group_by(MPG, TRT) %>%
  top_n(-1,spawn_yr) %>%
  select(spawn_yr, MPG, TRT, median)

pop_N <- allAbundDf %>%
  # mutate(TRT = factor(TRT, levels = trt_ord),
  #        MPG = factor(MPG, levels = mpg_ord)) %>%
  filter(!is.na(TRT)) %>%
  filter(!(TRT == 'CRLOC-s' & spawn_year <= 2016)) %>%
  group_by(TRT) %>%
  mutate(zscore = (median - mean(median))/sd(median)) %>%
  ggplot(aes(x = spawn_yr,
             y = median,
             #colour = TRT,
             group = TRT)) +
  geom_ribbon(aes(ymin = lowerCI,
                     ymax = upperCI),alpha = .2) +
  geom_line(aes(colour = MPG)) +
  geom_point(aes(colour = MPG), size = 2) +#, position = position_dodge(width = .5)) +
  #geom_label(data = ts_names, aes(x = spawn_yr, y = median, label = TRT), size = 2, position = position_nudge(y = 200, x=.25)) +
  scale_colour_viridis_d() +
  #scale_y_continuous(sec.axis = sec_axis(~ ., breaks = ts_names)) +
  # coord_flip() +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_wrap(~TRT, labeller = label_wrap_gen(width = 30), scales = 'free_y', ncol = 4) + #space = 'free_y', 
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Natural-Origin ',spp, ' Population Abundance'),
       x = 'Spawn Year',
       y = 'Abundance',
       colour = 'Major Population Group')

pop_N

ggsave(paste0('Figures/PopAbund_',spp,'.png'),
       pop_N,
       width = 8,
       height = 10)

# Sex ------------------------------------------------------------------------------
# make an observed vs predicted female proportion plot, by population
allSexDf = as.list(yr_range) %>%
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
             mutate_at(vars(MPG, TRT, Group),
                       list(fct_explicit_na),
                       na_level = 'Overall') %>%
             mutate(TRT = fct_reorder(TRT, mean),
                    TRT = fct_relevel(TRT, 'Overall',
                                      after = 0))
         })

sex_ObsVsPred_p = allSexDf %>%
  filter(!is.na(propF)) %>%
  ggplot(aes(x = propF,
             y = mean)) +
  geom_abline(linetype = 2) +
  #geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  #geom_errorbarh(aes(xmin = propF - 1.96*propF_se, xmax = propF + 1.96*propF_se)) +
  geom_point(aes(colour = spawn_year, size = nSexed)) +
  scale_color_viridis_d() +
  xlim(0,1) +
  ylim(0,1) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  facet_wrap(~ TRT) +
             #scales = 'free') +
  labs(x = 'Observed',
       y = 'Predicted',
       title = paste0(spp,' Female Proportions'),
       color = 'Spawn Year',
       size = '# Sexed')

sex_ObsVsPred_p

ggsave(paste0('Figures/ObsVsPred_SexProp_',spp,'.png'),
       sex_ObsVsPred_p,
       width = 8,
       height = 8)

sex_Pred <- allSexDf %>%
  #mutate(MPG = factor(MPG, levels = c(mpg_ord,'Overall'))) %>%
  filter(TRT != 'Overall') %>%
  #filter(spawn_year == '2018') %>%
  ggplot(aes(x = TRT,
             y = mean,
             color = MPG,
             group = TRT)) +
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
  facet_grid(spawn_year ~., scales = "free_y", space = 'free_y', labeller = label_wrap_gen(width = 10)) +
  labs(x = '',
       y = 'Female Proportion',
       title = paste0('Natural-Origin ',spp, ' Female Proportions'),
       color = 'Major Population Group')

sex_Pred

ggsave(paste0('Figures/SexProp_',spp,'.png'),
       sex_Pred,
       width = 7,
       height = 9)

# Age ------------------------------------------------------------------------------
# make an observed vs predicted age proportion plot, by population
allAgeDf = as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x) {
           load(paste0('Age_results/Population_AgeProp_', spp, '_', x[1], '.rda'))
           
           poss_ages <- modAgeDf %>%
             filter(!is.na(TRT)) %>%
             select(starts_with('age'))
           
           poss_ages <- as_tibble(na.omit(as.numeric(
             gsub('[[:alpha:]]','', names(poss_ages[,!colSums(poss_ages) == 0])))))  %>%
             mutate(age_fct = as.integer(as.factor(value))) %>%
             rename(age = value)
           
           age_mod$summary %>%
             as_tibble(rownames = 'param') %>%
             filter(grepl('^pi', param)) %>%
             mutate(popNum = str_extract(param, '[:digit:]+'),
                    age_fct = str_split(param, ',', simplify = T)[,2],
                    age_fct = str_remove(age_fct, '\\]')) %>%
             mutate_at(vars(popNum, age_fct),
                       list(as.integer)) %>%
            left_join(poss_ages, by = 'age_fct') %>%
            select(-age_fct) %>%
            # mutate(age_fct = age_fct + 1) %>%
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
  xlim(0,1) +
  ylim(0,1) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  #scale_color_brewer(palette = 'Set1') +
  scale_size_continuous(breaks = c(10,50,100,150)) +
  scale_colour_viridis_d() +
  facet_wrap(~ TRT) + #,
             #scales = 'free') +
  labs(x = 'Observed',
       y = 'Predicted',
       title = paste0('Natural-Origin ',spp, ' Age Proportions'),
       color = 'Age',
       size = '# Aged')

age_ObsVsPred_p

ggsave(paste0('Figures/ObsVsPred_AgeProp_',spp,'.png'),
       age_ObsVsPred_p,
       width = 8,
       height = 8)

# plot showing estimated age proportions (mean) for each TRT population, by year
ageProp_p = allAgeDf %>%
  select(spawn_year, MPG, TRT, age, mean) %>%
  distinct() %>%
  ggplot(aes(x = fct_rev(TRT),
             y = mean,
             fill = as.factor(age))) +
  geom_bar(stat = 'identity',
           position = position_stack(reverse = T)) +
  #scale_fill_brewer(palette = 'Set1') +
  scale_fill_viridis_d() +
  scale_y_continuous(breaks = c(0,.5,1), labels = c(0,.5,1)) +
  coord_flip() +
  theme_bw() +
  theme(#axis.text.y = element_text(size = 5),
        strip.text.y = element_text(angle = .45),
        legend.position = 'bottom') +
  facet_grid(MPG ~ spawn_year, scales = "free", space = 'free', labeller = label_wrap_gen(width = 10)) +
  #facet_wrap(~ spawn_year) +
  labs(title = paste0('Natural-Origin ',spp, ' Age Proportions'),
       fill = 'Age',
       x = 'TRT',
       y = 'Proportion')

ageProp_p

ggsave(paste0('Figures/AgePropEst_',spp,'.png'),
       ageProp_p,
       width = 8,
       height = 8)

ageProp2_p = allAgeDf %>%
  select(spawn_year, MPG, TRT, age, mean) %>%
  distinct() %>%
  ggplot(aes(x = spawn_year,
             y = mean,
             fill = as.factor(age))) +
  geom_bar(stat = 'identity',
           position = position_stack(reverse = T)) +
  #scale_fill_brewer(palette = 'Set1') +
  scale_y_continuous(breaks = c(0,.5,1), labels = c(0,.5,1)) +
  scale_fill_viridis_d() +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = 'bottom') +
  facet_wrap(~ fct_rev(TRT)) +
  labs(fill = 'Age',
       title = paste0('Natural-Origin ',spp, ' Age Proportions'),
       x = 'Spawn Year',
       y = 'Proportion')

ageProp2_p

ggsave(paste0('Figures/AgePropEst2_',spp,'.png'),
       ageProp2_p,
       width = 8,
       height = 8)

muVecDf = as.list(yr_range) %>%
  rlang::set_names() %>%
  map_df(.id = 'spawn_year',
         .f = function(x) {
           load(paste0('Age_results/Population_AgeProp_', spp, '_', x[1], '.rda'))

           age_mod$summary %>%
             as_tibble(rownames = 'param') %>%
             filter(grepl('^avgPi', param)) %>%
             mutate(runType = str_extract(param, '[:digit:]'),
                    runType = if_else(runType == 1, 'A-Run', 
                                      if_else(runType == 2, 
                                              'B-Run', as.character(NA))),
                    age = str_sub(param, -2, -2),
                    age = as.integer(age),
                    age = age + 2) %>%
             mutate_at(vars(runType),
                       list(as.factor)) %>%
             mutate_at(vars(age),
                       list(~as.factor(as.character(.)))) %>%
             mutate(TRT = paste0(runType, '-run'),
                    MPG = paste0(runType, '-run'))
         })

muProp_p = muVecDf %>%
  select(spawn_year, runType, age, mean) %>%
  mutate(runType = ifelse(spp == 'Chinook', '',runType)) %>%
  distinct() %>%
  ggplot(aes(x = fct_rev(runType),
             y = mean,
             fill = as.factor(age))) +
  geom_bar(stat = 'identity',
           position = position_stack(reverse = T)) +
  #scale_fill_brewer(palette = 'Set1') +
  scale_fill_viridis_d() +
  coord_flip() +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        legend.position = 'right') +
  facet_wrap(~ spawn_year) +
  labs(fill = 'Age',
       title = paste0('Mean Snake Basin Natural-Origin ',spp, ' Age Proportion'),
       x = '',
       y = 'Proportion')

muProp_p

ggsave(paste0('Figures/AvgAgePropEst_mu_',spp,'.png'),
       muProp_p,
       width = 8,
       height = 8)

# Spawner-Recruits----
library(readxl)

allSR_Df <- read_excel(paste0('Abundance_results/LGR_AllSummaries_',spp,'.xlsx'), sheet = 'Pop Stock Recruit') %>%
  mutate(TRT = gsub('\"',"",TRT)) %>%
  filter(median != 0)

lambda_p <- allSR_Df %>%
  filter(variable == 'lambda') %>%
  ggplot(aes(x = TRT,
             y = median,
             colour = as.factor(brood_yr),
             group = brood_yr)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_linerange(aes(ymin = lowerCI, ymax = upperCI),
                 position = position_dodge(width = .5)) +
  geom_point(size = 4, position = position_dodge(width = .5)) +
  scale_colour_viridis_d() +
  coord_flip() +
  #facet_wrap(~TRT, scales = 'free_y') +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(title = paste0('Natural-Origin ',spp, ' Productivity'),
       caption = '** Productivity estimates exclude hatchery-origin spawners but includes their progeny as recruits.',
       x = 'Spawn Year',
       y = 'Recruits/Spawner',
       colour = 'Brood Year')

lambda_p  

ggsave(paste0('Figures/lambda_',spp,'.png'),
       lambda_p,
       width = 8,
       height = 8)

SR_df <- allSR_Df %>%
  filter(variable == 'Spawners') %>%
  select(brood_yr, species, TRT, spawners = median, s_lower = lowerCI, s_upper = upperCI) %>%
  left_join(allSR_Df %>%
              filter(variable == 'Recruits') %>%
              select(brood_yr, species, TRT, recruits = median, r_lower = lowerCI, r_upper = upperCI)) %>%
  left_join(grp_df %>% select(MPG, TRT), by = 'TRT') #%>%
  # mutate(TRT = factor(TRT, levels = trt_ord),
  #        MPG = factor(MPG, levels = mpg_ord))

SR_plot <- SR_df %>%
  #filter(brood_yr <= 2012) %>%
  ggplot(aes(x = spawners,
             y = recruits,
             colour = MPG)) +
  geom_abline(linetype = 2) +
  geom_errorbar(aes(ymin = r_lower,
                     ymax = r_upper)) +
  geom_errorbarh(aes(xmin = s_lower,
                     xmax = s_upper)) +
  geom_point(size = 2) +
  #xlim(0,4000) +
  #ylim(0,4000) +
  coord_fixed() +
  scale_colour_viridis_d() +
  facet_wrap(~ TRT) +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Stock',
       y = 'Recruits',
       title = paste0('Natural-Origin ',spp, ' Stock-Recruit Points'),
       colour = 'Major Population Group')

SR_plot

ggsave(paste0('Figures/stock_recruit_',spp,'.png'),
       SR_plot,
       width = 8,
       height = 8)


rec <- SR_df %>%
  filter(brood_yr == 2010) %>%
  select(TRT, max_recruit = recruits)

recruits <- SR_df %>%
  # #left_join(rec, by = 'TRT') %>%
  # group_by(TRT) %>%
  # mutate(p_recruit = recruits/max_recruit,
  #        rel_diff = (recruits - mean(recruits))/mean(recruits)) %>%
  ggplot(aes(x=brood_yr, y = recruits, group = TRT)) + #recruits)) +
  geom_ribbon(aes(ymin = r_lower, ymax = r_upper), alpha = .25) +
  geom_line(aes(colour = MPG)) +
  geom_point(aes(colour = MPG), size = 2) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~TRT, scales = 'free_y') +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Brood Year',
       y = 'Recruits',
       title = paste0('Natural-Origin ',spp, ' Adult Recruits'),
       colour = 'Major Population Group')

recruits  

ggsave(paste0('Figures/adult_recruits_',spp,'.png'),
       recruits,
       width = 8,
       height = 8)

recruits_p <- SR_df %>%
  right_join(rec, by = 'TRT') %>%
  group_by(TRT) %>%
  mutate(p_recruit = recruits/max_recruit,
          rel_diff = (recruits - mean(recruits))/mean(recruits)) %>%
  filter(TRT != 'SNTUC') %>%
  ggplot(aes(x=brood_yr, y = p_recruit, group = TRT)) + #recruits)) +
  geom_line(aes(colour = MPG)) +
  geom_point(aes(colour = MPG), size = 2) +
  expand_limits(y = 0) +
  scale_y_continuous(expand = c(0,0)) +
  #facet_wrap(~TRT, scales = 'free_y') +
  scale_colour_viridis_d() +
  theme_bw() +
  theme(legend.position = 'bottom') +
  labs(x = 'Brood Year',
       y = 'Recruits',
       subtitle = 'Calculated as the proportion of brood year 2010 recruits.',
       title = paste0('Natural-Origin ',spp, ' Adult Recruits'),
       colour = 'Major Population Group')

recruits_p

ggsave(paste0('Figures/adult_recruits_p_',spp,'.png'),
       recruits_p,
       width = 8,
       height = 8)
