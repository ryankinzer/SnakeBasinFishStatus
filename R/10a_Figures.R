# New Figures....start from output or cdmsR

library(tidyverse)
library(lubridate)
library(PITcleanr)
library(scales)

source('./R/theme_rk.R')

spp <- 'Chinook'

stadem_df <- readxl::read_excel(paste0('./Abundance_results/LGR_AllSummaries_',spp,'_v2.xlsx'),
                                sheet = 'LGR Esc') %>%
  mutate(spawn_year = as.numeric(spawn_year))


pop_df <- readxl::read_excel(paste0('./Abundance_results/LGR_AllSummaries_',spp,'_v2.xlsx'),
                              sheet = 'Pop Total Esc') 

site_df <- readxl::read_excel(paste0('./Abundance_results/LGR_AllSummaries_',spp,'_v2.xlsx'),
                                sheet = 'Site Esc') 

stadem_fig <-  stadem_df %>%
  filter(origin != 'Total') %>%
  complete(spawn_year = full_seq(spawn_year, period = 1), nesting(origin)) %>%
  ggplot(aes(x = spawn_year, y = estimate, group = origin)) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha = .2) +
  geom_line() +
  geom_point() +
  scale_y_continuous(labels = comma, breaks = seq(0,300000, by = 50000))+
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_wrap(~origin) +
  labs(title = paste0('Natural-origin ', spp, ' Escapement to Lower Granite Dam'),
       subtitle = 'Grey shaded area represents the 95% confidence interval band.',
       caption = 'Data Source: STAte-space Dam Escapement Model (STADEM; See et al. 2021) estimates produced by the Nez Perce Tribe',
       x = 'Spawn Year',
       y = 'Escapement') +
  theme_rk()

stadem_fig

ggsave(paste0('./Figures/',spp,'_stadem_estimates.png'), stadem_fig, width = 11, height = 8.5)


pop_fig <- pop_df %>%
  filter(valid_est == 1) %>%
  complete(spawn_yr = full_seq(spawn_yr, period = 1), nesting(TRT)) %>%
  ggplot(aes(x = spawn_yr, y = median, group = TRT)) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha = .2) +
  geom_line() +
  geom_point() +
  #scale_color_viridis_d() +
  scale_y_continuous(labels = comma)+
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_wrap(~TRT, scales = 'free_y') +
  theme_rk() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom',
        axis.text = element_text(size=8)) +
  labs(title = paste0('Natural-origin ', spp, ' Escapement to Snake River Populations'),
       subtitle = 'Grey shaded area represents the 95% confidence interval band.',
       caption = 'Data Source: In-stream pit-tag array DABOM estimates produced by the Nez Perce Tribe',
       x = 'Spawn Year',
       y = 'Abundance',
       colour = 'Site')

pop_fig

ggsave(paste0('./Figures/',spp,'_dabom_pop_estimates.png'), pop_fig, width = 11, height = 8.5)


site_fig <- site_df %>%
  filter(site %in% c('SFG', 'ESS', 'KRS', 'ZEN')) %>%
  complete(spawn_yr = full_seq(spawn_yr, period = 1), nesting(site)) %>%
  ggplot(aes(x = spawn_yr, y = estimate, group = site)) +
  geom_ribbon(aes(ymin = lowerCI, ymax = upperCI), alpha = .2) +
  geom_line() +
  geom_point() +
  #scale_color_viridis_d() +
  scale_y_continuous(labels = comma)+
  scale_x_continuous(breaks = pretty_breaks()) +
  facet_wrap(~site) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = .45)) +
  theme(legend.position = 'bottom') +
  labs(title = paste0(spp, ' Abundance in Secesh River'),
       subtitle = 'Grey shaded area represents the 95% confidence interval band.',
       caption = 'Data Source: In-stream pit-tag array estimates produced by the Nez Perce Tribe',
       x = 'Spawn Year',
       y = 'Abundance',
       colour = 'Site')

site_fig

ggsave('./Figures/secesh_estimates.png', site_fig, width = 7, height = 5)
