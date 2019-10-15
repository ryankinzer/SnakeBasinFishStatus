# Author: Kevin See
# Purpose: diagnostics on DABOM MCMC runs
# Created: 10/7/2019
# Last Modified: 10/7/19
# Notes: postpack package can be installed by typing "devtools::install_github("bstaton1/postpack")"

# load needed libraries
library(tidyverse)
library(coda)
library(jagsUI)
library(postpack)
library(mcmcr)
library(ggmcmc)
library(ggpubr)

theme_set(theme_bw())

#---------------------------------------
# choose species and year combination
spp = c('Chinook', 'Steelhead')[2]
yr = 2013

#---------------------------------------
# where are results stored?
dabomFolder = 'DABOM_results'

# load DABOM model output
load(paste0(dabomFolder,'/LGR_DABOM_', spp, '_', yr, '.rda'))

# pull out mcmc.list object
my_mod = dabom_mod

#---------------------------------------
# using postpack
# what parameters were tracked?
get_p(my_mod,
      type = 'base')

# some summary statistics
post_summ(my_mod,
          '_p$') %>%
  t() %>%
  as_tibble(rownames = 'param') %>%
  filter(!grepl('p_pop_main', param))

#---------------------------------------
# using mcmcr
anyNA(my_mod)
my_mcmcr = as.mcmcr(my_mod)

# get Rhat statistics for all parameters
rhat_df = rhat(my_mcmcr,
               by = 'parameter',
               as_df = T) %>%
  mutate(type = if_else(grepl('_p$', parameter),
                        'Detection',
                        if_else(grepl('^p_pop', parameter) |
                                  grepl('^phi', parameter),
                                'Movement',
                                'Other')))

# plot histogram of Rhat statistics
rhat_df %>%
  ggplot(aes(x = rhat)) +
  geom_histogram(fill = 'blue',
                 bins = 40) +
  facet_wrap(~ type,
             scales = 'free')


# which parameters have converged and which haven't?
convg_df = converged(my_mcmcr, 
                     by = 'parameter',
                     as_df = T)

# look at parameters that have not converged
convg_df %>%
  filter(!converged) %>%
  left_join(rhat_df)

param_chk = convg_df %>%
  filter(!converged) %>%
  # don't look at p_pop_main, it contains > 1000 separate parameters
  filter(parameter != 'p_pop_main') %>%
  pull(parameter)


# use character vector of different parameters here if you want
 param_chk = 'ZENB0_p$'

# diagnostic plots with postpack package
# you can save these plots by adding an arguement save = T, and file = your file path
diag_plots(post = my_mod,
           p = param_chk)


# other diagnostic plots with ggmcmc
my_ggs1 = ggs(my_mod,
              param_chk[1])

my_ggs = param_chk %>%
  as.list() %>%
  map_df(.f = function(x) {
    ggs(my_mod,
        x[1])
  })
for(my_attr in c('nChains', 'nParameters', 'nIterations', 'nBurnin', 'nThin', 'description', 'class')) {
  attr(my_ggs, my_attr) = attr(my_ggs1, my_attr)
}

# save a file with lots of different diagnostic plots
ggmcmc(my_ggs,
       file = 'Figures/DABOM_diagnostic_plots.pdf',
       param_page = 10)


# a few specific plots to look at
dens_p = ggs_density(my_ggs)
trace_p = ggs_traceplot(my_ggs)
run_mean_p = ggs_running(my_ggs)
rhat_p = ggs_Rhat(my_ggs)
geweke_p = ggs_geweke(my_ggs)
ggs_autocorrelation(my_ggs)

ggarrange(plotlist = list(dens_p,
                          trace_p),
          ncol = 2,
          common.legend = T,
          legend = 'bottom')


#---------------------------------------
get_p(my_mod,
      type = 'base')

# family name of parameters to plot diagnostics for
param_fam = "p_pop"
param_fam = '_p$'
param_fam = "sigma"
param_fam = 'phi'
param_fam = 'p_pop_main'
param_fam = c('JOSEPC', 'JOC', 'phi_josepc')

diag_plots(my_mod,
           param_fam)

#---------------------------------------
# use tools from Shiny STAN
#---------------------------------------
library(shinystan)

# lauch shinystan, then go to Diagnose, then Rhat, n_eff, se_mean tab
my_mod %>%
  as.shinystan %>%
  launch_shinystan()
