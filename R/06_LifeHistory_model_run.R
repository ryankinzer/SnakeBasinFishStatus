# Author: Kevin See
# Purpose: Run life history models to estimate sex ratio and age structure.
# Created: 7/10/2019
# Last Modified: 7/10/19
# Notes: Rick Orme provided feedback on PITcleanr output on 7/1/19


# Load Packages
library(tidyverse)
library(readxl)
library(jagsUI)

# set up folder structure
SexFolder = 'Sex_results' # for processed files
if(!dir.exists(SexFolder)) {
  dir.create(SexFolder)
}

AgeFolder = 'Age_results' # for biologist cleaned files and dabom data
if(!dir.exists(AgeFolder)) {
  dir.create(AgeFolder)
}

#------------------------------------------------------------------------------
# File path to life history summaries
filepath <- './data/LifeHistoryData/'

# Set species and Spawn Year
spp = 'Steelhead'
#yr = 2018

#------------------------------------------------------------------------------
# Create JAGs model to estimate female proportion

sexModelNm = './ModelFiles/FemalePropJAGS.txt'

modelCode = '
model {
  
  for(i in 1:length(f)) {
    f[i] ~ dbin(p[popNum[i]], tags[i])
  }

  for(j in 1:max(popNum)) {
    p[j] <- ilogit(logit_p[j])
    logit_p[j] ~ dnorm(mu, tau)
  }
  # transform overall mean back to proportion scale
  mu_ilogit <- ilogit(mu)

  mu ~ dnorm(0, 0.001)
  sig ~ dunif(0, 100)
  tau <- pow(sig, -2)

}'

cat(modelCode,
    file = sexModelNm)

#------------------------------------------------------------------------------
# Run Sex Model
# read in data

for(yr in 2010:2018){
  cat(paste('Starting sex model for', spp, 'in run year', yr, '\n'))
  
  modSexDf = read_excel(paste0(filepath, 'LGR_', spp, '_', yr, '.xlsx'),
                        'SexRatio',
                        progress = F)
  
  # pull out relevant bits for JAGS, and name them appropriately
  sex_jagsData = modSexDf %>%
    filter(!is.na(modBranch)) %>%
    mutate(popNum = as.integer(as.factor(TRT))) %>%
    select(f = F,
           tags = nSexed,
           popNum) %>%
    as.list()
  
  # set parameters to save
  jagsParams = c('p', 'mu', 'sig', 'mu_ilogit')
  
  # run JAGS model
  sex_mod = jags(data = sex_jagsData,
                 parameters.to.save = jagsParams,
                 model.file = sexModelNm,
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 5000,
                 n.thin = 20,
                 verbose = F)
  
  save(sex_mod, sex_jagsData, modSexDf,
       file = paste0(SexFolder, '/Population_SexProp_', spp, '_', yr, '.rda'))
}



#------------------------------------------------------------------------------
# Create JAGs model to estimate age proportions

ageModelNm = './ModelFiles/AgePropJAGS_hier.txt'

modelCode = '
data {
  D <- dim(ageMat)
}

model {

  for(i in 1:D[1]) {
    ageMat[i,] ~ dmulti(pi[popNum[i],], tags[i])
  }
  
  # multivariate logistic normal transformation to make it hierarchical
  for(j in 1:max(popNum)) {
    p[j, 1] <- 0
    p[j,2:D[2]] ~ dmnorm(mu[1:(D[2] - 1)], Tau[1:(D[2] - 1), 1:(D[2] - 1)])
    
    sum_exp_p[j] <- sum(exp_p[j,])
    
    for(k in 1:D[2]) {
      exp_p[j,k] = exp(p[j, k])
      pi[j, k] <- exp(p[j, k]) / sum_exp_p[j]
    }
  }
  
  # transform mu back to proportions
  muProp[1] = 0
  for(i in 2:D[2]) {
    muProp[i] = mu[i-1]
  }
  sum_exp_mu = sum(exp_mu)
  for(i in 1:D[2]) {
    exp_mu[i] = exp(muProp[i])
    avgPi[i] = exp_mu[i] / sum_exp_mu
  }
  
  # Cauchy prior on the MVN mean vector
  for(i in 1:(D[2] - 1)) {
    mu[i] ~ dt(0, 0.001, 1)
  }
  # Priors on the precision matrix
  Tau ~ dwish(R, k)
  k <- D[2] + 1
}'

cat(modelCode,
    file = ageModelNm)

#------------------------------------------------------------------------------
# Run Age Model
# read in data

for(yr in 2010:2018){
  cat(paste('Starting age model for', spp, 'in run year', yr, '\n'))
  
  modAgeDf = read_excel(paste0(filepath, 'LGR_', spp, '_', yr, '.xlsx'),
                        'AgeFreq',
                        progress = F)
  
  # pull out relevant bits for JAGS, and name them appropriately
  age_jagsData = modAgeDf %>%
    filter(!is.na(modBranch)) %>%
    mutate(popNum = as.integer(as.factor(TRT))) %>%
    select(tags = nAged,
           popNum) %>%
    as.list()
  
  age_jagsData$ageMat = modAgeDf %>%
    filter(!is.na(modBranch)) %>%
    select(starts_with('age')) %>%
    as.matrix()
  
  # drop ages with no observed fish in them
  if(sum(colSums(age_jagsData$ageMat) == 0) > 0) {
    age_jagsData$ageMat = age_jagsData$ageMat[,!colSums(age_jagsData$ageMat) == 0]
  }
  
  age_jagsData$R = diag(1, ncol(age_jagsData$ageMat) - 1)
  
  
  # set parameters to save
  jagsParams = c('pi', 'mu', 'Tau', 'avgPi')
  
  # run JAGS model
  age_mod = jags(data = age_jagsData,
                 parameters.to.save = jagsParams,
                 model.file = ageModelNm,
                 n.chains = 3,
                 n.iter = 20000,
                 n.burnin = 10000,
                 n.thin = 20,
                 verbose = T)
  
  save(age_mod, age_jagsData, modAgeDf,
       file = paste0(AgeFolder, '/Population_AgeProp_', spp, '_', yr, '.rda'))
  
  rm(age_mod, age_jagsData, modAgeDf)
}

#------------------------------------------------------------------------------
# make an observed vs predicted age proportion plot, by population
allAgeDf = as.list(2010:2018) %>%
  rlang::set_names() %>%
  map_df(.id = 'Year',
         .f = function(x) {
           load(paste0(AgeFolder, '/Population_AgeProp_', spp, '_', x[1], '.rda'))
           
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
