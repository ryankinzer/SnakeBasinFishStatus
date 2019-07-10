# Author: Kevin See
# Purpose: Run life history models to estimate sex ratio and age structure.
# Created: 7/10/2019
# Last Modified: 7/10/19
# Notes: Rick Orme provided feedback on PITcleanr output on 7/1/19


# Load Packages
library(tidyverse)
library(readlx)
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
# Create JAGs model to estimate female proporiton

modelNm = './ModelFiles/FemalePropJAGS.txt'

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
    file = modelNm)

#------------------------------------------------------------------------------
# File path to life history summaries
filepath <- './data/LifeHistoryData/'

# Set species and Spawn Year
spp = 'Steelhead'
#yr = 2010

#------------------------------------------------------------------------------
# Run Sex Model
# read in data

for(yr in 2010:2018){

modSexDf = read_excel(paste0(filepath, 'LGR_', spp, '_', yr, '.xlsx'),
                      'SexRatio',
                      progress = F)

# pull out relevant bits for JAGS, and name them appropriately
jagsData = modSexDf %>%
  filter(!is.na(modBranch)) %>%
  mutate(popNum = as.integer(as.factor(TRT))) %>%
  select(f = F,
         tags = nSexed,
         popNum) %>%
  as.list()

# set parameters to save
jagsParams = c('p', 'mu', 'sig', 'mu_ilogit')

# run JAGS model
sex_mod = jags(data = jagsData,
           parameters.to.save = jagsParams,
           model.file = modelNm,
           n.chains = 3,
           n.iter = 10000,
           n.burnin = 5000,
           n.thin = 20,
           verbose = F)

save(spp, yr, sex_mod, modSexDf,
     file = paste0(SexFolder, '/Population_SexProp_', spp, '_', yr, '.rda'))
}

#------------------------------------------------------------------------------
# Run Sex Model
# read in data
