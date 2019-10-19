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

# Set species
spp = 'Steelhead'

# Set year range based on spp.
if(spp == 'Steelhead'){
  year_range = 2010:2019
} else {
  year_range = 2010:2018
}

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

for(yr in year_range){
  cat(paste('Starting sex model for', spp, 'in run year', yr, '\n'))
  
  modSexDf = read_excel(paste0(filepath, 'LGR_', spp, '_', yr, '.xlsx'),
                        'SexRatio',
                        progress = F)
  
  # pull out relevant bits for JAGS, and name them appropriately
  sex_jagsData = modSexDf %>%
    filter(!is.na(TRT)) %>%
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

ageModel_simp = './ModelFiles/AgePropJAGS_simp.txt'

modelCode = '
data {
  D <- dim(ageMat)
}

model {

  for(i in 1:D[1]) {
    ageMat[i,] ~ dmulti(pi[popNum[i],], tags[i])
  }

  for(k in 1:D[2]) {
    alpha[k] <- 1
  }

  for(j in 1:max(popNum)) {
    pi[j,1:D[2]] ~ ddirch(alpha[1:D[2]])
  }

 }'

cat(modelCode,
    file = ageModel_simp)

ageModel_hier = './ModelFiles/AgePropJAGS_hier.txt'

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
    p[j,2:D[2]] ~ dmnorm(mu[runType[j], 1:(D[2] - 1)], Tau[1:(D[2] - 1), 1:(D[2] - 1)])
    
    sum_exp_p[j] <- sum(exp_p[j,])
    
    for(k in 1:D[2]) {
      exp_p[j,k] = exp(p[j, k])
      pi[j, k] <- exp(p[j, k]) / sum_exp_p[j]
    }
  }
  
  # transform mu back to proportions
  for(j in 1:max(runType)) {
    muProp[j,1] = 0
    for(i in 2:D[2]) {
      muProp[j,i] = mu[j,i-1]
    }
    sum_exp_mu[j] = sum(exp_mu[j,])
    for(i in 1:D[2]) {
      exp_mu[j,i] = exp(muProp[j,i])
      avgPi[j,i] = exp_mu[j,i] / sum_exp_mu[j]
    }
  }
  
  # Cauchy prior on the MVN mean vector
  for(i in 1:(D[2] - 1)) {
    for(j in 1:max(runType)) {
      mu[j, i] ~ dt(0, 0.001, 1)
    }
  }
  # Priors on the precision matrix
  Tau ~ dwish(R, k)
  k <- D[2] + 1

}'

cat(modelCode,
    file =ageModel_hier)

#------------------------------------------------------------------------------
# Run Age Model
# set model as fixed population effect or random population with
# A and B run groups for steelhead.
model <- c('simple','hierarchical')[2]

if(model == 'simple'){
    ageModelNm = ageModel_simp
    jagsParams = 'pi'
  } else {
    ageModelNm = ageModel_hier
    jagsParams = c('pi', 'mu', 'Tau', 'avgPi')
  }

for(yr in year_range){  
  cat(paste('Starting age model for', spp, 'in run year', yr, '\n'))
  # read in data
  modAgeDf = read_excel(paste0(filepath, 'LGR_', spp, '_', yr, '.xlsx'),
                        'AgeFreq',
                        progress = F)
  
  # pull out relevant bits for JAGS, and name them appropriately
  age_jagsData = modAgeDf %>%
    filter(!is.na(TRT)) %>%
    mutate(popNum = as.integer(as.factor(TRT))) %>%
    select(tags = nAged,
           popNum) %>%
    as.list()
  
  age_jagsData$ageMat = modAgeDf %>%
    filter(!is.na(TRT)) %>%
    select(starts_with('age')) %>%
    as.matrix()
  
  # drop ages with no observed fish in them
  if(sum(colSums(age_jagsData$ageMat) == 0) > 0) {
    age_jagsData$ageMat = age_jagsData$ageMat[,!colSums(age_jagsData$ageMat) == 0]
  }
  
  # gather date specific to hierarchical model
  if(model == 'hierarchical'){
  
  # assign populations to groups (ex: A-run or B-run)
  age_jagsData$runType = modAgeDf %>%
    filter(Group != 'NotObserved') %>%
    select(TRT) %>%
    mutate(Run = if_else(TRT %in% c('CRLMA-s',
                                    'CRLOC-s',
                                    'CRLOL-s',
                                    'CRSEL-s',
                                    'CRSFC-s',
                                    'MFBIG-s',
                                    'SFMAI-s',
                                    'SFSEC-s'),
                         'B', 'A')) %>%
    distinct() %>%
    mutate(Run = as.factor(Run),
           runType = as.integer(Run)) %>%
    pull(runType)
  
  age_jagsData$R = diag(1, ncol(age_jagsData$ageMat) - 1)
  
  }

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
