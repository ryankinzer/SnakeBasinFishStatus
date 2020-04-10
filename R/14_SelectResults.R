# Author: Kevin See
# Purpose: pull out select results
# Created: 3/24/20
# Last Modified: 3/31/20
# Notes: Marika Dobos originally made this request (IDFG). 
# On 3/30, Ian Tattam (ODFW) asked for UGR results

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
library(readxl)
library(WriteXLS)

#-----------------------------------------------------------------
# which sites do we want to pull out?
my_sites = c("BBA",
             "KHS",
             "LRL",
             "LRU",
             "TAY",
             "LLR",
             "NFS",
             "SFG",
             "KRS",
             "LRW",
             "HYC")

my_sites = c('UGR')

#-----------------------------------------------------------------
excel_sheets('Abundance_results/LGR_AllSummaries_Steelhead.xlsx')

# pull out selected results
res_list = list('Pop Total Esc' = read_excel('Abundance_results/LGR_AllSummaries_Steelhead.xlsx',
                                             sheet = 'Pop Total Esc') %>%
                  filter(valid_est == 1) %>%
                  filter(grepl('Grande Ronde', MPG)) %>%
                  rename(estimate = median) %>%
                  mutate_at(vars(cv),
                            list(as.numeric)) %>%
                  select(-valid_est, -mean, -mode),
                'Site Esc' = read_excel('Abundance_results/LGR_AllSummaries_Steelhead.xlsx',
                                        sheet = 'Site Esc') %>%
                  mutate_at(vars(cv),
                            list(as.numeric)) %>%
                  filter(site %in% my_sites),
                'Node Detect Eff' = read_excel('Abundance_results/LGR_AllSummaries_Steelhead.xlsx',
                                               sheet = 'Node Detect Eff') %>%
                  mutate_at(vars(cv),
                            list(as.numeric)) %>%
                  mutate(site = str_remove(Node, 'A0$'),
                         site = str_remove(site, 'B0$')) %>%
                  filter(site %in% my_sites) %>%
                  select(-site))

WriteXLS(x = res_list,
         # ExcelFileName = paste0('outgoing/Sthd_Dobos_', format(Sys.Date(), '%Y%m%d'), '.xlsx'),
         ExcelFileName = paste0('outgoing/Sthd_GrandeRonde_', format(Sys.Date(), '%Y%m%d'), '.xlsx'),
         AdjWidth = T,
         BoldHeaderRow = T,
         FreezeRow = 1)
