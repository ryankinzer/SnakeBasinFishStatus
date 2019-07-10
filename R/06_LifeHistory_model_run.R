# Summarize life history data from cleaned capture history files output from
# PITcleanR.


# Load Packages
library(tidyverse)
library(readlx)

 filepath <- './data/CleanedProcHist'
# files <- list.files(filepath)
# n_files <- length(files)

# Set species and Spawn Year
spp = 'Steelhead'
yr = 2010

# read in data
proc_ch = read_excel(paste0(filepath,'/','LGR_',spp,'_EDITTED_',yr,'.xlsx'))




