#------------------------------------------------------------------------------
# Open connection to LGTrapping database and save table as .csv.
#
# Author: Ryan Kinzer
#------------------------------------------------------------------------------

library(tidyverse)

## Requires a copy of the Lower Granite Dam Trap Database and the odbc driver
if(.Platform$OS.type != 'unix') {
  source('./R/loadLGTrappingDBase.R')
  trap_filepath <- './data/TrappingDBase/LGTrappingExportJodyW.accdb'
  con <- loadLGTrappingDBase(trap_filepath)
  trap_dbase <- DBI::dbReadTable(con, 'tblLGDMasterCombineExportJodyW')
  DBI::dbDisconnect(con)
}

# this works for Mac, or any system really
if(.Platform$OS.type == 'unix') {
  # trap_filepath <- './data/TrappingDBase/LGTrappingExportJodyW.accdb'
  trap_filepath <- './data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv'
  trap_dbase = readLGRtrapDB(trap_filepath)
}

# save .csv of dbase for later use
write_csv(trap_dbase, file = './data/TrappingDBase/tblLGDMasterCombineExportJodyW.csv')
