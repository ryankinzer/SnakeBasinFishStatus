#' @title Define Popuation Groups - LGD
#'
#' @description Define which main branches or detection sites fall into TRT populations.
#'
#' @author Ryan Kinzer
#'
#' @param spp Species, either "Chinook" or "Steelhead"
#' @param node_order output of function \code{createNodeOrder}.
#'
#' @import dplyr stringr
#' @export
#' @return NULL
#' @examples definePopulations()

definePopulations = function(spp = c('Chinook', 'Steelhead')) {

  spp = match.arg(spp)

  if(spp == 'Chinook') {
## NEEDS FIXED
    # report_df = list('CRLAP' = 'LAP',
    #                  'CRPOT' = 'JUL',
    #                  'SCUMA' = 'SC1',
    #                  'IRBSH' = c('CMP', 'BSC', 'LSHEEF'),
    #                  'IRMAI' = c('IR1', 'COC'),
    #                  'GRLOS' = 'WR1', # WR1 includes Minam
    #                  'SNASO' = c('ACM', 'PENAWC', 'ALMOTC', 'ALPOWC', 'TENMC2'),
    #                  'SNTUC' = 'LTR',
    #                  'MFBIG' = 'TAY',
    #                  'SFEFS' = 'ESS',
    #                  'SFSEC' = 'ZEN',
    #                  'SFSMA' = c('KRS', 'SFG_bb'),
    #                  'SRLEM' = c('LLR','CRC'),
    #                  'SRNFS' = 'NFS',
    #                  'SRLMA' = 'USI_bb',
    #                  'SRUMA' = 'STL',
    #                  'SRVAL' = 'VC2',
    #                  'SRYFS' = 'YFK',
    #                  'CRLOC' = 'LRL',
    #                  'CRLOL' = 'LC1',
    #                  'SEUMA' = 'SW1',
    #                  'GRCAT' = 'CATHEW',
    #                  'GRUMA' = 'UGR', # need to subtract out CATHEW
    #                  'GRLOO' = 'LOOKGC',
    #                  'SRLSR' = 'RAPH',
    #                  'SREFS' = 'SALEFT',
    #                  'SRPAH' = 'PAHH') %>%
    #   stack() %>%
    #   tbl_df() %>%
    #   select(TRT_POPID = ind,
    #          SiteID = values)
  }

  if(spp == 'Steelhead') {

    report_df = list('CRLOC-s' = 'Lochsa',
                     'CRLOL-s' = 'Lolo',
                     'CRSEL-s' = 'Selway',
                     'CRSFC-s' = 'SFClearwater',
                     'GRUMA-s' = c('GrandeRonde', 'LookingGlass'),
                     'GRJOS-s' = 'JosephCreek',
                     'GRWAL-s' = 'Wallowa',
                     'IRMAI-s' = c('ImnahaRiver', 'CowCreek'),
                     'SNASO-s' = c('Asotin', 'Penawawa', 'Almota', 'Alpowa', 'TenMileCreek'),
                     'SREFS-s' = c('SALEFT'), # need to think about USI_bb
                     'SRLEM-s' = c('Lemhi','CarmenCreek'),
                     'SRNFS-s' = 'NFSalmon',
                     'SRPAH-s' = 'PAHH',
                     'SRUMA-s' = c('past_YFK', 'past_VC2', 'past_STL'),
                     'SFSEC-s' = 'past_ZEN',
                     'SFMAI-s' = c('past_KRS', 'past_ESS', 'SFSalmon_bb'),
                     'CRLMA-s' = c('Potlatch', 'Lapwai', 'ClearCreek'),
                     'SNTUC-s' = 'Tucannnon',
                     'MFBIG-s' = 'BigCreek',
                     'SRLSR-s' = 'RapidRiver') %>%
      stack() %>%
      tbl_df() %>%
      select(TRT = ind,
             area = values)

  }

  return(report_df)

}

