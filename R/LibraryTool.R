#' @title SplitInchiKey
#' @author Zhiwei Zhou
#' @description split inchikey to 3 parts
#' @return a data.frame with 4 columns
#' @export
#' @example
#' SplitInchiKey(inchikey = 'VGONTNSXDCQUGY-RRKCRQDMSA-N')

setGeneric(name = 'SplitInchiKey',
           def = function(
             inchikey = 'AEMRFAOFKBGASW-UHFFFAOYSA-N'
           ){
             result <- stringr::str_split_fixed(inchikey, pattern = '-', n = 3)
             result <- data.frame(result, stringsAsFactors = F)
             colnames(result) <- c('inchikey1', 'inchikey2', 'inchikey3')

             result <- data.frame(result,
                                  inchikey=inchikey,
                                  stringsAsFactors = F)

             return(result)
           }
)



#' @title GetMsMsInfo
#' @author Zhiwei Zhou
#' @description extract compound information with zhulib
#' @param database a database with zhulib info
#' @param id a character of lib id
#' @return a data.frame with compound info
#' @export
#' @example
#' load('F:/01 MetIMMS/00 data processing/190226 intact NIST and Metlin database/metlin/metlinlib.rdata')
#' GetMsMsInfo(database = metlinlib, id = 'M0001')

setGeneric(name = 'GetMsMsInfo',
           def = function(
             database,
             id='L0002'
             # polarity=c('pos', 'neg'),
             # instrument=c('ABSciex', 'Agilent'),
             # ce=c('10', '15', '20', '25', '30', '35', '35,15', '40', '45', '50', '55', '60', '65', '70')
           ){
             database_info <- database[[1]][[2]]$compound

             id_list <- database_info$labid
             id_idx <- which(id_list==id)

             if (length(id_idx)==0) {
               result <- rep(NA, ncol(database_info))
               result <- t(data.frame(result))
               colnames(result) <- colnames(database_info)
               rownames(result) <- NULL

               return(result)
             }

             result <- database_info[id_idx,,drop=F]
             rownames(result) <- NULL

             return(result)

           }
)




#' @title GetMsMs
#' @author Zhiwei Zhou
#' @description extract compound MS/MS spectra from zhulib
#' @param database a database with zhulib info
#' @param id a character of lib id
#' @param polarity 'pos', 'neg'. Default: 'pos'
#' @param instrument 'ABSciex', 'Agilent'
#' @param ce character. '10', '15', '20', '25', '30', '35', '35,15', '40', '45', '50', '55', '60', '65', '70'
#' @return a data.frame with compound info
#' @export
#' @example
#' load('F:/01 MetIMMS/00 data processing/190226 intact NIST and Metlin database/metlin/metlinlib.rdata')
#' GetMsMs(database = metlinlib, polarity = 'pos', id = 'M0001', instrument = 'Agilent', ce = '20')

setGeneric(name = 'GetMsMs',
           def = function(
             database,
             polarity=c('pos', 'neg'),
             id='L0002',
             instrument=c('ABSciex', 'Agilent'),
             ce=c('10', '15', '20', '25', '30', '35', '35,15', '40', '45', '50', '55', '60', '65', '70')
           ){
             polarity <- match.arg(polarity)
             instrument <- match.arg(instrument)
             ce <- match.arg(ce)

             database_msms <- database[[1]][[1]]

             switch(polarity,
                    'pos' = {database_msms <- database_msms[[1]]},
                    'neg' = {database_msms <- database_msms[[2]]}
             )

             # convert CEs
             switch (instrument,
                     'ABSciex' = {ce <- as.numeric(ce)},
                     'Agilent' = {
                       if (ce=='35,15') {
                         stop('Agilent data do not have ce 35,15')
                       }

                       ce <- as.numeric(ce)+10
                     }
             )

             # extract the msms data of input id
             #  if no id, return a dataframe, mz=NA, intensity=NA
             id_list <- names(database_msms)
             id_idx <- which(id_list==id)

             if (length(id_idx)==0) {
               result <- data.frame(mz=NA, intensity=NA, stringsAsFactors = F)
               return(result)
             }

             temp_msms <- database_msms[[id_idx]]


             # extract the msms data of input ce
             #  if no id, return a dataframe, mz=NA, intensity=NA
             ce_list <- names(temp_msms)

             ce_idx <- which(ce_list==as.character(ce))

             if (length(ce_idx)==0) {
               result <- data.frame(mz=NA, intensity=NA, stringsAsFactors = F)
               return(result)
             }

             result <- temp_msms[[ce_idx]]

             return(result)
           }
)
