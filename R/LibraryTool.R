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



#' @title GenerateMSP
#' @author Zhiwei Zhou
#' @description Convert MS/MS spectra to MSP files
#' @param file_name Required. The file name of MSP. The suffix of file name must be ".msp".
#' @param cmp_name Required
#' @param precusormz Required
#' @param spec Required. It should be in a dataframe, 1st: "mz", 2nd: "intensity"
#' @param adduct Default: NULL
#' @param instrument_type Default: NULL
#' @param instrument Default: NULL
#' @param smiles Default: NULL
#' @param inchikey Default: NULL
#' @param inchikey1 Default: NULL
#' @param formula Default: NULL
#' @param polarity 'Positive' or 'Negative'
#' @param ce Default: NULL
#' @param rt Default: NULL
#' @param ccs Default: NULL
#' @param zhulib_id Default: NULL
#' @param kegg_id Default: NULL
#' @param hmdb_id Default: NULL
#' @param pubchem_id Default: NULL
#' @param links Default: ''
#' @param comment Default: ''
#' @export
#' @example
#' setwd('F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/')
#'     GenerateMSP(file_name = 'zhumetlib_validation_pos_20v_190520.msp',
#'                 cmp_name = external_data_pos$compound_name[i],
#'                 precusormz = external_data_pos$mz[i],
#'                 adduct = external_data_pos$adducts[i],
#'                 instrument_type = 'LC-ESI-qTof',
#'                 instrument = 'LC-ESI-qTof',
#'                 smiles = external_data_pos$smiles[i],
#'                 inchikey = external_data_pos$inchikey[i],
#'                 inchikey1 = external_data_pos$inchikey1[i],
#'                 formula = external_data_pos$formula[i],
#'                 polarity = 'Positive',
#'                 ce = '20',
#'                 ccs = external_data_pos$CCS[i],
#'                 zhulib_id = external_data_pos$matched_zhulib_id[i],
#'                 pubchem_id = external_data_pos$pubchem_cid[i],
#'                 comment = paste(external_data_pos$id[i], external_data_pos$source[i], sep = ' '),
#'                 spec = temp_spec)



setGeneric(name = 'GenerateMSP',
           def = function(
             file_name = "./test.msp",
             cmp_name = NULL,
             precusormz = NULL,
             adduct=NULL,
             instrument_type=NULL,
             instrument=NULL,
             smiles=NULL,
             inchikey=NULL,
             inchikey1=NULL,
             formula=NULL,
             polarity=c('Positive', 'Negative'),
             ce=NULL,
             rt=NULL,
             ccs=NULL,
             zhulib_id=NULL,
             kegg_id=NULL,
             hmdb_id=NULL,
             pubchem_id=NULL,
             links='',
             comment='',
             spec=NULL
           ){

             if (!stringr::str_detect(file_name, '.msp')) {
               stop('The suffix of file name must be .msp\n')
             }

             if (is.null(cmp_name)) {
               stop('Please input cmp_name\n')
             }

             if (is.null(precusormz)) {
               stop('Please input precusormz\n')
             }

             if (is.null(spec)) {
               stop('Please input spec\n')
             }

             polarity = match.arg(polarity)


             if (ncol(spec)!=2) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }

             if (!all(colnames(spec)==c('mz', 'intensity'))) {
               stop('Please check the format of spectrum. It should be in a dataframe, 1st: "mz", 2nd: "intensity" \n')
             }



             # write into msp
             file_result <- file(description = file_name, open = "a")
             cat('NAME: ', cmp_name, '\n', sep = '', file = file_result)
             cat('PRECURSORMZ: ', precusormz, '\n', sep = '', file = file_result)

             if (!is.null(adduct)) {
               cat('PRECURSORTYPE: ', adduct, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument_type)) {
               cat('INSTRUMENTTYPE: ', instrument_type, '\n', sep = '', file = file_result)
             }

             if (!is.null(instrument)) {
               cat('INSTRUMENT: ', instrument, '\n', sep = '', file = file_result)
             }

             if (!is.null(smiles)) {
               cat('SMILES: ', smiles, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey)) {
               cat('INCHIKEY: ', inchikey, '\n', sep = '', file = file_result)
             }

             if (!is.null(inchikey1)) {
               cat('INCHIKEY1: ', inchikey1, '\n', sep = '', file = file_result)
             }

             if (!is.null(formula)) {
               cat('FORMULA: ', formula, '\n', sep = '', file = file_result)
             }

             if (!is.null(polarity)) {
               cat('IONMODE: ', polarity, '\n', sep = '', file = file_result)
             }

             if (!is.null(ce)) {
               cat('COLLISIONENERGY: ', ce, '\n', sep = '', file = file_result)
             }

             if (!is.null(rt)) {
               cat('RETENTIONTIME: ', rt, '\n', sep = '', file = file_result)
             }

             if (!is.null(ccs)) {
               cat('COLLISIONCROSSSECTION: ', ccs, '\n', sep = '', file = file_result)
             }

             if (!is.null(zhulib_id)) {
               cat('ZHULAB: ', zhulib_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(kegg_id)) {
               cat('KEGG: ', kegg_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(hmdb_id)) {
               cat('HMDB: ', hmdb_id, '\n', sep = '', file = file_result)
             }

             if (!is.null(pubchem_id)) {
               cat('PUBCHEM: ', pubchem_id, '\n', sep = '', file = file_result)
             }

             cat('Links: ', links, '\n', sep = '', file = file_result)
             cat('Comment: ', comment, "\n", sep = '', file = file_result)

             cat('Num Peaks: ',  nrow(spec),  '\n',  sep = '', file = file_result)

             for (i in 1:nrow(spec)) {
               cat(paste(as.numeric(round(spec[i,1], digits = 4)),
                         as.numeric(round(spec[i,2], digits = 2)),
                         collapse = ' '),
                   '\n', sep = '', file = file_result)
             }

             cat('\n', file = file_result)

             close(file_result)
           }
)
