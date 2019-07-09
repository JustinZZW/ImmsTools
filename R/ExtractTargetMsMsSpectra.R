#' @title ExtractTargetMsMs
#' @author Zhiwei Zhou
#' @param msms_file mzXML file format
#' @param polarity 'positive', 'negative'. Default: 'positive'
#' @param mode SPMCE: single precursor with multiple CEs; MPSCE: multiple precursors, but single CE for each precursor; SPSCE: single precursor with one CE;
#' @param ce_levels collison energy. Default: c("CE10", "CE20", "CE40")
#' @export
#' @examples
#' setwd('F:/01 MetIMMS/00 data processing/190702 extract msms spectra from targeted 6560 data/convert_msms')
#' test <- ExtractTargetMsMS(msms_file = "L0656_H_pos_10ppm.mzXML", polarity = 'positive', mode = 'SPMCE', ce_levels = c("CE10", "CE20", "CE40"))

setGeneric(name = 'ExtractTargetMsMS',
           def = function(
             msms_file,
             polarity=c('positive', 'negative'),
             mode = c('SPMCE', 'MPSCE', 'SPSCE'),
             ce_levels = c("CE10", "CE20", "CE40")
           ){
             polarity <- tolower(polarity)
             polarity <- match.arg(polarity)

             # extract ms1 scan spectrum ---------------------------------------
             # ce.num = length(ce.levels)
             xraw <- xcms::xcmsRaw(msms_file, includeMSn = TRUE)
             scantime <- xraw@scantime
             tic <- xraw@tic
             idx_tic_max <- which.max(tic)

             mz_precursors <- unique(xraw@msnPrecursorMz)

             # check parameters
             switch(mode,
                    'SPMCE' = {
                      if(length(mz_precursors)!=1) {
                        stop('It has multiple precursors, please check data\n')
                      }

                      if (length(ce_levels)==1) {
                        stop('Please check parameters, because only one ce_level was input\n')
                      }

                    },
                    'MPMCE' = {
                      if (length(mz_precursors)==1) {
                        stop('It only has one precursor, please check data\n')
                      }

                      if (length(ce_levels)==1) {
                        stop('Please check parameters, because only one ce_level was input\n')
                      }
                    },
                    'SPSCE' = {
                      if (length(mz_precursors)!=1) {
                        stop('It has multiple precursors, please check data\n')
                      }

                      if (length(ce_levels)!=1) {
                        stop('Please check parameters, because multiple ce_levels were input\n')
                      }
                    }

             )

             apex_ms1_scantime <- scantime[idx_tic_max]
             apex_ms1_tic <- tic[idx_tic_max]
             apex_ms1_scan <- xcms::getScan(xraw, idx_tic_max)

             # check ms1 quality
             ms1_quality <- lapply(mz_precursors, function(x){
               idx_ms1 <- which(apex_ms1_scan[,1] >= x-0.01 & apex_ms1_scan[,1] <= x + 0.01)
               if (length(idx_ms1)==0) {
                 ms1_check <- FALSE
                 measured_precursor <- NA
               } else {
                 ms1_check <- max(apex_ms1_scan[idx_ms1,2]) >= 1000

                 if (length(idx_ms1)>1) {
                   idx_ms1 <- idx_ms1[which.max(as.numeric(apex_ms1_scan[idx_ms1,2]))]
                   measured_precursor <- as.numeric(apex_ms1_scan[idx_ms1,1])
                 }

                 measured_precursor <- as.numeric(apex_ms1_scan[idx_ms1,1])

               }

               result <- list(ms1_check, measured_precursor)
             })

             ms1_check <- sapply(ms1_quality, function(x)x[[1]])
             measured_precursor <- sapply(ms1_quality, function(x)x[[2]])



             if (!all(ms1_check)) {
               info <- data.frame(mz_precursors=mz_precursors,
                                  measured_precursor=measured_precursor,
                                  polarity=polarity,
                                  ms1_check=ms1_check)
               # spec <- NA
               result <- list(ms1=info,
                              ms2=NA)
               return(result)
             }


             # extract ms2 spectra ---------------------------------------------
             xraw_ms2 <- xcms::msn2xcmsRaw(xraw)

             scantime_ms2 <- xraw_ms2@scantime

             # compare experiment RT with expect RT to avoid wrong scan order caused by invalid ms2 scan
             # experimental ms2 scan RT
             scantime_ms2_experiment <- scantime_ms2[scantime_ms2 > scantime[idx_tic_max] &
                                                       scantime_ms2 < scantime[idx_tic_max+1]
                                                     ]
             # predicted ms2 scan RT
             scantime_ms2_predict <- seq(scantime[idx_tic_max], scantime[idx_tic_max+1],
                                         length.out = c(length(ce_levels)+2))
             scantime_ms2_predict <- scantime_ms2_predict[-c(1, length(scantime_ms2_predict))]
             # names(scantime_ms2_predict) <- paste(ce_levels)

             idx_scan_ms2 <- GetIdxMsMsScan(scantime_ms2_experiment = scantime_ms2_experiment,
                                            scantime_ms2_predict = scantime_ms2_predict,
                                            xraw_ms2 = xraw_ms2)

             apex_ms2_result <- lapply(idx_scan_ms2, function(x) {
               # if(is.na(x)) return(NA)
               spec <- xcms::getScan(xraw_ms2, x)
               sn <- max(spec[, 2])/median(spec[, 2])

               precursor <- xraw_ms2@msnPrecursorMz[x]
               ce <- xraw_ms2@msnCollisionEnergy[x]

               purified_spec <- MetIMMS102::purifyMsMs(spec = spec,
                                                       mz_precursor = precursor,
                                                       is_include_precursor = TRUE,
                                                       mz_range_ms2 = c(25, 1700),
                                                       int_ms2_min_abs = 50,
                                                       int_ms2_min_relative = 0.01,
                                                       ppm_precursor_filter = 10)

               purified_spec[,2] <- as.numeric(purified_spec[, 2]/max(purified_spec[,2])*100)

               ms2_info <- data.frame(precursor=precursor,
                                      ce=ce,
                                      sn=sn,
                                      stringsAsFactors = F)

               result <- list(ms2_info=ms2_info,
                              ms2_spec=purified_spec)
             })

             # output result ---------------------------------------------------
             # info <- c(name, mz_precursor, concentrition, polarity, ms1_check)
             # spec <- NA
             info <- data.frame(mz_precursors=mz_precursors,
                                measured_precursor=measured_precursor,
                                polarity=polarity,
                                ms1_check=ms1_check)

             result <- list(ms1=info,
                            ms2=apex_ms2_result)

             return(result)
           }
)


# return ms2 scan idx according the error between experimental RT and predicted RT ------------
GetIdxMsMsScan <- function (scantime_ms2_experiment, scantime_ms2_predict, xraw_ms2)
{
  idx_cor_scantime_experiment <- sapply(scantime_ms2_predict, function(x){
    idx <- which.min(abs(x-scantime_ms2_experiment))
  })

  scantime_ms2_experiment <- scantime_ms2_experiment[idx_cor_scantime_experiment]

  idx_msms_scan <- lapply(scantime_ms2_experiment,
                          function(x) {idx <- which(xraw_ms2@scantime %in% x)})
  names(idx_msms_scan) <- names(scantime_ms2_experiment)
  return(idx_msms_scan)
}


#' @title TargetMsMsProcessing
#' @author Zhiwei Zhou
#' @param msms_file mzXML file format
#' @param root_path Default: '.'
#' @param polarity 'positive', 'negative'. If it is not given (NULL), it would be extracted from msms_file name. Default: NULL
#' @param adduct If it is not given (NULL), it would be extracted from msms_file name. Default: NULL
#' @param concentration If it is not given (NULL), it would be extracted from msms_file name. Default: NULL
#' @param mode SPMCE: single precursor with multiple CEs; MPSCE: multiple precursors, but single CE for each precursor; SPSCE: single precursor with one CE;
#' @param ce_levels collison energy. Default: c("CE10", "CE20", "CE40")
#' @export
#' @examples
#' setwd('F:/01 MetIMMS/00 data processing/190702 extract msms spectra from targeted 6560 data/convert_msms')
#' # Standard
#' TargetMsMsProcessing(msms_file = "L0656_H_pos_10ppm.mzXML", root_path = '..', mode = 'SPMCE', ce_levels = c("CE10", "CE20", "CE40"))
#' # QC
#' TargetMsMsProcessing(msms_file = "QC00.mzXML", root_path = '..', polarity = 'Positive', adduct = '[M+H]+', concentration = '10ppm',mode = 'MPMCE', ce_levels = c("CE20", "CE40"))


setGeneric(name = 'TargetMsMsProcessing',
           def = function(
             msms_file = "L0656_H_pos_10ppm.mzXML",
             root_path = '.',
             polarity = NULL,
             adduct = NULL,
             concentration = NULL,
             mode = c('SPMCE', 'MPSCE', 'SPSCE'),
             ce_levels = c("CE10", "CE20", "CE40")
           ){
             name <- stringr::str_split(msms_file, pattern = '\\.mzXML')[[1]][1]

             if (is.null(polarity)) {
               name <- stringr::str_split(msms_file, pattern = '_')[[1]][1]
               adduct <- stringr::str_split(msms_file, pattern = '_')[[1]][2]
               polarity <- stringr::str_split(msms_file, pattern = '_')[[1]][3]
               concentration <- stringr::str_extract(string = stringr::str_split(msms_file, pattern = '_')[[1]][4],
                                                     pattern = '\\d+ppm')

               switch (polarity,
                       'pos' = {
                         polarity <- 'Positive'
                         ifelse(adduct=='H', adduct <- '[M+H]+', adduct <- '[M+Na]+')
                       },
                       'neg' = {
                         polarity <- 'Negative'
                         ifelse(adduct=='H', adduct <- '[M-H]-', adduct <- '[M+Na-2H]-')
                       }
               )
             }



             extract_result <- ExtractTargetMsMS(msms_file = msms_file,
                                                 polarity = polarity,
                                                 mode = mode,
                                                 ce_levels = ce_levels)
             spec_info <- extract_result[[1]]
             spec_result <- extract_result[[2]]

             # if it meet the ms1 criteria, export ms/ms spectra
             if (sum(spec_info$ms1_check)>0) {
               output_path <- file.path(root_path,
                                        paste(name, polarity, concentration, sep = '_'))

               dir.create(output_path,
                          recursive = TRUE)

               temp <- lapply(spec_result, function(x){
                 temp_info <- x[[1]]
                 temp_spec <- x[[2]]
                 ImmsTools::GenerateMSP(file_name = file.path(output_path,
                                                              paste0(name, '_',
                                                                     polarity, '_',
                                                                     concentration, '.msp')),
                                        cmp_name = name,
                                        precusormz = temp_info$precursor,
                                        adduct = adduct,
                                        instrument_type = 'LC-IM-QTOF',
                                        instrument = 'Agilent6560',
                                        polarity = polarity,
                                        ce = temp_info$ce,
                                        zhulib_id = name,
                                        comment = paste0('sn_', temp_info$sn),
                                        spec = temp_spec)
               })
             }

             spec_info <- spec_info %>%
               dplyr::mutate(name=name, adduct=adduct)

             return(spec_info)
           }
)
