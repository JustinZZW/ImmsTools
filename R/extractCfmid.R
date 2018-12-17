#' @title extractCfmid
#' @description extract cfm_id predicted spec, .log
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param file_name the name of log files. e.g. C00182.log
#' @param mode 'spec_only' or 'spec_annotate'. spec_only: extract spec; spec_annotate: extract spec with annotate
#' @export
#' @examples
#' extractCfmid("C00182.log", mode = 'spec_annotate')

extractCfmid <- function(file_name,
                         mode = c('spec_only', 'spec_annotate')) {
  mode <- match.arg(mode)
  raw_data <- readr::read_lines(file_name)
  idx_energy0_start <- which(raw_data=="energy0") + 1
  idx_energy0_end <- which(raw_data=="energy1") - 1

  idx_energy1_start <- which(raw_data=='energy1') + 1
  idx_energy1_end <- which(raw_data=="energy2") - 1

  idx_energy2_start <- which(raw_data=='energy2') + 1
  idx_energy2_end <- which(raw_data=="") - 1

  idx_annotate_start <- which(raw_data=='')+1
  idx_annotate_end <- length(raw_data)

  spec_energy0 <- raw_data[idx_energy0_start:idx_energy0_end]
  spec_energy1 <- raw_data[idx_energy1_start:idx_energy1_end]
  spec_energy2 <- raw_data[idx_energy2_start:idx_energy2_end]
  info_annotate <- raw_data[idx_annotate_start:idx_annotate_end]

  spec_energy0 <- separate_spec(spec_energy0)
  spec_energy1 <- separate_spec(spec_energy1)
  spec_energy2 <- separate_spec(spec_energy2)

  info_annotate <- separate_info(info_annotate)

  if (mode=='spec_only') {
    result <- list(spec_energy0[, 1:2, drop=F],
                   spec_energy1[, 1:2, drop=F],
                   spec_energy2[, 1:2, drop=F])

    names(result) <- c(10, 20, 40)
  } else {
    result <- list(spec_energy0,
                   spec_energy1,
                   spec_energy2,
                   info_annotate)

    names(result) <- c(10, 20, 40, 'info_annotate')
  }

  return(result)
}


separate_spec <- function(raw_spec){
  raw_spec <- stringr::str_split(string = raw_spec, pattern = ' ')
  raw_spec <- do.call(rbind, raw_spec)[,1:3,drop=F]
  # spec_result <- raw_spec
  if (nrow(raw_spec) > 1) {
    spec_result <- apply(raw_spec, 2, function(x){
      as.numeric(x)
    })
  } else {
    spec_result <- apply(raw_spec, 2, function(x){
      as.numeric(x)
    })
    spec_result <- t(as.matrix(spec_result))
  }

  colnames(spec_result) <- c('mz', 'intensity', 'annotation_idx')

  spec_result
}


separate_info <- function(raw_info){
  raw_info <- stringr::str_split(string = raw_info, pattern = ' ')
  raw_info <- do.call(rbind, raw_info)
  raw_info <- data.frame(number=as.numeric(raw_info[,1]),
                         mz=as.numeric(raw_info[,2]),
                         annotation=as.character(raw_info[,3]),
                         stringsAsFactors = F)
}
