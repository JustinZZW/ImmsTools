#' @title readMSP
#' @description read MSP spectra files
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param file the file name
#' @param mode standard: extract name, mz and RT from the file, which fit for MSP data exported from QI software; all: extract all information in the file
#' @return
#' A list object. info: MS1 information, including ms1, rt etc.; spec: the MS/MS spectrum
#' @export
#' @examples
#' test <- readMSP(file = 'F:/01 MetIMMS/00 data processing/190515 external validation msms data extraction/zhumetlib_validation_pos_20v_190520.msp', mode = 'all')

setGeneric('readMSP', function(file,
                               mode=c('standard', 'all')) {
  # devtools::use_package('dplyr')

  mode <- match.arg(mode)
  msp.data.list <- ListDB(file)
  nr.num.pk <- grep('Num Peaks', stringr::str_to_title(msp.data.list[[1]]))
  info.spec <- lapply(msp.data.list, function(msp.data) {
    info.list <- strsplit(msp.data[1:nr.num.pk], split = ': ')
    info <- do.call(cbind, info.list)
    colnames(info) <- info[1, ]

    if (mode=='standard') {
      mz <- round(as.numeric(info[2,3]), digits = 4)
      rt <- round(as.numeric(strsplit(info[2, "Comment"], split = "_")[[1]][1])*60, digits = 0)
      name <- info[2, "Comment"]
      # info <- matrix(c(mz, rt), ncol = 2)
      info <- data.frame(mz=mz, rt=rt)
      rownames(info) <- name
      # colnames(info) <- c("mz", "rt")
    } else {
      info <- as.data.frame(tibble::as.tibble(info))
      info <- info[-1,,drop=F]
      rownames(info) <- NULL
    }

    spec.list <- strsplit(msp.data[(nr.num.pk+1):length(msp.data)], ' |\\t')
    spec.list <- lapply(spec.list, as.numeric)
    spec <- do.call(rbind, spec.list)[, c(1, 2), drop = FALSE]
    colnames(spec) <- c('mz', 'intensity')
    # spec <- list('spec' = spec)

    # list('info' = info[-1, , drop = FALSE],
    #      'spec' = spec)

    list('info' = info,
         'spec' = spec)

  })

  return(info.spec)
})

setGeneric('ListDB', function(file) {
  msp.data <- readLines(file)
  nl.db.new <- 1
  idx.db <- 1
  db.list <- list()
  len.data <- length(msp.data)
  for(nl in 1:len.data)
  {
    if(msp.data[nl]=="") {
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:(nl-1)])
      nl.db.new <- nl + 1
      idx.db <- idx.db + 1
    } else if (nl == len.data){
      db.list[idx.db] <- list(Compound=msp.data[nl.db.new:nl])
    }
  }
  db.list
})



#' @title mz_transform
#' @description calculate the m/z of adducts
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @usage
#' \code{mz_transform(M, adduct, polarity=c('pos', 'neg'))}
#' @param M a vector of exact mass.
#' @param adduct a vector of adducts. These adducts need to be
#' @param polarity 'pos' or 'neg'
#' @param formula sdf
#' @export
#' @examples
#' mz_transform(M = 123.2324, adduct = c('[M+H]+', '[M+Na]+'), polarity = 'pos')
#' mz_transform(formula="C2H5OH", adduct = c('M-', '[M-H]-'), polarity = 'neg')


# transform exact mass to mz according adduct form
mz_transform <- function(M=NULL,
                         adduct=NULL,
                         formula=NULL,
                         polarity=c("pos", 'neg')){

  polarity <- match.arg(polarity)

  if (all(is.null(M), is.null(formula))) {
    stop('Please input M or formula.')
  }


  if (!is.null(M)) {
    M <- as.numeric(M)
  } else {
    M <- Calcu_EM(formula)
  }


  if (polarity=="pos") {

    if (is.null(adduct)) {
      adduct <- adduct.table$pos$adduct
    } else {
      temp <- adduct %in% adduct.table$pos$adduct

      if (sum(temp)!=length(adduct)) {
        stop("Please check the adduct list.")
      }
    }

    adduct.table <- adduct.table$pos

  } else {
    if (is.null(adduct)) {
      adduct <- adduct.table$neg$adduct
    } else {
      temp <- adduct %in% adduct.table$neg$adduct

      if (sum(temp)!=length(adduct)) {
        stop("Please check the adduct list.")
      }
    }

    adduct.table <- adduct.table$neg

  }

  idx.adduct <- match(adduct, adduct.table$adduct)

  mz <- sapply(c(1:length(idx.adduct)), function(i){
    temp.idx <- idx.adduct[i]
    mz <- M + adduct.table$mz[temp.idx]
  })

  output <- data.frame(M=M,
                       adduct=adduct,
                       mz=mz,
                       stringsAsFactors = F)
  return(output)
}


#' @title Calcu_EM
#' @description Calculate exact mass using formula
#' @author Zhiwei Zhou
#' \email {zhouzw@@sioc.ac.cn}
#' @param formula
#' @example
#' Calcu_EM("C2H5OH")
#' @export

Calcu_EM <- function(formula="C2H5OH") {
  molecule <- Rdisop::getMolecule(formula)
  # getFormula(molecule)
  Rdisop::getMass(molecule)
}



#' @title impute
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @param raw_matrix the numeric matrix
#' @param label the label of each row
#' @param k please ref impute.knn. Default: 10
#' @param rowmax please ref impute.knn. Default: 0.5
#' @param colmax please ref impute.knn. Default: 0.8
#' @param maxp please ref impute.knn. Default: 1500
#' @param rng.seed please ref impute.knn. Default: 362436069
#' @export

ZZWImpute <- function(raw_matrix,
                      label,
                      k=10,
                      rowmax=0.5,
                      colmax=0.8,
                      maxp=1500,
                      rng.seed=362436069
){
  n_inf <- apply(raw_matrix, 2, function(x){
    sum(is.infinite(x))
  })

  n_na <- apply(raw_matrix, 2, function(x){
    temp <- is.na(x)
    sum(temp)
  })

  n_nan <- apply(raw_matrix, 2, function(x){
    temp <- is.nan(x)
    sum(temp)
  })

  cat(paste0('The number of inf: ', sum(n_inf), '\n'),
      paste0('The number of NA: ', sum(n_na), '\n'),
      paste0('The number of NAN: ', sum(n_nan), '\n')
  )

  # if (sum(n_inf) > 0 | sum(n_na) > 0| sum(n_nan) > 0) {
  #   stop("Inf/NA/NAN existed!\n")
  # }

  # if(exists(".Random.seed")) {rm(.Random.seed)}

  result <- impute::impute.knn(as.matrix(raw_matrix),
                               k = k,
                               rowmax = rowmax,
                               colmax = colmax,
                               maxp = maxp,
                               rng.seed = rng.seed)

  result <- result[[1]]
  result <- data.frame(name=label,
                       result,
                       stringsAsFactors = F)

  return(result)
}


#' @title CalcuRsquare
#' @description calcualte r2 for regression
#' @author Zhiwei Zhou
#' @param y numeric
#' @param y.pred numeric
#' @export

setGeneric('CalcuRsquare',
           def = function(y, y.pred){
             result <- summary(lm(y, y.pred))$r.squared
             return(result)
           }
)
