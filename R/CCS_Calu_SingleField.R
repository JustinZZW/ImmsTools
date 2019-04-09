#' @title CCS_Calu_SingleField
#' @description Calculate CCS values using single field method
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#'
#' @usage
#' \code{CCS_Calu_SingleField(polarity="pos", charge=1, Agilent.table="V1")}
#'
#' @param raw.data Character. The excel file name. The top 3 columns should be named as: 'M', 'Adducts', 'Arrival.time'
#' @param polarity Character. Ionzation polarity of experiments. Should be one of "pos" or "neg". "pos" denotes "positive mode", "neg" denotes "negative mode". Default is "pos".
#' @param charge Numeric. The charge of calibrants. Default is 1.
#' @param Agilent.table Character. The verison of Agilent tuning mix reference table is used to establish calibration Scurve. "V1" denotes the reference table in Agilent IM-MS browser B.07.01; "V2" denotes the reference table reported by Stow et. al. Anal. Chem. 2017
#' @return a csv table includs, m/z and CCS values etc.
#' @export
#' @examples
#' CCS_Calu_SingleField(polarity="pos", charge=1, Agilent.table="V1")
#'



CCS_Calu_SingleField <- function(raw.data='raw data file.xlsx',
                                 polarity="pos",
                                 charge=1,
                                 Agilent.table="V1",
                                 TFix=NULL,
                                 Beta=NULL,
                                 output.filename="CCS results.csv"){
  # import data=====================================================
  cat("Import feature table"); cat("\n")
  raw.data <- readxl::read_excel(raw.data)
  # load(system.file("extdata", "adduct_table.RData", package = "MetIMMS"))
  # load(system.file("list", "adduct_table.RData", package = "MetIMMS"))

  if (any(is.null(TFix), is.null(Beta))) {
    cat("Build calibration curve using single field method\n")
    coefficient <- Cal_SingleField(polarity = polarity,
                                   charge = charge,
                                   Agilent.table = Agilent.table)

    cat("Start calculate CCS values\n")

    trans.data <- mz_transform(raw.data = raw.data, polarity = polarity)
    gama <- sqrt(trans.data$mz/(trans.data$mz+28.0061))
    CCS <- (trans.data$Arrival.time-coefficient["TFix"])/(coefficient["Beta"]/charge*gama)

    cat("Export results\n")
    output.result <- data.frame(trans.data, CCS=CCS)
  } else {

    cat("Start calculate CCS values\n")
    coefficient <- c(as.numeric(TFix), as.numeric(Beta))
    names(coefficient) <- c('TFix', 'Beta')

    trans.data <- mz_transform(M = raw.data$mz, polarity = polarity)
    gama <- sqrt(trans.data$mz/(trans.data$mz+28.0061))
    CCS <- (trans.data$Arrival.time-coefficient["TFix"])/(coefficient["Beta"]/charge*gama)
  }

  temp.file.name <- paste(getwd(), output.filename, sep = "/")
  write.csv(output.result, temp.file.name, row.names = F)
  return(output.result)
}




#' @title Cal_SingleField
#' @description Establish calibration curve using single field method
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#'
#' @usage
#' \code{Cal_SingleField(polarity="pos", charge=1, Agilent.table="V1")}
#' @param raw.data the data.frame including following columns
#' @param polarity Character. Ionzation polarity of experiments. Should be one of "pos" or "neg". "pos" denotes "positive mode", "neg" denotes "negative mode". Default is "pos".
#' @param charge Numeric. The charge of calibrants. Default is 1.
#' @param Agilent.table Character. The verison of Agilent tuning mix reference table is used to establish calibration Scurve. "V1" denotes the reference table in Agilent IM-MS browser B.07.01; "V2" denotes the reference table reported by Stow et. al. Anal. Chem. 2017
#' @return "TFix" is the coefficient denoting the mobility-indenpend flight time, "Beta" is the coefficient denoting the instrument-dependent proportionality cofficient.
#' @export
#' @examples
#' Cal_SingleField(polarity="pos", charge=1, Agilent.table="V1")
#'


Cal_SingleField <- function(raw.data=NULL,
                            polarity="pos",
                            charge=1,
                            Agilent.table="V1",
                            file.name = "Single_Field_cal_result") {
  # Load data ==================================================================
  if (is.null(raw.data)) {
    if (!file.exists("calibration file table.xlsx")) {
      stop("There is no calibration file for CCS calibration.")
    }

    cal.table <- readxl::read_excel("calibration file table.xlsx")
  }

  ref.table <- get_ref_table(Agilent.table = Agilent.table, polarity = polarity)

  # Build calibration curve ====================================================
  idx.at <- which(!is.na(cal.table$`Arrival time`)) # index of arrive time
  cal.at <- cal.table$`Arrival time`[idx.at]

  ref.table <- ref.table[idx.at,]
  gama <- sqrt(ref.table$mz/(ref.table$mz+28.0061))
  ref.CCS <- ref.table$CCS
  X <- ref.CCS*gama

  reg.cal <- lm(cal.at~X)
  reg.cal.coefficient <- reg.cal$coefficients
  names(reg.cal.coefficient) <- c("TFix", "Beta")
  R.squred <- summary(reg.cal)$r.squared
  Residuals <- reg.cal$residuals
  CalCondition <- c(Agilent.table, polarity)
  names(CalCondition) <- c("Reference table verison", "Polarity")

  # Output result ==============================================================
  calibration.result <- list(coefficient=reg.cal.coefficient, Condition=CalCondition, r.squred=R.squred, residuals=Residuals, model=reg.cal)
  file.name <- paste(getwd(), "Single_Field_cal_result", sep = "/")
  save(calibration.result, file = file.name)

  return(calibration.result$coefficient)

}



get_ref_table <- function(Agilent.table="V1", polarity="pos"){
  if (Agilent.table=="V1") {
    if (polarity=="pos") {
      ref.table <- Agilent_RefTable_V1$pos
    } else {
      ref.table <- Agilent_RefTable_V1$neg
    }
  }

  if (Agilent.table=="V2") {
    if (polarity=="pos") {
      ref.table <- Agilent_RefTable_V2$pos
    } else {
      ref.table <- Agilent_RefTable_V2$neg
    }
  }

  return(ref.table)
}

