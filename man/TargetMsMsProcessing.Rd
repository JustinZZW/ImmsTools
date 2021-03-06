% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtractTargetMsMsSpectra.R
\name{TargetMsMsProcessing}
\alias{TargetMsMsProcessing}
\title{TargetMsMsProcessing}
\usage{
TargetMsMsProcessing(msms_file = "L0656_H_pos_10ppm.mzXML",
  root_path = ".", polarity = NULL, adduct = NULL,
  concentration = NULL, mode = c("SPMCE", "MPSCE", "SPSCE"),
  ce_levels = c("CE10", "CE20", "CE40"))
}
\arguments{
\item{msms_file}{mzXML file format}

\item{root_path}{Default: '.'}

\item{polarity}{'positive', 'negative'. If it is not given (NULL), it would be extracted from msms_file name. Default: NULL}

\item{adduct}{If it is not given (NULL), it would be extracted from msms_file name. Default: NULL}

\item{concentration}{If it is not given (NULL), it would be extracted from msms_file name. Default: NULL}

\item{mode}{SPMCE: single precursor with multiple CEs; MPSCE: multiple precursors, but single CE for each precursor; SPSCE: single precursor with one CE;}

\item{ce_levels}{collison energy. Default: c("CE10", "CE20", "CE40")}
}
\description{
TargetMsMsProcessing
}
\examples{
setwd('F:/01 MetIMMS/00 data processing/190702 extract msms spectra from targeted 6560 data/convert_msms')
# Standard
TargetMsMsProcessing(msms_file = "L0656_H_pos_10ppm.mzXML", root_path = '..', mode = 'SPMCE', ce_levels = c("CE10", "CE20", "CE40"))
# QC
TargetMsMsProcessing(msms_file = "QC00.mzXML", root_path = '..', polarity = 'Positive', adduct = '[M+H]+', concentration = '10ppm',mode = 'MPMCE', ce_levels = c("CE20", "CE40"))
}
\author{
Zhiwei Zhou
}
