#' @title MdCleaning
#' @author Zhiwei Zhou
#' @description Molecular Descriptor Cleaning to remove NA and same molecular descriptors
#' @param raw.data a data.frame of MD, each column represent one MD, and each row represent one compound
#' @param col.number.res the column index to reserve
#' @export
#' @examples MdCleaning

MdCleaning <- function(raw.data,
                       col.number.res){
  temp.data <- raw.data[, -col.number.res]

  # delete descriptors which na appear in 50% samples
  idx <- apply(temp.data, 2, function(x){
    number.na <- sum(is.na(x))
    if (number.na >= round(0.5*length(x), 0)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  idx <- which(idx)
  temp.data <- temp.data[,-idx]
  idx <- which(apply(temp.data, 2, sd) == 0) # remove same variables
  temp.data <- temp.data[,-idx]
  result <- data.frame(raw.data[,col.number.res], temp.data)
}
