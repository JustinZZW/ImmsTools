#' @title extractClyssyFire
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @description extract sdf file exported from ClassyFire
#' @param filename the sdf file name
#' @importFrom magrittr '%>%'
#' @export


extractClassyFire <- function(filename){
  temp <- readr::read_lines(filename)

  idx.kegg.id <- which(temp=='$$$$') + 1
  idx.kegg.id <- c(1, idx.kegg.id)
  idx.inchikey <- which(temp=='> <InChIKey>') + 1

  if (length(idx.kegg.id) == length(idx.inchikey)) {
    idx.inchikey <- which(temp=='> <InChIKey>')+1
    idx.smiles <- which(temp=='> <SMILES>')+1
    idx.kingdom <- which(temp== '> <Kingdom>')+1
    idx.superclass <- which(temp=="> <Superclass>") + 1
    idx.class <- which(temp=="> <Class>") + 1
    idx.subclass <- which(temp=="> <Subclass>") + 1
  } else {
    less.num <- length(idx.kegg.id) - length(idx.inchikey)
    NA.num <- rep(NA, less.num)

    idx.inchikey <- which(temp=='> <InChIKey>') + 1
    idx.inchikey <- c(idx.inchikey, NA.num)

    idx.smiles <- which(temp=='> <SMILES>') + 1
    idx.smiles <- c(idx.smiles, NA.num)

    idx.kingdom <- which(temp== '> <Kingdom>') + 1
    idx.kingdom <- c(idx.kingdom, NA.num)

    idx.superclass <- which(temp=="> <Superclass>") + 1
    idx.superclass <- c(idx.superclass, NA.num)

    idx.class <- which(temp=="> <Class>") + 1
    idx.class <- c(idx.class, NA.num)

    idx.subclass <- which(temp=="> <Subclass>") + 1
    idx.subclass <- c(idx.subclass, NA.num)
  }



  temp_result <- data.frame(kegg.id=temp[idx.kegg.id],
                            inchikey=temp[idx.inchikey],
                            smiles=temp[idx.smiles],
                            kingdom=temp[idx.kingdom],
                            superclass=temp[idx.superclass],
                            class=temp[idx.class],
                            subclass=temp[idx.subclass],
                            stringsAsFactors = F)

  idx.na <- apply(temp_result, 2, function(x){
    temp <- which(nchar(x)==0)
  })

  for (i in 1:6) {
    if (length(idx.na[[i]]) > 0 ) {
      temp_result[idx.na[[i]],i] <- NA
    }
  }

  temp_result
}
