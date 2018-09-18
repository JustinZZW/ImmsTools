#' @title msms_plot
#' @author Zhiwei Zhou
#' \email{zhouzw@@sioc.ac.cn}
#' @description plot MS/MS mirror plot.
#' @param raw.exp.spec Matrix or data frame. The experiment MS/MS spectrum. The 1st column is "mz", the 2nd column is "intensity".
#' @param raw.std.spec Matrix or data frame. The library MS/MS spectrum. The 1st column is "mz", the 2nd column is "intensity"
#' @param match.result output from SpecMatch
#' @param dir.path the file path of output
#' @param cmp.name Character. The name of compound
#' @param is.output Logical. Whether output the plot?
#' @param pg.format png or pdf format.
#' @examples
#' raw.exp.spec <- data.frame(mz=c(42.0364, 58.0664, 59.0745, 118.0858), intensity=c(0.01, 1, 0.4, 0.25))
#' raw.std.spec <- data.frame(mz=c(42.0364, 58.0664, 75.4568, 118.0858), intensity=c(0.01, 1, 0.2, 0.25))
#' msms_plot(raw.exp.spec, raw.std.spec, is.output=FALSE)
#' @export

msms_plot <- function(raw.exp.spec,
                      raw.std.spec,
                      # lib.meta,
                      match.result=NULL,
                      dir.path=NULL,
                      cmp.name=NULL,
                      # exp.spec.annotate=NULL,
                      # std.spec.annotate=NULL,
                      is.output=c(TRUE, FALSE),
                      pg.formate=c("png", "pdf")) {

  # is.output <- match.arg(is.output)
  # pg.formate <- match.arg(pg.formate)

  temp.exp.spec <- as.data.frame(raw.exp.spec)
  temp.std.spec <- as.data.frame(raw.std.spec)

  exp.max.int <- max(as.numeric(temp.exp.spec$intensity))
  std.max.int <- max(as.numeric(temp.std.spec$intensity))
  temp.exp.spec$intensity <- temp.exp.spec$intensity/exp.max.int
  temp.std.spec$intensity <- temp.std.spec$intensity/std.max.int

  mz <- c(temp.exp.spec$mz, temp.std.spec$mz)
  int <- c(temp.exp.spec$intensity, -temp.std.spec$intensity)

  # extract MS/MS matched information
  if (is.null(match.result)) {
    color <- c(rep("dodgerblue", length(c(temp.exp.spec$mz))), rep("firebrick1", length(c(temp.std.spec$mz))))
  } else {
    score <- round(as.numeric(match.result), digits = 4)
    information <- attributes(match.result)

    is.matched.exp <- as.character(temp.exp.spec$mz) %in%
      as.character(information$spec$exp[which(information$spec$exp[,2]>0 & information$spec$ref[,2]>0)
                                        ,1])
    # as.character(round(information$spec$exp[,1], 5))

    is.matched.std <- as.character(temp.std.spec$mz) %in%
      as.character(information$spec$ref[which(information$spec$exp[,2]>0 & information$spec$ref[,2]>0)
                                        ,1])

    color <- rep('gray', length(mz))
    color[which(is.matched.exp)] <- "dodgerblue"
    color[which(is.matched.std) + length(is.matched.exp)] <- "firebrick1"
  }

  # color[infor.annotate=="precursor"] <- "gray"

  x.range <- c(0.95*min(as.numeric(mz)), 1.05*max(as.numeric(mz)))


  if (is.output) {
    file.name <- gsub("\\:", "_", gsub("\\/", "_", cmp.name))
    # file.name <- paste(file.name, cmp.name, sep = "_")

    switch(pg.formate,
           "png"={
             file.name <- paste(file.name,".png", sep = "")
             file.name <- file.path(dir.path, file.name)
             png(filename = file.name, width = 800, height = 400)
           },
           "pdf"={
             file.name <- paste(file.name,".pdf", sep = "")
             file.name <- file.path(dir.path, file.name)
             pdf(file=file.name, width = 12, height = 6)
           })
  }


  par(mar=c(3, 3, 2, 2), oma=c(0.5, 0.5, 0.5, 0.25))

  plot(mz, int,
       type = "h",
       lwd = 1.5,
       xlim = x.range,
       ylim = c(-1.1, 1.1),
       col = color,
       xlab = "m/z",
       ylab = "Relative intensity",
       # main = information[1],
       mgp = c(1.8, 0.5, 0),
       cex.main = 1,
       cex.lab = 1,
       cex.axis = 0.8)

  abline(h=0, lwd=1)

  legend("topleft",
         legend = c("Matched peaks: experiment", "Matched peaks: library"),
         pch = 19, col = c("dodgerblue", "firebrick1"), bty="n", cex = 0.6)

  if (!is.null(match.result)) {
    frag.annotate <- c(is.matched.exp, is.matched.std)
    points(x = mz[frag.annotate],
           y = int[frag.annotate],
           pch = 19,
           cex = 0.5,
           # col = "black",
           col = c(rep("dodgerblue", sum(is.matched.exp)),
                   rep("firebrick1", sum(is.matched.std))))

    legend('bottomleft',
           legend = c(paste('Score:', round(as.numeric(score), 4), sep = ''),
                      paste('Name:', cmp.name, sep = '')),
           bty = 'n',
           cex = 0.6)

  }

  if (is.output) {
    dev.off()
  }

}
