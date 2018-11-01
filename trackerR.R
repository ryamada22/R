###
#
# Rscript --slave --vanilla trackerR.R args0 args1 args2
#
###
library(EBImage)
library(png)
library(jpeg)
args <- commandArgs(TRUE)
if(length(args) < 2){
  cat(
"command line args are
args[0]: image sequence folder. format is currently .png
args[1]: output filename. space separated txt. ID,t,x,y
args[2]: optional. filepath of temporaly result of tracking. you can restart from this file.\n"
      )
  do <- FALSE
} else if (length(args) > 3){
  cat(
"args is less than 4.\n"
)
  do <- FALSE
} else {
  do <- TRUE
}

input_dir <- args[1]
out_name  <- args[2]
if(length(args) == 3){
  tmp_res <- read.table(args[3], header=TRUE)
  final <- split(tmp_res[,-1], tmp_res$ID)
} else {
  final <- list()
}

txt0 <- "Do you continue next tracking ?\nYes: Left click\nNo: Right click"
txt1 <- "Do you accept your tracking ?\nYes: Left click\nNo: Right click"

imgs <- list.files(input_dir, full.names=TRUE)
ans <- TRUE
accept <- FALSE
if(do){
  while(ans){
    while(!accept){
      res <- matrix(0, length(imgs), 3)
      for(i in seq(imgs)){
        f <- readImage(imgs[i])
        xleft <- 1
        ybottom <- nrow(f)
        xright <- ncol(f)
        ytop <- 1
        par(mar=rep(0, 4))
        image(t(0), xlim=c(xleft, xright), ylim=c(ybottom, ytop), col="white")
        rasterImage(f, xleft, ybottom, xright, ytop)
        legend("topright", legend=paste0("t = ", i), text.col="yellow", cex=1.5)
        if(length(final) > 0){
          for(j in seq(final)){
            points(final[[j]][i, 2], final[[j]][i, 3], col="yellow", pch=4, cex=1, font=2)
          }
        }
        xy <- unlist(locator(1))
        res[i, ] <- switch(is.numeric(xy)+1, c(i, NA, NA), c(i, xy))
      }
      text((xleft+xright)/2, (ybottom+ytop)/2, txt1, col="yellow", pos=3, cex=2, font=2)
      accept <- is.list(locator(1))
    }
    final <- c(final, list(res))
    image(t(0), xlim=c(xleft, xright), ylim=c(ybottom, ytop), col="white")
    rasterImage(f, xleft, ybottom, xright, ytop)
    text((xleft+xright)/2, (ybottom+ytop)/2, txt0, col="yellow", pos=3, cex=2, font=2)
    ans <- is.list(locator(1))
    accept <- FALSE
  }
  output <- cbind(rep(seq(final), each=length(imgs)), do.call(rbind, final))
  colnames(output) <- c("ID", "t", "x", "y")
  write.table(output, out_name, quote=FALSE, row.names=FALSE)
}

