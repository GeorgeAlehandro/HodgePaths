#' Make dataframe from Whintney interpolation
#'
#' @param out Output from fucntion WhitneyMap1f_BC
#'
#' @return returns dataframe that one can put into plotting function vplot
#' @export
#'
#' @examples
makedf <- function(out) {
  df <- data.frame(x=out$BC[,1],y=out$BC[,2],dx=out$oneform[,1],dy=out$oneform[,2])
  return(df)
}




#' Plot vector field
#'
#' @param df dataframe with $x $y $dx, $dy being a vector field
#' @param seg #into how many segments do you want to average out a vector field
#' @param normalize TRUE or FALSE if you wish to normalize vectors when they are displayed
#'
#' @return
#' @export
#'
#' @examples
vplot <- function(df,seg,normalize=FALSE) {
  library(ggplot2)
  library(plyr)
  minimum<-min(df$x)
  maximum<-max(df$x)
  spacing <- (maximum-minimum)/(seg-1)

  df$cutX <- cut(df$x, breaks = seq(minimum-1,maximum+1,spacing))
  df$cutY <- cut(df$y, breaks = seq(min(df$y)-1,max(df$y)+1,spacing))



  mVal = ddply(df, .(cutX, cutY), summarize, mean = mean(dx))
  mVal$cutX = as.character(mVal$cutX)
  mVal$cutX = unlist(strsplit(mVal$cutX, ","))[c(T,F)]
  mVal$cutX = as.numeric(sub("\\(", "", mVal$cutX))+spacing/2
  mVal$cutY = as.character(mVal$cutY)
  mVal$cutY = unlist(strsplit(mVal$cutY, ","))[c(T,F)]
  mVal$cutY = as.numeric(sub("\\(", "", mVal$cutY))+spacing/2
  #coordinates(mVal) = ~cutX +cutY
  #gridded(mVal) = TRUE
  #mergedx <- plyr::ddply(df, .(cutX,cutY), summarize, mean=mean(dx))
  mergedy <- ddply(df, .(cutX,cutY), summarize, mean=mean(dy))

  final=data.frame(x=mVal$cutX,y=mVal$cutY, dx=mVal$mean,dy=mergedy$mean)
  final$vel=final$dx^2+final$dy^2
  if (normalize==TRUE) {
    final$dx <- (final$dx/sqrt(final$vel))*spacing
    final$dy <- (final$dy/sqrt(final$vel))*spacing
  }

  ggplot(final, aes(x = x, y = y))+scale_colour_continuous(low = "grey80", high = "darkred")+geom_raster(aes(fill=vel))+geom_segment(aes(xend = x + dx, yend = y + dy),size = 0.3, arrow = arrow(length = unit(0.2, "cm")))

  #ggplot(final, aes(x=x,y=y)) +geom_raster(aes(fill=vel))+ metR::geom_streamline(data = final, aes(dx = dx, dy = dy),L = 1.75, res = 1, n = 60, jitter =150)

}
