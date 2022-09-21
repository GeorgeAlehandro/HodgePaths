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
#' @param fc function for averiging vectors in cell ("mean"/"median")
#'
#' @return
#' @export
#'
#' @examples
vplot <- function(df,seg,normalize=FALSE,fc="mean") {
  minimum<-min(df$x)
  maximum<-max(df$x)
  spacing <- (maximum-minimum)/(seg-1)

  df$cutX <- cut(df$x, breaks = seq(minimum-1,maximum+1,spacing))
  df$cutY <- cut(df$y, breaks = seq(min(df$y)-1,max(df$y)+1,spacing))

  if (fc=="mean") {
  mergedy <- ddply(df, .(cutX,cutY), summarize, mean=mean(dy))
  mVal = ddply(df, .(cutX, cutY), summarize, mean = mean(dx)) }
  if (fc=="median") {
    mergedy <- ddply(df, .(cutX,cutY), summarize, mean=median(dy))
    mVal = ddply(df, .(cutX, cutY), summarize, mean = median(dx)) }



  mVal$cutX = as.character(mVal$cutX)
  mVal$cutX = unlist(strsplit(mVal$cutX, ","))[c(T,F)]
  mVal$cutX = as.numeric(sub("\\(", "", mVal$cutX))+spacing/2
  mVal$cutY = as.character(mVal$cutY)
  mVal$cutY = unlist(strsplit(mVal$cutY, ","))[c(T,F)]
  mVal$cutY = as.numeric(sub("\\(", "", mVal$cutY))+spacing/2
  #coordinates(mVal) = ~cutX +cutY
  #gridded(mVal) = TRUE
  #mergedx <- plyr::ddply(df, .(cutX,cutY), summarize, mean=mean(dx))


  final=data.frame(x=mVal$cutX,y=mVal$cutY, dx=mVal$mean,dy=mergedy$mean)
  final$vel=final$dx^2+final$dy^2
  if (normalize==TRUE) {
    final$dx <- (final$dx/sqrt(final$vel))*spacing
    final$dy <- (final$dy/sqrt(final$vel))*spacing
  }

  ggplot(final, aes(x = x, y = y))+scale_fill_gradientn(colours=c("darkblue","darkgreen","yellow"))+geom_raster(aes(fill=vel))+geom_segment(aes(xend = x + dx, yend = y + dy),size = 0.3, arrow = arrow(length = unit(0.2, "cm")))

}







#' Title
#'
#' @param df dataframe with $x $y $dx, $dy being a vector field
#' @param seg #into how many segments do you want to average out a vector field
#' @param normalize TRUE or FALSE if you wish to normalize vectors when they are displayed
#' @param fc function for averiging vectors in cell ("mean"/"median")
#' @param jitter amount of jitter
#' @param n amount of lines
#' @param mask Draws vector field over it (T or F)
#'
#' @return
#' @export
#'
#' @examples
vplot2 <- function(df,seg,normalize=FALSE,fc="mean",jitter=130,n=70,mask=F) {
  minimum<-min(df$x)
  maximum<-max(df$x)
  spacing <- (maximum-minimum)/(seg-1)

  df$cutX <- cut(df$x, breaks = seq(minimum-1,maximum+1,spacing))
  df$cutY <- cut(df$y, breaks = seq(min(df$y)-1,max(df$y)+1,spacing))

  if (fc=="mean") {
    mergedy <- ddply(df, .(cutX,cutY), summarize, mean=mean(dy))
    mVal = ddply(df, .(cutX, cutY), summarize, mean = mean(dx)) }
  if (fc=="median") {
    mergedy <- ddply(df, .(cutX,cutY), summarize, mean=median(dy))
    mVal = ddply(df, .(cutX, cutY), summarize, mean = median(dx)) }



  mVal$cutX = as.character(mVal$cutX)
  mVal$cutX = unlist(strsplit(mVal$cutX, ","))[c(T,F)]
  mVal$cutX = as.numeric(sub("\\(", "", mVal$cutX))+spacing/2
  mVal$cutY = as.character(mVal$cutY)
  mVal$cutY = unlist(strsplit(mVal$cutY, ","))[c(T,F)]
  mVal$cutY = as.numeric(sub("\\(", "", mVal$cutY))+spacing/2


  final=data.frame(x=mVal$cutX,y=mVal$cutY, dx=mVal$mean,dy=mergedy$mean)
  final$vel=final$dx^2+final$dy^2
  if (normalize==TRUE) {
    final$dx <- (final$dx/sqrt(final$vel))*spacing
    final$dy <- (final$dy/sqrt(final$vel))*spacing
  }

  expanze <- expand.grid(x=unique(final$x),y=unique(final$y))

  vysledek<-merge(expanze, final, by.x=c("x", "y"),by.y=c("x", "y"),all.x=TRUE)
  vysledek[is.na(vysledek)] <- 0

  plot1<-ggplot(vysledek, aes(x=x,y=y))+scale_fill_gradientn(colours=c("darkblue","darkgreen","yellow"))+geom_raster(aes(fill=vel))+ metR::geom_streamline(data = vysledek, aes(dx = dx, dy = dy,color = vel),L = 5, res = 1.5, n =n, jitter =jitter,arrow = arrow(length = unit(0.15, "cm")))

  if (mask==T) {
  plot2<-ggplot(vysledek, aes(x=x, y=y)) +metR::geom_streamline(aes(dx = dx, dy = dy, alpha = ..step..,color = sqrt(..dx..^2 + ..dy..^2)),n=n,L = 5,jitter=jitter, res = 1.5, arrow=NULL,lineend = "round") +scale_size(range = c(0.2, 1.2))+theme_bw()+geom_segment(data=final, aes(xend = x + dx, yend = y + dy),size = 0.1, arrow = arrow(length = unit(0.2, "cm"),ends="first", type = "open"))
  }
  if (mask==F) {
    plot2<-ggplot(vysledek, aes(x=x, y=y)) +metR::geom_streamline(aes(dx = dx, dy = dy, alpha = ..step..,color = sqrt(..dx..^2 + ..dy..^2)),n=n,L = 5,jitter=jitter, res = 1.5, arrow=NULL,lineend = "round") +scale_size(range = c(0.2, 1.2))+theme_bw()
  }
  grid.arrange(plot1, plot2, ncol=2)

}

