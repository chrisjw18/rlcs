#'@title plot_rlcs
#'@description Plot multiple rlcs on one plot
#'@param dataframe dataframe object containing the data
#'@param par character string, dataframe column name containing par levels
#'@param parameter character string, dataframe column name containing desired parameter data, e.g. etr, npq etc.
#'@param treatment character string, dataframe column name containing treatment levels
#'@param error character string, dataframe column name containing SE data
#'@return plot showing all treatments parameter (+/- SE)
#'@import magrittr
#'@import randomcoloR
#'@export

plot_rlcs <- function(dataframe=NA, par='par', parameter='etr', treatment=NA, error=NA){

  lev <- dataframe[,treatment] %>% levels
  n.lev <- lev %>% length
  length.each <- nrow(dataframe) / n.lev
  max.val <- max(dataframe[,parameter])
  my.cols <- randomcoloR::distinctColorPalette(length(lev))

  x <- subset(dataframe, dataframe[,treatment] == lev[1])

  par(mar=c(3,3,1,1), mgp=c(1.8,0.3,0), tck=-0.01, las=1)
  plot(x[,par], x[,parameter], ylim=c(0,max.val*1.2), xlab=expression(paste('PAR (',mu,'mol photons ',m^2,' ',s^-1,')')), ylab=toupper(parameter), col=my.cols[1], pch=19, type='o')
  if(!is.na(error)){
    arrows(x[,par], x[,parameter], x[,par], x[,parameter]+x[,error], length=0.05, angle=90, col=my.cols[1])
    arrows(x[,par], x[,parameter], x[,par], x[,parameter]-x[,error], length=0.05, angle=90, col=my.cols[1])
  }#end of if clause

  for(i in 2:length(lev)){
    x <- subset(dataframe, dataframe[,treatment] == lev[i])
    points(x[,par], x[,parameter], type='o', pch=19, col=my.cols[i])
    if(!is.na(error)){
    arrows(x[,par], x[,parameter], x[,par], x[,parameter]+x[,error], length=0.05, angle=90, col=my.cols[i])
    arrows(x[,par], x[,parameter], x[,par], x[,parameter]-x[,error], length=0.05, angle=90, col=my.cols[i])
    }#end of if clause
  }#end of for loop

  legend('bottomright', lty=1, pch=19, col=my.cols, legend=paste('Treatment', lev), bty='n', horiz=F)

}#end of function
