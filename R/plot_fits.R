#'@title plot_fits
#'@description Plot multiple e and p fits with u/l bounds on one plot
#'@param dataframe dataframe object containing the data
#'@param par character string, dataframe column name containing par levels
#'@return plot showing all fits contained in the dataframe (+/- upper and lower CI)
#'@export


plot_fits <- function(dataframe=NA, par='par'){

  fits <- colnames(dataframe)[which(grepl('fit', colnames(dataframe)))]
  upper <- colnames(dataframe)[which(grepl('upper', colnames(dataframe)))]
  lower <- colnames(dataframe)[which(grepl('lower', colnames(dataframe)))]
  my.cols <- randomcoloR::distinctColorPalette(length(fits))

  plot.par <- dataframe[,par]
  y <- dataframe[, fits[1]]
  up <- dataframe[,upper[1]]
  low <- dataframe[,lower[1]]

  max.val <- dataframe[,-which(colnames(dataframe)==par)] %>% max

  par(mar=c(3,3,1,1), mgp=c(1.8,0.3,0), tck=-0.01, las=1, mfrow=c(1,1))
  plot(plot.par, y, ylim=c(0,max.val*1.2), xlab=expression(paste('PAR (',mu,'mol photons ',m^2,' ',s^-1,')')), ylab='rETR', col=my.cols[1], pch=19, type='l', lwd=1)
  polygon(c(plot.par, rev(plot.par)), c(up, rev(low)), col=addTrans(my.cols[1], 50), border=NA)


  for(i in 2:length(fits)){
    y <- dataframe[, fits[i]]
    up <- dataframe[,upper[i]]
    low <- dataframe[,lower[i]]
    points(plot.par, y, type='l', lwd=1, col=my.cols[i])
    polygon(c(plot.par, rev(plot.par)), c(up, rev(low)), col=addTrans(my.cols[i], 50), border=NA)
  }#end of for loop

  key <- gsub('_fit', replacement = '', x = fits)
  legend('bottomright', lty=1, col=my.cols, legend=key, bty='n', horiz=F)

}#end of function
