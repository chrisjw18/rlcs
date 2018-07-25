#'@title plot_params
#'@description Plot RLC average parameters (fvfm, rETRmax, alpha, ek) +/- SE after calculation with rlc_fit
#'@param dataframe av_params dataframe derived from fit_rlcs
#'@param treatment character string, dataframe column name containing treatment levels
#'@return 4 panel plot showing average Fv/Fm, rETRmax, alpha and Ek (+/- SE) per treatment
#'@export

plot_params <- function(dataframe = NA, treatment = NA){

  my.cols <- randomcoloR::distinctColorPalette(nrow(dataframe))

  par(mar=c(3,3,1,1), mgp=c(1.8,0.3,0), tck=-0.01, las=1, mfrow=c(2,2))

  p1 <- barplot(as.matrix(dataframe[,'fvfm']), beside=T, space=0.1, names.arg=dataframe[,treatment], ylim=c(0, (max(dataframe[,'fvfm'])+max(dataframe[,'fvfm_se']))*1.2), col = my.cols, ylab='Fv/Fm')
  arrows(p1, dataframe[,'fvfm'], p1, dataframe[,'fvfm']+dataframe[,'fvfm_se'], length=0.05, angle=90)

  p1 <- barplot(as.matrix(dataframe[,'etrmax']), beside=T, space=0.1, names.arg=dataframe[,treatment], ylim=c(0, (max(dataframe[,'etrmax'])+max(dataframe[,'etrmax_se']))*1.2), col = my.cols, ylab='rETRmax')
  arrows(p1, dataframe[,'etrmax'], p1, dataframe[,'etrmax']+dataframe[,'etrmax_se'], length=0.05, angle=90)

  p1 <- barplot(as.matrix(dataframe[,'alpha']), beside=T, space=0.1, names.arg=dataframe[,treatment], ylim=c(0, (max(dataframe[,'alpha'])+max(dataframe[,'alpha_se']))*1.2), col = my.cols, ylab=expression(alpha))
  arrows(p1, dataframe[,'alpha'], p1, dataframe[,'alpha']+dataframe[,'alpha_se'], length=0.05, angle=90)

  p1 <- barplot(as.matrix(dataframe[,'ek']), beside=T, space=0.1, names.arg=dataframe[,treatment], ylim=c(0, (max(dataframe[,'ek'])+max(dataframe[,'ek_se']))*1.2), col = my.cols, ylab='Ek')
  arrows(p1, dataframe[,'ek'], p1, dataframe[,'ek']+dataframe[,'ek_se'], length=0.05, angle=90)

}
