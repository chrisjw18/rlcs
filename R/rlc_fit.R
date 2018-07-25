#'@title rlc_fit
#'@description simple fit of e and p model to PAM data.
#'@param dataframe name of the dataframe containing the data
#'@param par character string of dataframe column containing PAM par levels
#'@param yield character string, dataframe column containing Y(PSII) data
#'@param f character string, dataframe column containing F data
#'@param fm. character string, dataframe column containing Fm' data
#'@param sample character string, dataframe column containing sample names per light curve
#'@param treatment character string, dataframe column containing treatment levels if applicable. If provided, averages are made per level of treatment for variables (rETR, NPQ, YNPQ, YNO) and parameters (ETRmax, alpha, Ek, Fv/Fm)
#'@param predict_ul verbose (T/F), provide predictions of fit, upper and lower bounds for e and p models
#'@return rlc_results list containing the following dataframes:
#'    1. 'rlcs' containing the original data, with calculate variables (rETR, NPQ, YNPQ, YNO)
#'       and e and p model parameters (a,b,c) and derived parameters (rETRmax, alpha, Ek)
#'
#'    If predict_ul ==T:
#'    2. 'fits' If predict_ul == Tcontaining predicted fit, upper and lower bounds for each
#'       light curve e and p model calculated using predictNLS function
#'
#'    If 'treatment' is provided:
#'    3. 'av_variables' containing averages of rETR, NPQ, YNPQ and YNO per treatment,
#'        including standard deviation and standard error
#'    4. 'av_params' containing averages of rETRmax, alpha, ek and fvfm per treatment,
#'        including standard deviation and standard error
#'    5. 'av_fits' containing e and p fit of etr ~ par with upper and lower bounds per
#'        treatment.
#'@import magrittr
#'@export


rlc_fit<-function(dataframe=NA,par = 'par', yield = 'y', f = 'f', fm. = 'fm.', sample = NA, treatment = NA, predict_ul = T){

  source('./R/predictNLS.R')


  #calculate number of light curves in dataframe
  no.curves <- nrow(dataframe) / (dataframe[,par] %>% unique %>% length)
  dataframe$lc <- rep(1:no.curves, each = (dataframe[,par] %>% unique %>% length)) %>% factor
  lev <- dataframe$lc %>% levels

  #calculate rETR
  dataframe$etr <- dataframe[, yield] * dataframe[,par] * 0.5

  #blank columns for below
  dataframe$npq <- NA
  dataframe$ynpq <- NA
  dataframe$yno <- NA
  dataframe$a <- NA
  dataframe$b <- NA
  dataframe$c <- NA
  dataframe$etrmax <- NA
  dataframe$alpha <- NA
  dataframe$ek <- NA

  if(predict_ul ==T){
    max.par <- max(dataframe[,par])
    pred.catch <- data.frame(par = seq(0,max.par,10))
  }

  #calculate NPQs and e and ps
  for(i in 1:length(lev)){

    dat <- subset(dataframe, lc == lev[i])

    if(is.na(sample)){
      current.sample <- paste('lc', lev[i], sep='_')
    } else {
      current.sample <- dat[,sample][1]
    }

    #NPQ = (Fm - Fm')/Fm'
    npq<-(dat[,fm.][1]-dat[,fm.])/dat[,fm.]
    dataframe$npq[dataframe$lc == lev[i]]<-npq
    rm(npq)

    #YNO = F (actual)/fm (first)
    yno<-dat[,f]/dat[,fm.][1]
    dataframe$yno[dataframe$lc == lev[i]]<-yno

    #YNPQ = 1 - Yield - YNO
    ynpq<-1-dat[,yield]-yno
    ynpq[ynpq<0] <- 0
    dataframe$ynpq[dataframe$lc == lev[i]]<-ynpq
    rm(yno)
    rm(ynpq)

    #Eileers and Peters model
    #E and P model
    a0<-0.00003; b0<--0.003; c0<-4
    x <- dat[,par]
    y <- dat[,'etr']
    m1 <- nls(y~x/(A*x^2+B*x+C),start=list(A=a0,B=b0,C=c0))
    a <- summary(m1)$parameters[1][1]
    b <- summary(m1)$parameters[2][1]
    c <- summary(m1)$parameters[3][1]
    dataframe$a[dataframe$lc==lev[i]][1] <- a
    dataframe$b[dataframe$lc==lev[i]][1] <- b
    dataframe$c[dataframe$lc==lev[i]][1] <- c
    dataframe$etrmax[dataframe$lc==lev[i]][1] <- 1 / (b + 2 * sqrt(a*c))
    dataframe$alpha[dataframe$lc==lev[i]][1] <- 1 / c
    dataframe$ek[dataframe$lc==lev[i]][1] <- (1 / (b + 2 * sqrt(a*c))) / (1/c)

    #model predictions
    if(predict_ul ==T){
      new.data<-data.frame(x=seq(0,max.par,10))
      pred<-predictNLS(m1, new.data) %>% as.data.frame
      pred$par<-new.data$x
      colnames(pred)<-c("fit","mean","sd","median","mad","lower","upper","par")
      pred.catch[,paste(current.sample, 'fit',sep='_')] <- pred$fit
      pred.catch[,paste(current.sample, 'upper',sep='_')] <- pred$upper
      pred.catch[,paste(current.sample, 'lower',sep='_')] <- pred$lower
      rm(pred)
    }

  }#end of i loop


  #make averages if a treatment column is provided
  if(!is.na(treatment)){
    dataframe[,treatment] %<>% factor

    #average variables
    av <- aggregate(dataframe$etr ~ dataframe[,par] + dataframe[,treatment], FUN=mean)
    colnames(av) <- c(par, treatment, 'etr')
    av$etr_sd <- aggregate(dataframe$etr ~ dataframe[,par] + dataframe[,treatment], FUN=sd)[,3]
    av$count <- aggregate(dataframe$etr ~ dataframe[,par] + dataframe[,treatment], FUN=length)[,3]
    av$etr_se <- av$etr_sd / sqrt(av$count)

    av$npq <- aggregate(dataframe$npq ~ dataframe[,par] + dataframe[,treatment], FUN=mean)[,3]
    av$npq_sd <- aggregate(dataframe$npq ~ dataframe[,par] + dataframe[,treatment], FUN=sd)[,3]
    av$npq_se <- av$npq_sd / sqrt(av$count)

    av$ynpq <- aggregate(dataframe$ynpq ~ dataframe[,par] + dataframe[,treatment], FUN=mean)[,3]
    av$ynpq_sd <- aggregate(dataframe$ynpq ~ dataframe[,par] + dataframe[,treatment], FUN=sd)[,3]
    av$ynpq_se <- av$ynpq_sd / sqrt(av$count)

    av$yno <- aggregate(dataframe$yno ~ dataframe[,par] + dataframe[,treatment], FUN=mean)[,3]
    av$yno_sd <- aggregate(dataframe$yno ~ dataframe[,par] + dataframe[,treatment], FUN=sd)[,3]
    av$yno_se <- av$yno_sd / sqrt(av$count)

    #fit cumulative e and p model per treatment
    av$etrmax <- NA
    av$alpha <- NA
    av$ek <- NA
    treat.catch <- data.frame(par = seq(0, max.par, 10))
    treat.lev <- levels(av[,treatment])
    for(j in 1:length(treat.lev)){
      t.dat <- subset(av, av[,treatment] == treat.lev[j])
      a0<-0.00003; b0<--0.003; c0<-4
      x <- t.dat$par; y <- t.dat$etr
      m1 <- nls(y~x/(A*x^2+B*x+C),start=list(A=a0,B=b0,C=c0))
      a <- summary(m1)$parameters[1][1]
      b <- summary(m1)$parameters[2][1]
      c <- summary(m1)$parameters[3][1]
      av$etrmax[av[,treatment]==treat.lev[j]][1] <- 1 / (b + 2 * sqrt(a*c))
      av$alpha[av[,treatment]==treat.lev[j]][1] <- 1 / c
      av$ek[av[,treatment]==treat.lev[j]][1] <- (1 / (b + 2 * sqrt(a*c))) / (1/c)

      new.data<-data.frame(x=seq(0,max.par,10))
      pred<-predictNLS(m1, new.data) %>% as.data.frame
      pred$par<-new.data$x
      colnames(pred)<-c("fit","mean","sd","median","mad","lower","upper","par")
      treat.catch[,paste('treatment', treat.lev[j], 'fit',sep='_')] <- pred$fit
      treat.catch[,paste('treatment', treat.lev[j], 'upper',sep='_')] <- pred$upper
      treat.catch[,paste('treatment', treat.lev[j], 'lower',sep='_')] <- pred$lower
      rm(pred)

    }#end of av etr e and p fitting


    #average parameters
    p.dat <- dataframe[dataframe[,par]==0,]
    av.param <- aggregate(p.dat$etrmax ~ p.dat[,treatment], FUN=mean)
    colnames(av.param) <- c(treatment, 'etrmax')
    av.param$etrmax_sd <- aggregate(p.dat$etrmax ~ p.dat[,treatment], FUN=sd)[,2]
    av.param$count <- aggregate(p.dat$etrmax ~ p.dat[,treatment], FUN=length)[,2]
    av.param$etrmax_se <- av.param$etrmax_sd / sqrt(av.param$count)

    av.param$alpha <- aggregate(p.dat$alpha ~ p.dat[,treatment], FUN=mean)[,2]
    av.param$alpha_sd <- aggregate(p.dat$alpha ~ p.dat[,treatment], FUN=sd)[,2]
    av.param$alpha_se <- av.param$alpha_sd / sqrt(av.param$count)

    av.param$ek <- aggregate(p.dat$ek ~ p.dat[,treatment], FUN=mean)[,2]
    av.param$ek_sd <- aggregate(p.dat$ek ~ p.dat[,treatment], FUN=sd)[,2]
    av.param$ek_se <- av.param$ek_sd / sqrt(av.param$count)

    av.param$fvfm <- aggregate(p.dat[,yield] ~ p.dat[,treatment], FUN=mean)[,2]
    av.param$fvfm_sd <- aggregate(p.dat[,yield] ~ p.dat[,treatment], FUN=sd)[,2]
    av.param$fvfm_se <- av.param$fvfm_sd / sqrt(av.param$count)

  }#end of treatment averaging if clause

#return results dataframe for individual light curves, predictions dataframe (if selected),
#and average variables and parameters if a treatment column is provided.

  rlc_results<-list()
  rlc_results$rlcs <- dataframe

  if(predict_ul == T){
   rlc_results$fits <- pred.catch
  }

  if(!is.na(treatment)){
   rlc_results$av_variables <- av
   rlc_results$av_params <- av.param
   rlc_results$av_fits <- treat.catch
  }

  rlc_results<<-rlc_results
  return(rlc_results)

}#end of function
