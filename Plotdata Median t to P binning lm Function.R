#Calculating data for plots (binning destruction and calculating median lamda.M for percent loss
#, using median because it reflects that the data distribution is heavily skewed,
#, confidence intervals claculated using bootstrapping)
##############################################################################################################
plotdata.median.t.eq.by.lm<-function(input.data){
  input.data<-input.data
  input.data$binned.lm.r <- cut(input.data$lambda.M.r/20, c(-Inf,seq(0, 1, 0.05)), labels=seq(0,1,0.05))
  lm.r<-levels(droplevels(input.data$binned.lm.r))
  bins<-lm.r
  destruction.median<-rep(NA, length(bins))
  destruction.upper.CI<-rep(NA, length(bins))
  destruction.lower.CI<-rep(NA, length(bins))
  for (i in 1:length(bins)) {
    x<-subset(input.data, input.data$binned.lm.r==bins[i])
    #destruction
    destruction.median[i] <- median(x$time.to.eq.r*x$lambda.M.r/20)
    bootobject <- boot(x$time.to.eq.r*x$lambda.M.r/20, function(u,j) median(u[j]), R=1000)
    CI<-boot.ci(bootobject, conf=0.95, type="basic")
    if(is.null(CI)==TRUE){destruction.upper.CI[i]<-destruction.median[i]
    destruction.lower.CI[i]<-destruction.median[i]}
    else{
      destruction.upper.CI[i]<-CI$basic[5]
      destruction.lower.CI[i]<-CI$basic[4]}}
    #degradation
    input.data$binned.lm.e <- cut(input.data$lambda.M.e/20, c(-Inf, seq(0, 1, 0.05)), labels=seq(0,1,0.05))
    lm.e<-levels(droplevels(input.data$binned.lm.e))
    bins<-lm.e
    degradation.upper.CI<-rep(NA, length(bins))
    degradation.lower.CI<-rep(NA, length(bins))
    degradation.median<-rep(NA, length(bins))
    for (i in 1:length(bins)) {
      x<-subset(input.data, input.data$binned.lm.e==bins[i])
    degradation.median[i] <- median(x$time.to.eq.e*x$lambda.M.e/20)
    bootobject <- boot(x$time.to.eq.e*x$lambda.M.e/20, function(u,j) median(u[j]), R=1000)
    CI<-boot.ci(bootobject, conf=0.95, type="basic")
  if(is.null(CI)==TRUE){degradation.upper.CI[i]<-degradation.median[i]
  degradation.lower.CI[i]<-degradation.median[i]}
  else{
    degradation.upper.CI[i]<-CI$basic[5]
    degradation.lower.CI[i]<-CI$basic[4]}}
    
  plot.data<-data.frame(lm=bins, destruction.median, destruction.upper.CI, destruction.lower.CI, degradation.median, degradation.upper.CI, degradation.lower.CI)
  return(plot.data)
}
##############################################################################################################a=c(0.067, 0.105, 0.429, 0.571, 0.96, 100)
#CHECK
