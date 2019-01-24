#Calculating data for plots (binning destruction and calculating median lamda.M for percent loss
#, using median because it reflects that the data distribution is heavily skewed,
#, confidence intervals claculated using bootstrapping)
##############################################################################################################
plotdata.totalext<-function(input.data){
  input.data<-input.data
  input.data$binned.percent.loss <- cut(input.data$percent.habitatloss, c(-Inf,seq(0, 1, 0.05)), labels=seq(0,1,0.05))
  percent.loss<-levels(droplevels(input.data$binned.percent.loss))
  bins<-percent.loss
  destruction.totalext<-rep(NA, length(bins))
  degradation.totalext<-rep(NA, length(bins))
  for (i in 1:length(bins)) { #for each bin
    x<-subset(input.data, input.data$binned.percent.loss==bins[i])
    for (j in 1:2000){ #for each metapop in each of those bins
      #destruction
      if (is.na(destruction.totalext[i]) == TRUE) { #if no metapops have been recorded as going extinct
        if (x[x[x$sim.eq.size.r==0,] & x[x$rep.no==j,]]) { #if the data has a row where an extinction occured for this metapop
          destruction.totalext[i] <- 1} #record an metapop extinction as having occured within that bin
      }
      else if (x[x[x$sim.eq.size.r==0,] & x[x$rep.no==j,]]) { #otherwise, if the data has a row where an extinction occured for this metapop
        destruction.totalext[i] <- destruction.totalext[i] + 1} #record another metapop as having gone extinct
      #degradation
      if (is.na(degradation.totalext[i]) == TRUE) { #if no metapops have been recorded as going extinct
        if (x[x[x$sim.eq.size.e==0,] & x[x$rep.no==j,]]) { #if the data has a row where an extinction occured for that run
          degradation.totalext[i] <- 1} #record an metapop extinction as having occured within that bin
      }
      else if (x[x[x$sim.eq.size.e==0,] & x[x$rep.no==j,]]) { #otherwise, if the data has a row where an extinction occured for this metapop
        degradation.totalext[i] <- degradation.totalext[i] + 1} #record another metapop as having gone extinct
    }
  plot.data<-data.frame(percent.loss=seq(0,1,0.05), destruction.totalext, degradation.totalext)}
  return(plot.data)
}
##############################################################################################################a=c(0.067, 0.105, 0.429, 0.571, 0.96, 100)
#CHECK
