#NNDIST SUMMARY
#################################################################################################################
#calculate the summary stats of the distribution of nearest neighbour distances across multiple created landscapes of a given type/clustering level

#landscape.data = data of multiple created landscapes
#no.replicates = number of landscapes that were created
#n.patches = number of patches within those landscapes created

#NOTE: only enter landscapes of one type/level of clustering
#use output data from eg. create multiple landscapes function selected for one given landscape type
#############################################################################################################
distribution.nndists.function<-function(landscape.data, no.replicates, n.patches){
  for (i in 1:no.replicates) {
    nndists<-nndist(landscape.data$x.coord[1*i:n.patches*i], landscape.data$y.coord[1*i:n.patches*i])
    #look at distribution of nndists
    max<-max(nndists)
    third<-quantile(nndists, prob = 0.75)
    median<-median(nndists)
    mean<-mean(nndists)
    first<-quantile(nndists, prob = 0.25)
    min<-min(nndists)
    dists.sum<-data.frame(min, first, median, mean, third, max)
    if (i==1) {dist.data<-dists.sum}
    else{dist.data<-rbind(dist.data, dists.sum)}
  }
  dist.data
  
  min.min<-min(dist.data$min)
  min<-mean(dist.data$min)
  min.se<-sd(dist.data$min/sqrt(no.replicates))
  first<-mean(dist.data$first)
  first.se<-sd(dist.data$first/sqrt(no.replicates))
  median<-mean(dist.data$median)
  median.se<-sd(dist.data$median/sqrt(no.replicates))
  mean<-mean(dist.data$mean)
  mean.se<-sd(dist.data$mean/sqrt(no.replicates))
  third<-mean(dist.data$third)
  third.se<-sd(dist.data$third/sqrt(no.replicates))
  max<-mean(dist.data$max)
  max.se<-sd(dist.data$max/sqrt(no.replicates))
  max.max<-max(dist.data$max)
  dist.summary<-data.frame(min.min, min, min.se, first, first.se, median, median.se, mean, mean.se, third, third.se, max, max.se, max.max)
  return(dist.summary)
}
#END NNDIST SUMMARY FUNCTION
############################################################################################################
