#CHOOSING ALPHA VALUES 
#################################################################################################################
#choose alpha values for simulation based on each landscape type's distribution of nndists
#thus, levels of alpha represent diff species but of equivalent dispersal ability within the given landscape
#controls for the fact that clustering will have smaller nearest neighbour distances by way of the clustering algorithm, while regular will have larger nndists
#isolates for the effect of regularity in patch spacing

#nndists.distribution = data.frame of summary stats of nndists for all three landscape types 
#(aka output from NNDISTS Distribution Summary function)
#landscape.limit = size of landscape (max x and y coordinate that was possible in landscape generation)

#############################################################################################################
choose.alphas<-function(nndists.distribution, landscape.limit){
  a.16min<-nndists.distribution[2]/1000 #average dispersal distance 16x smaller than mean(min) nndist
  a.8min<-nndists.distribution[2]/100 #average dispersal distance 8x smaller than mean(min) nndist 
  a.4min<-nndists.distribution[2]/10 #average dispersal distance 4x smaller than mean(min) nndist
  a.2min<-nndists.distribution[2] #average dispersal distance = mean(min) nndist
  min.a.first<-nndists.distribution[2]+((nndists.distribution[4]-nndists.distribution[2])/2) #mean(min) nndist < average dispersal distance < mean(1stq) nndist
  first.a.median<-nndists.distribution[4]+((nndists.distribution[6]-nndists.distribution[4])/2) #mean(1stq) nndist < average dispersal distance < mean(median) nndist
  median.a.third<-nndists.distribution[6]+((nndists.distribution[10]-nndists.distribution[6])/2) #mean(median) nndist < average dispersal distance < mean(3rd) nndist
  third.a.max<-nndists.distribution[10]+((nndists.distribution[12]-nndists.distribution[10])/2) #mean(3rd) nndist < average dispersal distance < mean(max) nndist
  a.limit<-landscape.limit #average dispersal distance = landscape limit (global dispersal)
  alphas<-data.frame(a.8min, a.4min, a.2min, min.a.first, first.a.median, median.a.third, third.a.max, a.limit)
  alphas<-1/alphas
  return(alphas)
}
###########################################################################################################
#CHECK
#alphas<-choose.alphas(nndists.distribution=clustered.nndists.dist, landscape.limit=100)
#alphas