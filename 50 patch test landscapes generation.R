rm(list=ls()) #clear the workspace

#LANDSCAPE SPECIFICS ***Define the specifics of the landscapes you wish to generate here***
#############################################################################################################
#############################################################################################################
n.patches=50 #define the number of patches desired in the landscape
landscape.limit=100 #define the extent of the landscape, provides the x and y lim
no.runs=500 #define the number of iterations for the clustering procedure (higher the number the more 
#clustered/uniform the landscape for clustered/uniform landscapes respectively)
no.replicates=10 #number of landscapes to be generated of each specified type
directory="C:/Users/abuga/Desktop/Metapopulation Manuscript Output" #define where the 
#output should go
landscape.types=c("regular", "random", "clustered") #Define the type or types of landscapes to be generated
############################################################################################################
############################################################################################################

#LOAD REQUIRED PAKAGES
############################################################################################################
#nearest neighbours function comes from spatstat package which is unavailable for the R studio installed 
#on the server
library("spatstat") #contains nndist function
library("ggplot2") #for plotting landscapes (provides an image of what the landscape looks like)
library("boot") #for bootstrapping 95% CI's

setwd("C:/Users/abuga/Desktop/Metapopulation-Simulation-Package-master") #access working directory that all metapop 
#functions are stored in
source("Create Landscape Function.r") #load create.landscape.function (fuction creates a single landscape)
source("Rep Create Landscape Function.r") #load rep.create.landscape.function (function creates multiple 
#landscapes creating a dataset for them)
source("NNDist Distribution Summary Function.r") #load nndist.summary function (function calculates summary
#stats of the distribution of Nearest Neighbour distances within a set of landscapes) 
source("Choose Alphas Function.r") #load choose alphas function (function chooses average dispersal 
#distances based on the distribution of nearest neighbour distances within the set of landscapes)
############################################################################################################

#GENERATE LANDSCAPES
############################################################################################################
landscapes.data<-rep.create.landscape(n.patches=n.patches, landscape.limit=landscape.limit, no.runs=no.runs, 
                                      no.replicates=no.replicates, directory=directory) 
#creates the landscapes
#############################################################################################################

#SUBSETTING DATA BY LANDSCAPE TYPE
##############################################################################################################
clustered.data<-subset(landscapes.data, landscape.type=="clustered") 
#creates seperate clustered landscapes created dataframe
random.data<-subset(landscapes.data, landscape.type=="random")
#creates seperate random landscapes created dataframe
regular.data<-subset(landscapes.data, landscape.type=="regular")
#creates seperate regular landscapes created dataframe
##############################################################################################################

#CALCULATING SUMMARY STATS ON THE DISTRIBUTION OF NEAREST NEIGHBOUR DISTANCES FOR EACH LANDSCAPE TYPE
##############################################################################################################
regular.nndists.dist<-distribution.nndists.function(regular.data, no.replicates, n.patches)
#calculates summary stats for the distribution of nearest neighbour distances within regular landscapes created
random.nndists.dist<-distribution.nndists.function(random.data, no.replicates, n.patches)
#calculates summary stats for the distribution of nearest neighbour distances within random landscapes created
clustered.nndists.dist<-distribution.nndists.function(clustered.data, no.replicates, n.patches)
#calculates summary stats for the distribution of nearest neighbour distances within clustered landscapes created
all.nndists.dists<-rbind(regular.nndists.dist, random.nndists.dist, clustered.nndists.dist) 
#puts these all into one table
all.nndists.dists$landscape.type<-landscape.types #adds the respective landscape types to this table
write.csv(all.nndists.dists, "nndists distribution for each landscape type.csv") #create a csv file of this 
#table of summary stats for the distribution of nearest neighbour distances within each landscape type
##############################################################################################################

#CHOOSING ALPHA VALUES FOR SIMULATION BASED ON EACH LANDSCAPE TYPE'S DISTRIBUTION OF NNDISTS
#Thus, levels of alpha represent diff species but of equivalent dispersal ability within the given landscape
#controls for the fact that clustering will have smaller nearest neighbour distances by way of the clustering 
#algorithm, while regular will have larger nndists
#Isolates for the effect of regularity in patch spacing
#############################################################################################################
#CHECK
regular.alphas<-choose.alphas(nndists.distribution=regular.nndists.dist, landscape.limit=landscape.limit)
random.alphas<-choose.alphas(nndists.distribution=random.nndists.dist, landscape.limit=landscape.limit)
clustered.alphas<-choose.alphas(nndists.distribution=clustered.nndists.dist, landscape.limit=landscape.limit)
#putting it all into one table
all.alphas<-rbind(regular.alphas, random.alphas, clustered.alphas) #creates a table of alpha values chosen 
#for each landscape type
all.alphas$landscape.type<-landscape.types #adds the respective landscape types to this dataframe
write.csv(all.alphas, "alphas.csv") #create a csv file with the values for alpha chosen for each landscape type
#############################################################################################################
