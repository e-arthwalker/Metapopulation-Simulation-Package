rm(list=ls()) #clear workspace

setwd("C:/Users/Administrator/Desktop/Manuscript/Functions") #access working directory that all my metapop functions are stored in
source("Create Landscape Function.r") #load create.landscape.function
source("Lambda M Function.r") #load the calculate lamda.M.function
source("Pstar Function.r") #load the persistence.function
source("time to eq3.r") #load the SRLM.sim function
source("Destroy and Degrade a Landscape Function2 (P and LM and SimEq)2.r") #destroy.vs.degrade.function only calculating Pstar and Lambda.M 
source("Degrade and Destroy Multiple Landscapes Function.r") #replicates landscape creation and destroy vs. degrade specified metapop parameters a specified number of times
source("Vary Alpha Destroy and Degrade Function for parallel.r") #provides replicates for a range of alphas across a landscape type




#CHECK landscape creation function
landscape<-create.landscape(n.patches=50, landscape.limit=100, landscape.type="clustered", no.runs=500)
head(landscape)
plot(landscape$y.coord, landscape$x.coord)

#CHECK lambda.M function
lambda.M<-lambda.M.function(landscape=landscape, a=0.1, delta=1)
lambda.M

#CHECK scaling of lambda.M
delta<-(lambda.M=lambda.M.function(landscape=landscape, a=0.1, delta=1))/20
delta
lambda.M.function(landscape=landscape, a=0.1, delta=delta)

#CHECK pstar function
p.star<-pstar.function(landscape=landscape, a=0.1, delta=delta, iterations=1000)
p.star

#CHECK SRLM function
data<-SRLM.sim(landscape=landscape, a=100, delta=delta, timesteps=1000, p.initial=p.star, avg.p=sum(p.star)/50)
data

#CHECK destroy.vs.degrade.p.n.Lm function
regular.landscape<-create.landscape(n.patches=50, landscape.limit=100, landscape.type="regular", no.runs=500)
#scaling delta starting with initial lamda.M of 20
delta<-20/(lambda.M=lambda.M.function(landscape=regular.landscape, a=100, delta=1))
regular.output<-destroy.vs.degrade(landscape=regular.landscape, a=100, delta=delta)
landscape.type<-rep("regular", 2)
regular.output<-data.frame(landscape.type, regular.output)
head(regular.output)
plot(regular.output$percent.habitatloss, regular.output$lambda.M.e)
plot(regular.output$percent.habitatloss, regular.output$lambda.M.r)
regular.output
data
