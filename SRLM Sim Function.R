#F4: SIMULATION OF SRLM MODEL FUNCTION (METAPOPULATION SIZE AND TIME TO EXTINCTION)
#################################################################################################################
#simulates colinization and extinction of patches within a specified landscape by a species with a specified 
#average dispersal distance and ratio of extinction to colonization for a specified number of timesteps from 
#a specified abundance

#landscape = a dataframe created by the create landscape function which specifies patch locations, areas, 
            #and and an interpatch distance matrix
#a = 1/(average dispersal distance) for a species
#delta = the ratio of the extinction to colonization rate of a species
#timesteps = # number of iterations of colonizations and extinctions of patches
#p.initial = initial probability a patch is occupied by a species

#NOTE: delta can be factored out of this calculation by letting delta = 1
#########################################################################################
SRLM.sim<-function(landscape, a, delta, timesteps, p.initial){
  
  landscape<-landscape
  a=a
  col.rate=1
  ext.rate=delta
  timesteps=timesteps
  
  n.patches<-length(landscape$patch.ID)
  A<-landscape$A
  d<-as.matrix(landscape[,5:(dim(landscape)+4)])
  
  initial.occupied<-rbinom(n.patches, size=1, prob=p.initial)
  #draw from the binomial distribution with pr(occupied)=p.star to 
  #determine which patches are occupied initially starting at equilibrium
  
  occupied<-rep(NA, n.patches*timesteps); dim(occupied)<-c(timesteps, n.patches)
  occupied[1,]<-initial.occupied 
  #initially all patches have a probability of 1 they are occupied, 
  #therefore they are all occupied
  for (t in 2:timesteps){ #for all future timesteps in the simulation
    P<-occupied[t-1,] #the probability a patch is occupied is the probability 
    #it was occupied in the previous iteration of the simulation
    connectivity.net<-rep(NA, n.patches*n.patches); dim(connectivity.net)<-c(n.patches,n.patches)
    for (i in 1:n.patches){ #for each row in the matrix
      for (j in 1:n.patches){ #for each column in the matrix
        if (i==j){connectivity.net[i,j]<-0} #let the connectivity of a patch to itself be 0
        else{connectivity.net[i,j]<-A[j]*exp(-a*d[i,j])*P[j]} 
        #let the connectivity of other patches to that patch be Aje^(-adij)Pj
      }
    }
    S<-rowSums(connectivity.net) #The patch connectivity (S) = the sum of a patches 
    #connectivity to all other patches (don't have to exclude itself because it is simply 0)
    C<-rep(NA, n.patches)
    exp.col<-rep(NA,n.patches)
    E<-rep(NA, n.patches)
    exp.ext<-rep(NA, n.patches)
    for (i in 1:n.patches) { #for each patch
      C[i]<-S[i]/(S[i]+(1/col.rate)) #discrete time version of colonization for SRLM
      exp.col[i]<-C[i]*(1-P[i]) #probability we expect a patch will be colonized = 
      #pr(colonization)*pr(unnoccupied (1-P)) (therefore, pr(empty and colonized))
    } 
    real.col<-rbinom(n.patches, size=1, prob=exp.col) 
    #number of patches actually colonized, 
    #drawn from the normal distribution with the expected probability of colonization
    for (i in 1:n.patches) { #for each patch
      E[i]<-ext.rate/A[i]*(1-C[i]) #pr(extinction) is e/Ai*(1-Ci), 
      #rescue effect which is the probability that a pop establishes within the timestep
      exp.ext[i]<-E[i]*P[i] #probability we expect a patch will go extinct = 
      #pr(occupied)*pr(extinction) (therefore, the pr(occupied and goes extinct)
    } 
    exp.ext<-ifelse(exp.ext>1,1,exp.ext) #pr(extinction) is never greater than 1 
    real.ext<-rbinom(n.patches, size=1, prob=exp.ext) 
    #number of patches actually gone extinct, drawn from the normal distribution 
    #with the expected probability of extinction
    occupied[t,]<-P-real.ext+real.col #number of occupied patches = *****
    #print(head(occupied))
  } 
  metapop<-rowSums(occupied) #size of metapop (given by total number of patches occupied)
  #plot(1:timesteps, metapop, type="b") #plotting the size of metapop over time
  colMeans(occupied[500:timesteps,])
  eq.size<-mean(metapop[(timesteps-50):timesteps]) 
  #finds average metapop size over the last 50 timesteps (assuming that equilibrium is reached by 950 iterations this should = p.star of the current landscape)
  if (any(metapop == 0)) {
    time.extinct<-which.max(metapop < 1) #finds first metapop entry with 0 patches occupied
  } else { time.extinct<-NA }#time.extinct=NA if never hits 0
  data<-list(#metapop=metapop, #include metapop to check
    eq.size=eq.size, time.extinct=time.extinct)
  return(data) 
} 
#END SRLM FUNCTION
############################################################################################################################################################################
