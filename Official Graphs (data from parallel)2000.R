rm(list=ls()) #clear the workspace

library("gridExtra")
library("grid")
library("lattice")

#nearest neighbours function comes from spatstat package which is unavailable for the R studio installed on the server
library("spatstat")
library("ggplot2")
#library("rcompanion")
#library("fitdistrplus")
library("boot")
library("ggthemes")

setwd("C:/Users/abuga/Desktop/Metapopulation-Simulation-Package-master")#access working directory that all my metapop functions are stored in
source("Plotdata Median Lm Function.r")
source("Plotdata Median P Function.r")
source("Plotdata Median t extinct Function.r")
source("Plotdata Median t to P Function.r")
source("Plotdata Median t to P by lm Function.r")
source("Plotdata Median t to P binning lm Function.r")
source("Plotdata Median tex Function.r")

library("RColorBrewer")
library("colorspace")
library("wesanderson") #SUPER NICE COLOUR SCHEMES!
library("ggthemes")
pal <- choose_palette()

#READING IN DATA
################################################################################################################
data.set<-read.csv("C:/Users/abuga/Desktop/Metapopulation Manuscript Output/50reps_PnLMnSimEq_3_50p.csv")
################################################################################################################

#ADD IN A SCALED TIME TO EQ COLUMN
#######################################################################################################
data.set$t.to.eq.by.lm.r<-data.set$time.to.eq.r/data.set$lambda.M.r
data.set$t.to.eq.by.lm.e<-data.set$time.to.eq.e/data.set$lambda.M.e
#######################################################################################################

#SUBSETTING THE DATA FOR A GIVEN LANDSCAPE TYPE
################################################################################################################
subset.landscape<-function(data, landscape.type){
data<-data[data$landscape.type==landscape.type,]
return(data)}
###########################################################################################################
#CHECK data<-subset.landscape(data.set, landscape.type)

#GET UNIQUE ALPHAS
##########################################################################################################
get.a<-function(data){
a<-unique(data$alpha)
return(a)}
################################################################################################################
#CHECK a<-get.a(data)

#CLEANING DATA
################################################################################################################
clean.data<-function(a, data){
data<-data[data$delta>1/100000,] #removing landscape data in which the only reason the species is persisting is because of the base colinization probablity
data<-data[complete.cases(data[ , 8]),]  #getting rid of rows with NA's in all but the time extinct column (based on values in column 8)
#exclude any data for which there are less than 20 replicates
df2<-NULL
for (i in 1:length(a)){
  input.data<-data[data$alpha==a[i],]
  df<-input.data[length(unique(input.data$rep.no[input.data$alpha==a[i]]))>20]
  df2<-rbind(df2, df)}
data<-df2}
################################################################################################################
#CHECK data<-clean.data(a, data)

#GET DISTRIBUTION OF DELTAS
###########################################################################################################
get.deltas<-function(a, data, landscape.type){
deltas<-rep(NA, 6*length(a)); dim(deltas)<-c(length(a),6)
for (i in 1:length(a)){
  a.group<-data[data$alpha==a[i],]
  deltas[i,]<-(summary(a.group$delta))
}
deltas<-as.data.frame(deltas)
deltas$a<-a
deltas$avg.disp<-c("1/8x Avg Min Nearest Neighbour Distance","1/4x Avg Min Nearest Neighbour Distance", "1/2x Avg Min Nearest Neighbour Distance", "1x Avg Min Nearest Neighbour Distance", "2x Avg Min Nearest Neighbour Distance", "4x Avg Min Nearest Neighbour Distance", "8x Avg Min Nearest Neighbour Distance", "Global Dispersal")
colnames(deltas)<-c("Min", "1st.Q", "Med", "Mean", "3rd.Q", "Max", "a", "avg.disp")
deltas
setwd("C:/Users/abuga/Desktop/Manuscript/Data and Figures/50 patches")
write.csv(deltas, file = paste0(landscape.type, "50p deltas2000 sim new.csv"))
return(deltas)
}
########################################################################################################################################
#CHECK deltas<-get.deltas(a, data, landscape.type)


#LM GRAPH
LMplot.data<-function(deltas, data){
deltas<-na.omit(deltas)
plot.data<-NULL
for (i in 1:length(deltas$a)){
  df3<-plotdata.median.Lm(data[data$alpha==deltas$a[i], ], scaled.lm=20)
  df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
  plot.data<-rbind(plot.data, df3)
}
  return(plot.data)
}
#CHECK plot.data<-LMplot.data(deltas, data)

plot.LM<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data) + theme_classic()
        + geom_ribbon(aes(x=percent.loss, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = percent.loss, y = destruction.median, group=(1/a), colour=1/(a)), size=1) 
        + geom_line(aes(x = percent.loss, y = degradation.lambda.M, colour="Degradation"), size=1,  color="red", linetype="dashed") 
        + geom_hline(aes(yintercept=0.05), color="black", size=1) 
        + labs(x = "Percent Habitat Loss", y = "Persistence Capacity Final/Intial", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", 
                                trans="log", 
                                breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), 
                                labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal")),
                                guide="legend") #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=1))
        #+ theme(legend.position="none") #to remove legend
        + theme(legend.position = c(0.9, 0.9))
        + guides(color=guide_legend(keywidth=0.1, keyheight=0.2, default.unit="inch", override.aes = list(size = 4)))
        )
}
#CHECK plot.LM(plot.data, landscape.type)

landscape.type="random"
title="More Uniform"
graph.LM<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-LMplot.data(deltas, data)
  plot.LM(plot.data, landscape.type, a, title)
}

plot1<-graph.LM(data.set, "regular", "More Uniform")
plot2<-graph.LM(data.set, "random", "Random")
plot3<-graph.LM(data.set, "clustered", "More Clustered")

#grid.arrange(plot1, plot2, plot3, ncol=2)

#generating all required initial persistence plots for every alpha value and landscape type
##############################################################################################################a=c(0.067, 0.105, 0.429, 0.571, 0.96, 100)
plot.Pi<-function(plot.data, landscape.type, title){
  plot.data$degradation.p.initial<-(1-plot.data$degradation.lambda.M)
  plot.data$destruction.p.initial<-(1-plot.data$destruction.median)
  print(ggplot(data=plot.data) + theme_classic() +
          geom_line(aes(x = percent.loss, y = degradation.p.initial, colour="Degradation"), size=1,  color="red", linetype="dashed") +
          #geom_ribbon(aes(x=percent.loss, ymin=0, ymax=degradation.p.initial, group=1, fill="Degradation"), alpha=0.3) +
          geom_ribbon(aes(x=percent.loss, ymin=0, ymax=destruction.p.initial, group=(1/a), fill=(1/a)), alpha=0.3) +
          labs(x = "Percent Habitat Loss", y = "Porportion of initial habitat occupied to persist", 
               title = paste0(title, " Landscapes"))
        #+ scale_fill_manual(values=rev(pal(8)), labels=c(rev(deltas$avg.disp), "Degradation"))) 
         + scale_fill_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
  + theme(text = element_text(size=15))
  + theme(legend.text=element_text(size=12))
  + theme(legend.key.size = unit(1.5, "cm"))
  + theme(legend.title = element_text(size=17, vjust=20))
  + theme(legend.position="none") #to remove legend
  )
  }

graph.Pi<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-LMplot.data(deltas, data)
  plot.Pi(plot.data, landscape.type, title)
}

plot1<-graph.Pi(data.set, "regular", "More Uniform")
plot2<-graph.Pi(data.set, "random", "Random")
plot3<-graph.Pi(data.set, "clustered", "More Clustered")

#P* GRAPH
Pplot.data<-function(deltas, data){
  deltas<-na.omit(deltas)
  plot.data<-NULL
  for (i in 1:length(deltas$a)){
    df3<-plotdata.median.p(data[data$alpha==deltas$a[i], ])
    df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
    plot.data<-rbind(plot.data, df3)
  }
  return(plot.data)
}
#CHECK plot.data<-Pplot.data(deltas, data)

plot.P<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data) + theme_classic()
        #+ geom_ribbon(aes(x=percent.loss, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        #+ geom_line(aes(x = percent.loss, y = destruction.median, group=(1/a), colour=1/(a)), size=1, linetype="dashed") 
        + geom_ribbon(aes(x=percent.loss, ymin= degradation.lower.CI, ymax=degradation.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = percent.loss, y = degradation.median, group=(1/a), colour=1/(a)), size=1, linetype="dashed") 
        + labs(x = "Percent Habitat Loss", y = "Avg. Patch Occupancy (P*)", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=20))
        + theme(legend.position="none") #to remove legend
        + ylim(0,1) + xlim(0,1)
)
}
#CHECK plot.P(plot.data, landscape.type, a, title)

graph.P<-function(data.set, landscape.type, title){
data<-subset.landscape(data.set, landscape.type)
a<-get.a(data)
data<-clean.data(a, data)
deltas<-get.deltas(a, data, landscape.type)
plot.data<-Pplot.data(deltas, data)
plot.P(plot.data, landscape.type, a, title)
}

plot1<-graph.P(data.set, "regular", "More Uniform")
plot2<-graph.P(data.set, "random", "Random")
plot3<-graph.P(data.set, "clustered", "More Clustered")

####
#P1000
simplot.data<-function(deltas, data){
  deltas<-na.omit(deltas)
  plot.data<-NULL
  for (i in 1:length(deltas$a)){
    df3<-plotdata.median.sim.eq(data[data$alpha==deltas$a[i], ])
    df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
    plot.data<-rbind(plot.data, df3)
  }
  return(plot.data)
}
#CHECK plot.data<-Pplot.data(deltas, data)

plot.sim<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data) + theme_classic()
        + geom_ribbon(aes(x=percent.loss, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = percent.loss, y = destruction.median, group=(1/a), colour=1/(a)), size=1) 
        #+ geom_ribbon(aes(x=percent.loss, ymin= degradation.lower.CI, ymax=degradation.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        #+ geom_line(aes(x = percent.loss, y = degradation.median, group=(1/a), colour=1/(a)), size=1) 
        + labs(x = "Percent Habitat Loss", y = "Avg. Patch Occupancy", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=20))
        + xlim(0,1) + ylim(0,1)
        + theme(legend.position="none") #to remove legend
  )
}
#CHECK plot.P(plot.data, landscape.type, a, title)

landscape.type="regular"
title="More Uniform"
graph.sim<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-simplot.data(deltas, data)
  plot.sim(plot.data, landscape.type, a, title)
}

plot1<-graph.sim(data.set, "regular", "More Uniform")
plot2<-graph.sim(data.set, "random", "Random")
plot3<-graph.sim(data.set, "clustered", "More Clustered")

####
#time to expected eq GRAPH
tplot.data<-function(deltas, data){
  deltas<-na.omit(deltas)
  plot.data<-NULL
  for (i in 1:length(deltas$a)){
    df3<-plotdata.median.t.eq(data[data$alpha==deltas$a[i], ])
    df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
    plot.data<-rbind(plot.data, df3)
  }
  return(plot.data)
}
#CHECK plot.data<-tplot.data(deltas, data)

plot.t<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data[complete.cases(plot.data[,1:7]),]) + theme_classic()
        #+ geom_ribbon(aes(x=percent.loss, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        #+ geom_line(aes(x = percent.loss, y = destruction.median, group=(1/a), colour=1/(a)), size=1) 
        + geom_ribbon(aes(x=percent.loss, ymin= degradation.lower.CI, ymax=degradation.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = percent.loss, y = degradation.median, group=(1/a), colour=1/(a)), size=1, linetype="dashed") 
        + labs(x = "Percent Habitat Loss", y = "time to P1000", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=20))
        + xlim(0,1) + ylim(0,100)
        + theme(legend.position="none") #to remove legend
  )
}
#CHECK plot.P(plot.data, landscape.type, a, title)

landscape.type="clustered"
title="More Clustered"
graph.t<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-tplot.data(deltas, data)
  plot.t(plot.data, landscape.type, a, title)
}

plot1<-graph.t(data.set, "regular", "More Uniform")
plot2<-graph.t(data.set, "random", "Random")
plot3<-graph.t(data.set, "clustered", "More Clustered")

####
#time extinct only GRAPH ***95% CI's not working but should be able to fix
texplot.data<-function(deltas, data){
  data<-data[data$sim.eq.size.r==0,] #subset the data to only include entries where extinction occured for destruction
  deltas<-na.omit(deltas)
  plot.data<-NULL
  for (i in 1:length(deltas$a)){
    df3<-plotdata.median.tex(data[data$alpha==deltas$a[i], ])
    df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
    plot.data<-rbind(plot.data, df3)
  }
  return(plot.data)
}
#CHECK plot.data<-texplot.data(deltas, data)

plot.tex<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data[complete.cases(plot.data[,1:7]),]) + theme_classic()
        #+ geom_ribbon(aes(x=percent.loss, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = percent.loss, y = destruction.median, group=(1/a), colour=1/(a)), size=1) 
        #+ geom_ribbon(aes(x=percent.loss, ymin= degradation.lower.CI, ymax=degradation.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        #+ geom_line(aes(x = percent.loss, y = degradation.median, group=(1/a), colour=1/(a)), size=1, linetype="dashed") 
        + labs(x = "Percent Habitat Loss", y = "time extinct", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=20))
        #+ xlim(0,1) 
        + ylim(0,1000)
        + theme(legend.position="none") #to remove legend
  )
}
#CHECK plot.P(plot.data, landscape.type, a, title)

landscape.type="clustered"
title="More Clustered"
graph.tex<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-texplot.data(deltas, data)
  plot.tex(plot.data, landscape.type, a, title)
}

plot1<-graph.tex(data.set, "regular", "More Uniform")
plot2<-graph.tex(data.set, "random", "Random")
plot3<-graph.tex(data.set, "clustered", "More Clustered")









####
#time extinct only GRAPH
texplot.data<-function(deltas, data){
  data<-data[data$sim.eq.size.e==0,] #subset the data to only include entries where extinction occured for destruction
  deltas<-na.omit(deltas)
  plot.data<-NULL
  for (i in 1:length(deltas$a)){
    df3<-plotdata.median.tex(data[data$alpha==deltas$a[i], ])
    df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
    plot.data<-rbind(plot.data, df3)
  }
  return(plot.data)
}
#CHECK 
plot.data<-texplot.data(deltas, data)

plot.tex<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data[complete.cases(plot.data[,1:7]),]) + theme_minimal()
        #+ geom_ribbon(aes(x=percent.loss, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        #+ geom_line(aes(x = percent.loss, y = destruction.median, group=(1/a), colour=1/(a)), size=1) 
        + geom_ribbon(aes(x=percent.loss, ymin= degradation.lower.CI, ymax=degradation.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = percent.loss, y = degradation.median, group=(1/a), colour=1/(a)), size=1) 
        + labs(x = "Percent Habitat Loss", y = "time extinct", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=20))
        + xlim(0,1) + ylim(0,1000)
        + theme(legend.position="none") #to remove legend
  )
}
#CHECK plot.P(plot.data, landscape.type, a, title)

landscape.type="regular"
title="More Uniform"
graph.tex<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-texplot.data(deltas, data)
  plot.tex(plot.data, landscape.type, a, title)
}

plot1<-graph.tex(data.set, "regular", "More Uniform")
plot2<-graph.tex(data.set, "random", "Random")
plot3<-graph.tex(data.set, "clustered", "More Clustered")


#pal <- wes_palette("Zissou", 8, type = "continuous") #pretty
#wes_palette("Zissou")
####################################################################

####
#time to expected eq / lm GRAPH
t.by.lm.plot.data<-function(deltas, data){
  deltas<-na.omit(deltas)
  plot.data<-NULL
  for (i in 1:length(deltas$a)){
    df3<-plotdata.median.t.eq.by.lm(data[data$alpha==deltas$a[i], ])
    df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
    plot.data<-rbind(plot.data, df3)
  }
  return(plot.data)
}
#CHECK plot.data<-tplot.data(deltas, data)

plot.t.by.lm<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data[complete.cases(plot.data[,1:7]),]) + theme_classic()
        #+ geom_ribbon(aes(x=percent.loss, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        #+ geom_line(aes(x = percent.loss, y = destruction.median, group=(1/a), colour=1/(a)), size=1) 
        + geom_ribbon(aes(x=percent.loss, ymin= degradation.lower.CI, ymax=degradation.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = percent.loss, y = degradation.median, group=(1/a), colour=1/(a)), size=1, linetype="dashed") 
        + labs(x = "Percent Habitat Loss", y = "time to P", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=20))
        + xlim(0,1)
        + theme(legend.position="none") #to remove legend
  )
}
#CHECK plot.P(plot.data, landscape.type, a, title)

landscape.type="clustered"
title="More Clustered"
graph.t.by.lm<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-t.by.lm.plot.data(deltas, data)
  plot.t.by.lm(plot.data, landscape.type, a, title)
}

plot1<-graph.t.by.lm(data.set, "regular", "More Uniform")
plot2<-graph.t.by.lm(data.set, "random", "Random")
plot3<-graph.t.by.lm(data.set, "clustered", "More Clustered")

####
#time to expected eq vs. lm GRAPH
t.vs.lm.plot.data<-function(deltas, data){
  deltas<-na.omit(deltas)
  plot.data<-NULL
  for (i in 1:length(deltas$a)){
    df3<-plotdata.median.t.eq.by.lm(data[data$alpha==deltas$a[i], ])
    df3<-data.frame(a=rep(deltas$a[i], nrow(df3)), avg.disp=rep(deltas$avg.disp[i], nrow(df3)), df3)
    plot.data<-rbind(plot.data, df3)
  }
  return(plot.data)
}
#CHECK plot.data<-tplot.data(deltas, data)


plot.t.vs.lm<-function(plot.data, landscape.type, a, title){
  print(ggplot(data=plot.data[complete.cases(plot.data[,1:7]),]) + theme_classic()
        #+ geom_ribbon(aes(x=lm, ymin= destruction.lower.CI, ymax=destruction.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        #+ geom_line(aes(x = lm, y = destruction.median, group=(1/a), colour=1/(a)), size=1) 
        + geom_ribbon(aes(x=lm, ymin= degradation.lower.CI, ymax=degradation.upper.CI, group=(1/a), fill=factor(1/a)), alpha=0.3)
        + geom_line(aes(x = lm, y = degradation.median, group=(1/a), colour=1/(a)), size=1, linetype="dashed") 
        + labs(x = "Lambda M Final/Intial", y = "time to P", 
               title = paste0(title, " Landscapes"))
        + scale_colour_gradient("Dispersal Ability",low = pal(8),high = "light grey", trans="log", breaks = c(100, 8, 4, 2, 1, 0.5, 0.25, 0.125), labels = rev(c("1/8x", "1/4x", "1/2x", "Avg Min Nearest Neighbour", "2x", "4x", "8x", "Global Dispersal"))) #
        + scale_fill_manual(values=pal(8), guide="none")
        + theme(text = element_text(size=15))
        + theme(legend.text=element_text(size=12))
        + theme(legend.key.size = unit(1.5, "cm"))
        + theme(legend.title = element_text(size=17, vjust=20))
        + xlim(0,1) + ylim(0,1000)
        + theme(legend.position="none") #to remove legend
  )
}
#CHECK plot.t.by.lm(plot.data, landscape.type, a, title)
##########<-----------------------------------Error: Discrete value supplied to a continuous scale

landscape.type="clustered"
title="More Clustered"
graph.t.vs.lm<-function(data.set, landscape.type, title){
  data<-subset.landscape(data.set, landscape.type)
  a<-get.a(data)
  data<-clean.data(a, data)
  deltas<-get.deltas(a, data, landscape.type)
  plot.data<-t.vs.lm.plot.data(deltas, data)
  plot.t.vs.lm(plot.data, landscape.type, a, title)
}

plot1<-graph.t.vs.lm(data.set, "regular", "More Uniform")
plot2<-graph.t.vs.lm(data.set, "random", "Random")
plot3<-graph.t.vs.lm(data.set, "clustered", "More Clustered")

