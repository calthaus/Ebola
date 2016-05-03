###################################################################
#   Monte Carlo simulations *WITHOUT* Convalescent Transmission   #
#                                                                 #
#  Accompanying:                                                  #
#                                                                 #
#  Potential impact of sexual transmission on Ebola virus         #
#  epidemiology: Sierra Leone as a case study.                    #
#  JL Abbate, C-L Murall, H Richner, CL Althaus                   #
#  2016 PLoS Neglected Tropical Diseases                          #
#                                                                 #
#  PrePrint available bioRxiv http://dx.doi.org/10.1101/031880    #
#  December 11, 2015                                              #        
###################################################################

# required libraries
library(MASS)
library(plotrix)

###################################################################

#    Sierra Leone Ebola incidence data based on raw case counts from:
#    WHO Ebola data and statistics: Data on new cases per epi week for Sierra Leone
#    18 Nov 2015 http://apps.who.int/gho/data/node.ebola-sitrep.ebola-country-SLE-20151118

thefile <- read.csv(file="data_WHO_SierraLeone_updated.csv",header=T)

###################################################################

set.seed(451) 
runs <- 1000  
tmax<-5000   

#create an array for the results
results <- array(NA,c(runs,tmax,10)) #pop copies here
names(results)<-c("time","S","E","I","R","D","beta","tot","currcases","disevents")


#create an array for summary results
outs<-array(NA,c(runs,7))

colnames(outs)<-c("test","R0", "TotalNumCases", "NoMoreEIs","DayLastCase",
                  "EpidemicPeakNew","DayEpidemicPeakNew")


#give this run a descriptive name. what is being tested?
test <- "base parameters"

#set parameter values
R0 <- 2.133485        # fitted

#store initial parameter values
outs[,1]=test
outs[,2]=R0

#loop on the number of runs
for(i in c(1:runs))
{
  print(i)                              #track progress
  
  #initialise values
  N0=6.316e6         
  R0=R0
  R1=0.67                               # lowest R0 with control measures, fitted
  InitialNumCases<-1  
  gamma=1/3.9       
  sigma = 1/11.4
  f = 0.69
  k = 0.01103038                        # control measure transmission decay, fitted 
  tau1 = 51.00493200                    # day control measures are implemented, fitted
  beta0=R0*gamma/N0
  beta1=R1*gamma/N0
  beta=beta0
  tot <- 1
  S <- N0-InitialNumCases
  E <-0
  I <- InitialNumCases
  C <-0
  R<-0
  D<-0
  t <- 1
  currcases<-E+I    
  disevents<-I                         # current number of new acute disesae cases (E -> I) 
  simtime<-t
  N<-N0
  
  #initialise results vector
  pop <- c(t,S,E,I,R,D,beta,tot,currcases,disevents)
  names(pop)<-c("time","S","E","I","R","D","beta","tot","currcases","disevents")
  results[i,t,]<-pop
  
  #loop until max time reached or extinction
  while((t<=tmax)&(I+E>0))  
  { 
    N<-S+E+I+R
    beta <- ifelse(simtime<tau1,beta0,beta1 + (beta0-beta1)*exp(-k*(simtime-tau1)))
    
    #determine when next event happens (drawn from an exponential distribution)
    events<-c(beta*S*I,sigma*E,gamma*f*I,(1-f)*gamma*I)
    sevents<-sum(events)
    cevents<-cumsum(events)
    deltat <- 1/sevents*log(1/runif(1))
    
    #determine which event happens
    
    #transmission events S->E
    rand <- runif(1,0,sevents)
    if(rand<cevents[1])
    {
      S=S-1;
      E=E+1;
    }
    
    #incubation period E->I
    else if(rand<cevents[2])
    {
      E=E-1;
      I=I+1;
      tot=tot+1;
    }
    
    #death of infecteds I->D
    else if(rand<cevents[3])
    {
      I=I-1;
      D=D+1;
    }
    
    #recovery of infecteds I->R
    else if(rand<=cevents[4])
    {
      I=I-1;
      R=R+1;
    }
    
    currcases<-E+I          
    
    #update time and store results
    simtime<-simtime+deltat
    if(floor(simtime)>t){t<- floor(simtime)}
    pop <- c(t,S,E,I,R,D,beta,tot,currcases,disevents) #c("time","S","E","I","R","D","beta","tot","currcases","disevents")
    results[i,t,]<-pop
  }
  
  # pop <- c([1]t,[2]S,[3]E,[4]I,[5]R,[6]D,[7]beta,[8]tot,
  #          [9]currcases,[10]disevents)
  
  colnames(outs)<-c("test","R0", "TotalNumCases", "NoMoreEIs","DayLastCase",
                    "EpidemicPeakNew","DayEpidemicPeakNew")
  
  
  #              Cumulative Cases
  outs[i,3]<-max(results[i,,8],na.rm=TRUE)
  
  #              End of Epidemic  (No more E or I individuals)
  outs[i,4]<-max(which(results[i,,9]==0)) #latest day no. current cases (E+I) = 0
  
  #              Day of Last Symptomatic Case
  if(max(results[i,,4],na.rm=TRUE)==0){outs[i,5]<-1}else{
  outs[i,5]<-max(which(results[i,,4]>0))} #Last day a symptomatic case was seen
  
  #              New Disease Events (incidence, E --> I) per time step (t)
  for (a in 2:tmax){
    if(is.na(results[i,a,8])){results[i,a,8]<-results[i,a-1,8]};
    results[i,a,10]<-results[i,a,8]-results[i,a-1,8]
  }
  
  #              Peak Epidemic Size (max Incidence)
  outs[i,6]<-max(results[i,,10],na.rm=TRUE)   #the peak of the epidemic (max #disevents at any one time)    
  
  #              Day Epidemic Peaks (new cases)
  outs[i,7]<-min(which(results[i,,10]==as.integer(outs[i,6])),na.rm = T)
}

write.table(outs,"Stoch_noSTI_summaries_1000.csv")

# plots ######################################################

par(mfrow=c(1,1), lwd=1)

dayIndexCase<-109  # April 23, 2014
weekIndexCase<-17  # week ends April 27, 2014

#######################################
# Cumulative number of cases
#######################################
 pdf("stoch_1000_cumul_noSTI.pdf", width = 6, height = 4)

xmax<-max(as.integer(outs[,4],rm.na=TRUE)) 
ymax<-max(as.integer(outs[,3],rm.na=TRUE))
xmaxdays<-1500
spaces<-xmaxdays/5
plot(1,xlim=c(1,xmaxdays),ylim=c(0,30000),xlab="Days after index case",ylab="Cumulative number of cases",pch=20,main="Without STI")

for(i in 1:runs)
{ #cumulative number of cases
  lines(results[i,,1],results[i,,8],col="red",lwd=0.25)  
}

# zeromeans #######################
meanCumuls<-colMeans(results[,,8],na.rm = TRUE)
SEMcumul<-apply(results[,,8],2,std.error,na.rm=TRUE)
meanplus<-meanCumuls+SEMcumul
meanminus<-meanCumuls-SEMcumul
lines(c(1:tmax),meanplus,col="blue",lwd=0.75)  #stoch SEM
lines(c(1:tmax),meanminus,col="blue",lwd=0.75)  #stoch SEM
points(c(1:tmax),meanCumuls, col="blue",pch=16,cex=0.25)  #stoch average
# end zeromeans ###################

# Sierra Leone data (WHO)
xpts<-thefile[weekIndexCase:200,3]-dayIndexCase #(index case)
points((thefile[weekIndexCase:200,9])~xpts,pch=20,cex=0.5,col="black")

legend(800, 30000, c("Simulated Cases (I)","Avg. simulated cases","standard error","Sierra Leone data"), 
       pch=c(16,16,10,16),col = c("red","blue","blue","black"),  adj = c(0, 0.6))

 dev.off()

#######################################
# Weekly Incidence
#######################################
 pdf("stoch_200_peak_noSTI.pdf", width = 5, height = 4)


xmaxweeks<-150
xspaces<-(135-weekIndexCase+1)/3  #135 is the week of the last date that prints (31 July 2016)
plot(1,xlim=c(1,xmaxweeks),ylim=c(0,1000),xaxt="n",xlab="Epidemic week",ylab="Weekly incidence",pch=20,main="Without STI")
axis(1,at=seq(1,xmaxweeks,by=xspaces),labels=FALSE,par(tcl=-0.5))
axis(1,c((1+9):(xmaxweeks+9)),thefile[(1+weekIndexCase-1):(xmaxweeks+weekIndexCase-1),11],par(tcl=0),lty=0,hadj=0.5)

display_runs<-200
avgweekly<-array(0,c(display_runs,715))
for(i in 1:display_runs)
{ weeklies<-matrix(c(results[i,,10],NA,NA,NA,NA,NA),nr=7,ncol=715) #matrix of new infecteds cut into weeks
  weeklies[which(is.na(weeklies))]<-0   # include zeros in the sum for each week
  yres<-colSums(weeklies,na.rm=TRUE)               # for each run, get the sum of new cases per week over the course of the epidemic
  avgweekly[i,]<-yres                   # store the weekly new case sum for each run
  xres<-c(1:715)                        # set up x axis by weeks
  lines(xres,yres,col="pink",lwd=0.1)   # plot each run as weekly new cases 
}
meanweekly<-colMeans(avgweekly, na.rm = TRUE)
lines(xres,meanweekly,col="darkorchid3",lwd=2)  #stoch average
SEMweekly<-apply(avgweekly,2,std.error, na.rm = TRUE)
meanplus<-meanweekly+SEMweekly
meanminus<-meanweekly-SEMweekly
lines(xres,meanplus,col="darkorchid1",lwd=0.75)  #stoch SEM
lines(xres,meanminus,col="darkorchid1",lwd=0.75)  #stoch SEM

# Sierra Leone data (WHO)
xpts<-(thefile[weekIndexCase:200,4])-weekIndexCase
points((thefile[weekIndexCase:200,8])~xpts,pch=20,cex=0.5,col="black")

legend(60, 1000, c("Each stochastic run","Mean of stochastic runs","Standard error","Sierra Leone data"), 
       pch=c(16,16,10,16),col = c("pink","darkorchid3","darkorchid1","black"),adj = c(0, 0.6))

 dev.off()

####################################
# tail incidence
####################################
 pdf("stoch_200_tail_noSTI.pdf", width = 6, height = 4)


startweek<-56 # 17-May-2015 is week 72 (sim week 56); 2-July-2017 is the last date that prints - week 183 (sim week 167)
endweek<-190
nseps<-(167-startweek)/4   # end
plot(1,xlim=c(startweek,endweek),ylim=c(0,15),xaxt="n",xlab="Epidemic week",ylab="Weekly incidence",pch=20,main="Without STI")
axis(1,at=seq(startweek,endweek,by=nseps),labels=FALSE,par(tcl=-0.5))
axis(1,c((startweek):(endweek)),thefile[(startweek+weekIndexCase-1):(endweek+weekIndexCase-1),11],par(tcl=0),lty=0)

display_runs<-200
avgweekly<-array(0,c(display_runs,715))
for(i in 1:display_runs)
{ weeklies<-matrix(c(results[i,,10],NA,NA,NA,NA,NA),nr=7,ncol=715) #matrix of new infecteds cut into weeks
  weeklies[which(is.na(weeklies))]<-0   # include zeros in the sum for each week
  yres<-colSums(weeklies,na.rm = TRUE)               # for each run, get the sum of new cases per week over the course of the epidemic
  avgweekly[i,]<-yres                   # store the weekly new case sum for each run
  xres<-c(1:715)                        # set up x axis by weeks
  lines(xres,yres,col="pink",lwd=0.1)   # plot each run as weekly new cases 
}

avgweekly[,50]     # find a good run with a tail distribution
lines(xres,avgweekly[8,],col="red",lwd=1.5) #one run
meanweekly<-colMeans(avgweekly, na.rm = TRUE)
lines(xres,meanweekly,col="darkorchid3",lwd=1)  #stoch average
SEMweekly<-apply(avgweekly,2,std.error, na.rm = TRUE)
meanplus<-meanweekly+SEMweekly
meanminus<-meanweekly-SEMweekly
lines(xres,meanplus,col="darkorchid1",lwd=0.5)  #stoch SEM
lines(xres,meanminus,col="darkorchid1",lwd=0.5)  #stoch SEM

xpts<-(thefile[,4])-weekIndexCase
points((thefile[,8])~xpts,pch=20,cex=0.5,col="black")

legend(120, 15, c("Example stochastic run"), pch=16,col = c("red"),adj = c(0, 0.6))

 dev.off()



