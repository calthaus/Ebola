###################################################################
#   Sensitivity analysis for deterministic models                 #
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

rm(list=ls(all=TRUE))

#required libraries
library(deSolve)
library(pse)

####################################################################
#  A. ODE deterministic model  & Sensitivity Analysis Code         #
#  B. Sensitivity analysis from saved files                        #
####################################################################

####################################################################
#
#           A. ODE & Sensitivity Analysis (de novo)
#
####################################################################

####################################################################
#  Generate latin hypercube                                        #
#  package 'pse', function 'LHS'                                   #
#                                                                  #
#  Input parameters to vary:                                       #
#   - R0                                                           #
#   - Sexual transmission rate (eta)                               #  
#   - Sexually-active males (proportion of pop,p)                  # 
#   - Recovery rate of Convalescent cases (alpha)                  # 
####################################################################
set.seed(451)

#latin hypercube of input parameters
factors <- c("R0","eta", "p","alpha","q")
q <- c("qunif", "qunif", "qunif","qunif","qunif")
q.arg <- list(list(min=1.26,max=2.53),          # R0     (from 1.26 (Alizon et al. 2015) to 2.53 (Althaus 2014))
              list(min=0.69897,max=1.30103),    # eta    (0.001, from log10(0.001*10000/2) - log10(0.001*2*10000))
              list(min=1.238666, max=1.840726), # p      (0.474 male*0.731 ages 15-44, from log10(0.474*0.731*100/2) to log10(0.474*0.731*100*2)) 
              list(min=1.640233, max=2.242293), # 1/alpha(83.84 days, from log10(83.84/2) to log10(83.84*2))
              list(min=1.133539, max=1.735599)) # q      (8.27 coital acts per month, from log10(8.27*12*2*100/365) to log10(8.27*12*100/365*2))

hypercube.size<-1000   # N gives number of parameter value combinations to try
nboot<-50              #nboot gives confidence interval for sensitivity measure

uncoupledLHS <- LHS(model=NULL, factors, N=hypercube.size, q, q.arg,nboot=nboot)        

#save hypercube:
write.table(get.data(uncoupledLHS), file="pse_cube_pqre.csv",col.names = F, row.names = F)

#get input value combinations (R0, eta, p, alpha) from hypercube
lh<-read.table("pse_cube_pqre.csv")


####################################################################
#   Model Output:                                                  #
#        1) Cumulative number of symptomatic cases (I)             #
#        2) End of the epidemic (first day that E+I(+C) = 0)       #
#        3) Number of symptomatic cases at peak of the epidemic    #
#            (max #I at any one time)                              #
#        4) Date of epidemic peak                                  #
#        5) Date last active case was seen                         #
####################################################################

#set number of iterations (number of parameter combinations to test)
imax<-nrow(lh)

# set output array
outs<-array(0,c(imax,16))
colnames(outs)<-c("i","R0","eta", "p","alpha", 
                  "woTotalNumCases", "woEndOfEpidemic","woEpidemicPeak","woDayEpidemicPeak","woDayLastCase", 
                  "stiTotalNumCases", "stiEndOfEpidemic","stiEpidemicPeak","stiDayEpidemicPeak","stiDayLastCase","q")

####################################################################
#   Run Models                                                     #
#  (with and without sexual transmission)                          #
#  for each set of input values                                    #
####################################################################

for (i in 1:imax)
{   
  #set parameter values for each run
  R0 <- lh[i,1]
  eta <- (10^(lh[i,2]))/10000
  p <- (10^(lh[i,3]))/100
  alpha <- 1/(10^(lh[i,4]))
  q <- (10^(lh[i,5]))/100
  
  
 ##################################################################
 # EVD transmission model * WITH * convalescent transmission (SEICR)
 ##################################################################
 model.wt <- function(t, x, parms) {
   with(as.list(c(parms,x)),{
     
     # calculate transmission rate (beta) to account for control measure implementation 
     beta <- ifelse(t<tau1,beta0,(beta1 + (beta0-beta1)*exp(-k*(t-tau1)))) 
     
     N <- S + E + I + C + R                     ## Total population size
     dS <- - beta*S*I - (eta*p*q*S*C)/N         ## Susceptible
     dE <- beta*S*I + (eta*p*q*S*C)/N - sigma*E ## Exposed (incubating)
     dI <- sigma*E - gamma*I                    ## symptomatic infections
     dC <- (1-f)*gamma*I - alpha*C              ## convalescent infections
     dR <- alpha*C                              ## recovered, immune, non-infectious
     dD <- f*gamma*I                            ## dead
     CaseCount<-sigma*E                         ## running count of number of total cases
     dST <- (eta*p*q*C*S)/N                     ## cumm no of sex transmission cases
     der <- c(dS,dE,dI,dR,dD,dC,CaseCount,beta,dST)
     list(der)
   })
 }

 ##################################################################
 # Set initial conditions
 ##################################################################
 par(mfrow=c(2,2), lwd=2)
 N0<-6.316e6    # Total population size of Sierra Leone
 N <- N0      # Initial population size, including infected individual(s)
 tmax<-5000   
 init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 0,CaseCount = 1, beta<-0,dST=0)

 timepoints <- seq(1,tmax,1) #INDEX CASE(S) DISEASED ON DAY 1

 pars <- c(R0 = R0, 
           R1 = 0.67,                   # lowest Ro with control measures
           beta0 = R0/(3.9*N0),         # initial transmission rate (before control measures)
           beta1 = 0.5/(3.9*N0),        # lowest transmission rate (with control measures)
           sigma = 1/11.4,              # duration of incubation period (rate of moving E -> I)
           gamma = 1/3.9,               # duration of infectious period (rate of moving I -> D or C)
           f = 0.69,                    # case fatality rate (probability of dying vs. surviving)
           eta = eta,                   # sexual transmission probability
           p = p,                       # proportion of C who are infectious (sexually active men)
           q = q,                       # daily number of sexual contacts
           alpha= alpha,                # 1/alpha, the rate of convalescent recovery      
           k= 0.01103038,               # control measure transmission decay
           tau1= 51.00493200)           # day control measures begin
 

 ##################################################################
 # run models and store outputs
 ##################################################################

 #store initial parameter values
 outs[i,1]=i
 outs[i,2]=R0
 outs[i,3]=eta
 outs[i,4]=p
 outs[i,5]=alpha
 outs[i,16]=q
  
 
 ##################################################
 # * WITH * convalescent transmission
 ##################################################
 
  #calculate transmission rate at each timestep according to Althaus et al. 2014
  pars["beta0"] <- pars["R0"]*pars["gamma"]/N
  simulation <- as.data.frame(ode(init,timepoints,model.wt,parms=pars))
  
  #store outputs
  #              Cumulative Cases
  TotalNumCases<-round(simulation$CaseCount[tmax])  #cumulative number of symptomatic cases (I) at the end of epidemic
  outs[i,11]<-TotalNumCases
  
  #              End of Epidemic 
  simulation$CurrentCases<-round(simulation$E+simulation$I+(p*simulation$C))  #all active cases for each day (E+I+C)
  EndOfEpidemic<-simulation$time[min(which(simulation$CurrentCases==0))]  #first day there are no more active cases (end of the epidemic)
  outs[i,12]<-EndOfEpidemic
  
  #              Peak Epidemic Size
  SizeEpidemicPeak<-round(max(simulation$I))   #the peak of the epidemic (max #I at any one time) ; should make this a percentage for reporting    
  outs[i,13]<-SizeEpidemicPeak
  
  #              Day Epidemic Peaks
  DayEpidemicPeak<-simulation$time[which(simulation$I==max(simulation$I))] #day at which the epidemic peaks
  outs[i,14]<-DayEpidemicPeak
  
  #              Day of Last Symptomatic Case
  lastcase<-which(round(simulation$I)>0.5)
  DayLastCase<-simulation$time[max(lastcase)] #Last day a symptomatic case was seen
  outs[i,15]<-DayLastCase
  
 } 

print("end of loops for sensitivity analysis")

##################################################################
# Save stored outputs for sensitivity analysis                   #
##################################################################
rewrite<-F  # T or F
if(rewrite==T){write.table(outs,file="Det_sensitivity.csv")}

##################################################################
#  Sensitivity Analysis                                          #
#  package 'pse', function 'tell'                                #
#                                                                # 
# requires from above:                                           #
#   - the hypercube (uncoupledLHS)                               #
#   - results(outs)                                              #
##################################################################

# load outputs from either run (outs) or saved file (Det_sensitivity.csv)
load.from.saved.file<-F  #T or F
if(load.from.saved.file==T){sensitivity.output<-read.table(file="Det_sensitivity.csv",header=T)}
if(load.from.saved.file==F){sensitivity.output<-outs}

# output value to test =     Length of Epidemic (with sexual transmission)
output.to.test<-sensitivity.output[,"stiEndOfEpidemic"]

#sensitivity anlaysis
coupledLHS <- tell(uncoupledLHS, output.to.test,nboot=nboot)  #nboot - std dev for each parameter

#plots
plotprcc(coupledLHS, ylab="Sensitivity of End of Epidemic (E+I+pC=0)")

#stats
coupledLHS

# Day of Last Case
output.to.test<-sensitivity.output[,"stiDayLastCase"]
coupledLHS <- tell(uncoupledLHS, output.to.test,nboot=nboot)  #nboot - std dev for each parameter
plotprcc(coupledLHS,ylab="Sensitivity of Day of Last Case (with Convalescent Transmission)")
coupledLHS

# Cumulative Number of Cases
output.to.test<-sensitivity.output[,"stiTotalNumCases"]
coupledLHS <- tell(uncoupledLHS, output.to.test,nboot=nboot)  #nboot - std dev for each parameter
plotprcc(coupledLHS,ylab="Sensitivity of Total Size of the Epidemic (with Convalescent Transmission)")
coupledLHS

# Day of Epidemic Peak
output.to.test<-sensitivity.output[,"stiDayEpidemicPeak"]
coupledLHS <- tell(uncoupledLHS, output.to.test,nboot=nboot)  #nboot - std dev for each parameter
plotprcc(coupledLHS,ylab="Sensitivity of Day of Epidemic Peak (with Convalescent Transmission)")
coupledLHS

# Size of Epidemic Peak
output.to.test<-sensitivity.output[,"stiEpidemicPeak"] #size epidemic peak (#I)
coupledLHS <- tell(uncoupledLHS, output.to.test,nboot=nboot)  #nboot - std dev for each parameter
plotprcc(coupledLHS,ylab="Sensitivity of Size (current I) of Epidemic Peak (with Convalescent Transmission)")
coupledLHS

