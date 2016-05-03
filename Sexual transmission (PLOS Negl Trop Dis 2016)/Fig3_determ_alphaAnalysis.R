###################################################################
#                                                                 #
#  Deterministic model simulations - Alpha parameter analysis     #
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

#Load Libraries

library(deSolve)
library(pse)

###################################################################


#set initial conditions
N0<-6.316e6    
tmax <-2000   # day
tmax2 <-1400
init <- c(S = N0 - 1, E = 0, I = 1, R = 0, D = 0, C = 0,CaseCount = 1, ST=0)
initWO <- c(S = N0 - 1, E = 0, I = 1, R = 0, D = 0,CaseCount = 1)

timepoints2 <- 0:tmax
timepoints <- 0:tmax2

#set parameters
R0 = 2.133485 
R1 = 0.67                   # lowest Ro with control measures
beta0 = R0/(3.9*N0)         # initial transmission rate (before control measures)
beta1 = R1/(3.9*N0)         # lowest transmission rate (with control measures)
sigma = 1/11.4              # duration of incubation period (rate of moving E -> I)
gamma = 1/3.9               # duration of infectious period (rate of moving I -> D or C)
f = 0.69                    # case fatality rate (probability of dying vs. surviving)
p = 0.474*0.731             # proportion of C able to transmit (sexually-active men)
q = ((8.27*12)/365)         # daily rate of sexual contact (avg. sex acts per month/days in year)
eta = 0.001                 # per sex act rate of transmission
k= 0.01103038               # control measure transmission decay 
tau= 51.00493200            # day control measures are implemented

parsWO <- c(R0 = R0, R1 = R1,beta0 = beta0, beta1 = beta1, sigma = sigma,     
            gamma = gamma, f = f, k= k,tau= tau)


# EVD transmission model * WITHOUT * convalescent transmission 
#-------------------------------------------------------------
model.wot <- function(t, x, parms) {
  with(as.list(c(parms,x)),{ 
    
    # calculate transmission rate (beta) to account for control measure implementation 
    beta <- ifelse(t<tau,beta0,(beta1 + (beta0-beta1)*exp(-k*(t-tau))))  
    
    NN <- S + E + I  + R            ## Total population size
    dS <- - beta*S*I                ## Susceptible
    dE <- beta*S*I - sigma*E        ## Exposed (incubating)
    dI <- sigma*E - gamma*I         ## symptomatic infections
    dR <- (1-f)*gamma*I             ## recovered, immune, non-infectious
    dD <- f*gamma*I                 ## dead
    CaseCount<-sigma*E              ## running count of number of total cases
    der <- c(dS,dE,dI,dR,dD,CaseCount)
    list(der)
  })
}

#run model.wot
simulationWO <- as.data.frame(ode(initWO,timepoints2,model.wot,parms=parsWO))


# EVD transmission model * WITH * convalescent transmission
# ---------------------------------------------------------
model.wt <- function(t, x, parms) {
  with(as.list(c(parms,x)),{
    
    # calculate transmission rate (beta) to account for control measure implementation 
    beta <- ifelse(t<tau,beta0,(beta1 + (beta0-beta1)*exp(-k*(t-tau)))) 
    
    N <- S + E + I + C + R                     ## Total population size
    dS <- - beta*S*I - (eta*p*q*C*S)/N         ## Susceptible
    dE <- beta*S*I + (eta*p*q*C*S)/N - sigma*E ## Exposed (incubating)
    dI <- sigma*E - gamma*I                    ## symptomatic infections
    dC <- (1-f)*gamma*I - alpha*C              ## convalescent infections
    dR <- alpha*C                              ## recovered, immune, non-infectious
    dD <- f*gamma*I                            ## dead
    CaseCount<-sigma*E                         ## running count of number of total cases
    dST <- (eta*p*q*C*S)/N                     ## cumm no of sex transmission cases
    der <- c(dS,dE,dI,dR,dD,dC,CaseCount,dST)
    list(der)
  })
}


alphaVals <- c(1/(83.84), 1/(175.09), 1/(266.34))   # range of alpha values (3,6,9 months)

simulation <- list()                                # output of loop

#loop on the number of runs
for(i in 1:length(alphaVals))
{
  pars <-  c(R0 = R0, R1 = R1,beta0 = beta0, beta1 = beta1, sigma = sigma,     
             gamma = gamma, f = f,  p = p, q = q, alpha= alphaVals[i], k= k,tau= tau,
             eta = eta) #list of parameters
  
  
  simulation <- as.data.frame(ode(init,timepoints2,model.wt,parms=pars))  #run model.wt
  
  var.names=c("C","I","E","R","D","ST","CaseCount")   #list of all variables of interest
  obj.names1=paste("sim",var.names,sep="")            #building label: combine sim to letters of var.names
  obj.names=paste(obj.names1,toString(i),sep="")      #building label: combine simLetter with the value of the iteration
  
  mat=subset(simulation,select=var.names)    #makes a matrix of all the variable of interest with their outputs in columns
  
  for (j in 1:length(var.names)){     # assigns new object names to each variable's run in order to save and plot later
    assign(obj.names[[j]],mat[[j]])
  }
}


# Plots       Generates Fig. 3
# -----------------------------
 pdf("deterministic alpha sensitivity.pdf")

par(mfrow=c(2,2),oma = c(0,0,2,0))

int2 <- 5  #for lty
int3 <- 4
int4 <- 3
wd <- 1.8  #for lwd
endt<- tmax2+1 #for timepoints = 100 use 1001

# get day of last case for each model
endsim1<-max(which(simI1[c(1:endt)]>0.5))
endsim2<-max(which(simI2[c(1:endt)]>0.5))
endsim3<-max(which(simI3[c(1:endt)]>0.5))
endsimWO<-max(which(simulationWO$I>0.5))

plot(timepoints, simC3[c(1:endt)], col="blue", lty=int4,ylim=c(0,2500), lwd=wd, type="l", xlab="Time (days)", ylab="Number of individuals",frame=FALSE)
lines(timepoints, simE3[c(1:endt)], col="orange", lty=int4,lwd=wd)
lines(timepoints, simI3[c(1:endt)], col="red", lty=int4,lwd=wd)
lines(timepoints, simC1[c(1:endt)], col="blue",lty=int2,lwd=wd)
lines(timepoints, simI1[c(1:endt)], col="red",lty=int2,lwd=wd)
lines(timepoints, simE1[c(1:endt)], col="orange", lty=int2,lwd=wd)
lines(timepoints, simC2[c(1:endt)], col="blue", lty=int3,lwd=wd)
lines(timepoints, simE2[c(1:endt)], col="orange", lty=int3,lwd=wd)
lines(timepoints, simI2[c(1:endt)], col="red", lty=int3,lwd=wd)
lines(timepoints, simulationWO$E[c(1:endt)], col="orange",lwd=wd)
lines(timepoints, simulationWO$I[c(1:endt)], col="red",lwd=wd)
lines(rep.int(timepoints[endsim1], 600),c(1:600),col="black",lty=int2,lwd=wd )
lines(rep.int(timepoints[endsim2], 600),c(1:600),col="black",lty=int3,lwd=wd )
lines(rep.int(timepoints[endsim3], 600),c(1:600),col="black",lty=int4,lwd=wd )
lines(rep.int(timepoints[endsimWO], 600),c(1:600),col="black",lwd=wd )
legend("topright",legend=c("Exposed","Infected","Convalescent","Day of last case"), lty=1, lwd=wd,col=c("orange","red","blue","black"),bty="n",cex=0.65)


plot(timepoints, simCaseCount1[c(1:endt)], lty=int2,ylim=c(0,14000), col="purple",type="l",  xlab="Time (days)", ylab="Number of individuals",frame=FALSE)
lines(timepoints, simR1[c(1:endt)], lty=int2, col="darkgreen",lwd=wd)
lines(timepoints, simD1[c(1:endt)], lty=int2, col="brown",lwd=wd)
lines(timepoints, simCaseCount2[c(1:endt)], col="purple", lty=int3,lwd=wd)
lines(timepoints, simR2[c(1:endt)], col="darkgreen", lty=int3,lwd=wd)
lines(timepoints, simD2[c(1:endt)], col="brown", lty=int3,lwd=wd)
lines(timepoints, simCaseCount3[c(1:endt)], col="purple", lty=int4,lwd=wd)
lines(timepoints, simR3[c(1:endt)], col="darkgreen", lty=int4,lwd=wd)
lines(timepoints, simD3[c(1:endt)], col="brown", lty=int4,lwd=wd)
lines(timepoints, simulationWO$CaseCount[c(1:endt)], col="purple",lwd=wd)
lines(timepoints, simulationWO$R[c(1:endt)], col="darkgreen",lwd=wd)
lines(timepoints, simulationWO$D[c(1:endt)], col="brown",lwd=wd)
legend("top",legend=c("Recovered","Convalescent","Cum. no. cases"), lty=1, horiz=TRUE, lwd=wd,col=c("darkgreen","brown","purple"),bty="n",cex=0.55)


plot(timepoints2, c(0,diff(simST3)), lty=int4,lwd=wd, ylim=c(0,0.3),type="l",xlab="Time (days)", ylab="Daily incidence of sexual transmission",frame=FALSE)
lines(timepoints2, c(0,diff(simST2)), lty=int3,lwd=wd)
lines(timepoints2, c(0,diff(simST1)), lty=int2,lwd=wd)
legend("topright",legend=c("3 months", "6 months", "9 months"),lwd=wd, lty=c(int2,int3,int4), title="1/alpha", cex=0.65)

plot(timepoints, simST3[c(1:endt)], lty=int4,lwd=wd, ylim=c(1e0,1e2),type="l",xlab="Time (days)", ylab="Cumulative number of sexual transmissions",frame=FALSE)
lines(timepoints, simST2[c(1:endt)], lty=int3,lwd=wd)
lines(timepoints, simST1[c(1:endt)], lty=int2,lwd=wd)

mtext("Sensitivity to convalescent recovery period", outer = TRUE, cex = 1) #title outside the margins

 dev.off()

