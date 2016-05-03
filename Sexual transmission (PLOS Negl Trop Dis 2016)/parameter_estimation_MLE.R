###################################################################
#                                                                 #
#  Maximum-likelihood estimation of model parameters              #
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

# init
rm(list = ls())

# Necessary libraries
library(stats4)
library(deSolve)
library(chron)
library(bbmle)

# Read the data
ebola <- read.csv("data_WHO_SierraLeone_mle.csv",skip=3)

# Select only incidence counts between these two time points
begin <- chron("25 May 2014", format=c(dates = "day mon year"))
end <- chron("13 Sep 2015", format=c(dates = "day mon year"))
timepoints <- seq(begin,end,7)
ebola <- as.data.frame(cbind(timepoints,ebola[21:89,2:5]))

# Add last two entries of situation report to patient data base
ebola[68:69,4:5] <- ebola[68:69,2:3]
ebola <- ebola[,-2:-3]
ebola[,2] <- ebola[,2]+ebola[,3]
ebola[,3] <- cumsum(ebola[,2])
names(ebola) <- c("Date","Cases","Cumulative")

# Definition of the SEIR model
SEIR <- function(t, x, parms) {
    with(as.list(c(parms,x)),{
		# Exponential reduction in transmission rate
        ifelse(t < tau1, beta <- beta0, beta <- beta1 + (beta0-beta1)*exp(-k*(t-tau1)))
        # Density-dependent transmission
		N <- S + E + I + R + D
		dS <- - beta/N*S*I
		dE <- beta/N*S*I - sigma*E
		dI <- sigma*E - gamma*I
		dR <- (1-f)*gamma*I
		dD <- f*gamma*I
		dC <- sigma*E
		der <- c(dS,dE,dI,dR,dD,dC)
		list(der)
	})
}

# Negative log-likelihood
nll <- function(beta0,beta1,k,f,tau0,tau1,sigma,gamma,disp) {
	pars <- c(beta0=beta0,beta1=beta1,k=k,f=f,tau0=tau0,tau1=tau1,sigma=sigma,gamma=gamma,disp=disp)
	pars <- trans(pars)
	times <- c(0,min(data$times)+pars["tau0"]-7,data$times+pars["tau0"])
	simulation <- as.data.frame(ode(init,times,SEIR,parms=pars))
	ll <- sum(dnbinom(data$cases,size=pars["disp"]*diff(simulation$C)[-1],mu=diff(simulation$C)[-1],log=TRUE))
	return(-ll)
}

# Parameter transformation
trans <- function(pars) {
	pars["beta0"] <- exp(pars["beta0"])
	pars["beta1"] <- exp(pars["beta1"])
	pars["k"] <- exp(pars["k"])
	pars["f"] <- plogis(pars["f"])
	pars["tau0"] <- exp(pars["tau0"])
	pars["tau1"] <- exp(pars["tau1"])
	pars["disp"] <- exp(pars["disp"])
	return(pars)
}

# Prepare the data and set the initial values
data <- na.omit(ebola[c("Date","Cases")])
names(data) <- c("times","cases")
data$times <- data$times - data$times[1]
N <- 6.316e6      
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 0)

# Fit the model to the data
fixed <- c(sigma = 1/11.4, gamma = 1/(15.3-11.4), f = qlogis(0.69), tau0 = log(32))
free <- c(beta0 = log(2.02/(15.3-11.4)), beta1 = log(0.5/(15.3-11.4)), tau1 = log(32), k = log(0.01), disp = log(0.18))
fit <- mle2(nll,start=as.list(free),fixed=as.list(fixed),method="Nelder-Mead",control=list(maxit=1e3))
summary(fit)
trans(coef(fit))

# Plot
begin <- chron("25 May 2014", format=c(dates = "day mon year")) 
end <- chron("13 Sep 2015", format=c(dates = "day mon year"))
months <- seq(chron("01 Jun 2014", format=c(dates = "day mon year")),chron("01 Oct 2015", format=c(dates = "day mon year")),by="months")
pars <- trans(coef(fit))
times <- c(0,min(data$times)+pars["tau0"]-7,data$times+pars["tau0"])
simulation <- as.data.frame(ode(init,times,SEIR,parms=pars))
plot(data$times, diff(simulation$C)[-1], ylim=c(0,700), ty="l", lwd = 1, col = "red",xlab=NA,ylab="Weekly incidence of new cases",frame=FALSE,axes=FALSE)
polygon(x = c(data$times, rev(data$times)), y = c(qnbinom(0.025,size=pars["disp"]*diff(simulation$C)[-1],mu=diff(simulation$C)[-1]), rev(qnbinom(0.975,size=pars["disp"]*diff(simulation$C)[-1],mu=diff(simulation$C)[-1]))), col = rgb(1, 0, 0, alpha=0.1), border = NA)
points(data$times,data$cases,col="red")
axis(1, months - begin, months)
axis(2)
