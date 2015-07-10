###########################################################################
#                                                                         #
#   EBOLA VIRUS DISEASE OUTBREAK IN NIGERIA:                              #
#   TRANSMISSION DYNAMICS AND RAPID CONTROL                               #
#                                                                         #
#   REFERENCE: ALTHAUS CL ET AL. EPIDEMICS. 2015;11:80-4.                 #
#   DOI: 10.1016/j.epidem.2015.03.001                                     #
#                                                                         #
#   Copyright (C) 2015 by Christian L. Althaus                            #
#   (christian.althaus@alumni.ethz.ch)                                    #
#                                                                         #
#   This program is free software; you can redistribute it and/or modify  #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 2 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
#   This program is distributed in the hope that it will be useful,       #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
#   You should have received a copy of the GNU General Public License     #
#   along with this program; if not, write to the                         #
#   Free Software Foundation, Inc.,                                       #
#   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             #
#                                                                         #
###########################################################################

# init
date()
rm(list = ls())

# Necessary libraries
library(deSolve)
library(chron)
library(bbmle)
library(mvtnorm)

# Initialize random number generator
set.seed(239847)

# Read the data
ebola <- read.csv("Ebola_outbreak_Nigeria_data.csv")
ebola$Date <- chron(as.character(ebola$Date),format=c(dates = "day mon year"))

# Definition of the intervention function
exponential <- function(t,beta0,beta1,k,tau1) {
  if(t < tau1) {
    	beta <- beta0
	} else {
		beta <- beta1 + (beta0-beta1)*exp(-k*(t-tau1))
	}
	return(beta)
}
Exponential <- Vectorize(exponential, vectorize.args = 't')

# Definition of the SEIR model
SEIR <- function(t, x, parms) {
	with(as.list(c(parms,x)),{
		beta <- exponential(t,beta0,beta1,k,tau1)
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
nll <- function(beta0,beta1,k,f,tau0,tau1,sigma,gamma) {
	pars <- c(beta0=beta0,beta1=beta1,k=k,f=f,tau0=tau0,tau1=tau1,sigma=sigma,gamma=gamma)
	pars <- trans(pars)
	times <- c(data$times+pars["tau0"],max(data$times+pars["tau0"])+1)
	simulation <- as.data.frame(ode(init,times,SEIR,parms=pars))
	ll <- sum(dpois(data$cases,diff(simulation$C),log=TRUE)) + sum(dpois(data$deaths,diff(simulation$D),log=TRUE))
	return(-ll)
}

# Parameter transformation
trans <- function(pars) {
	pars["beta0"] <- exp(pars["beta0"])
	pars["beta1"] <- exp(pars["beta1"])
	pars["k"] <- exp(pars["k"])
	pars["f"] <- plogis(pars["f"])
	pars["tau0"] <- exp(pars["tau0"])
	if(is.na(pars["tau1"])) {
		pars["tau1"] <- pars["tau0"]
	}
	return(pars)
}

# Prepare the data and set the initial values
data <- na.omit(ebola[c("Date","Cases","Deaths")])
names(data) <- c("times","cases","deaths")
begin <- chron("20 Jul 2014", format=c(dates = "day mon year")) 
data$times <- data$times - data$times[1]
N <- 1e6      
init <- c(S = N - 1, E = 0, I = 1, R = 0, D = 0, C = 0)

# Fit the model to the data
fixed <- c(sigma = 1/9.312799, gamma = 1/7.411374, tau0 = -Inf, tau1 = 3, beta1 = -Inf)
free <- c(beta0 = log(1), k = log(1), f = 0)
fit <- mle2(nll,start=as.list(free),fixed=as.list(fixed),method="Nelder-Mead",control=list(maxit=1e3))
fit

# Parameter sampling (for CIs and prediction interval) 
n_sim <- 1e4
m <- coef(fit,exclude.fixed=TRUE)
sigma <- vcov(fit)
sim_coef <- data.frame(rmvnorm(n_sim, mean = m, sigma = sigma))

# Timepoints
pars <- trans(coef(fit))
end <- chron("14 Sep 2014", format=c(dates = "day mon year"))
months <- chron(c("20 Jul 2014","20 Aug 2014","20 Sep 2014"), format=c(dates = "day mon year"))
weeks <- seq.dates(from = chron("20 Jul 2014", format = c(dates = "day mon year")), to = chron("14 Sep 2014", format = c(dates = "day mon year")), by = 14)
timepoints <- seq(0,end-begin,1)

# Plot the daily incidence and cumulative number of new cases and deaths
par(mfrow=c(2,2))

# Simulations
mod_sim_C_Inc <- numeric()
mod_sim_D_Inc <- numeric()
obs_sim_C_Inc <- numeric()
obs_sim_D_Inc <- numeric()
mod_sim_C_Cum <- numeric()
mod_sim_D_Cum <- numeric()
obs_sim_C_Cum <- numeric()
obs_sim_D_Cum <- numeric()
for(i in 1:n_sim) {
	pars <- trans(c(unlist(sim_coef[i,]),fixed))
	simulation <- as.data.frame(ode(init,timepoints,SEIR,parms=pars))
	mod_sim_C_Inc <- cbind(mod_sim_C_Inc, c(0,diff(simulation$C)))
	mod_sim_D_Inc <- cbind(mod_sim_D_Inc, c(0,diff(simulation$D)))
	mod_sim_C_Cum <- cbind(mod_sim_C_Cum, simulation$C)
	mod_sim_D_Cum <- cbind(mod_sim_D_Cum, simulation$D)
	obs_sim_dC <- rpois(length(timepoints) - 1, lambda = diff(simulation$C))
	obs_sim_dD <- rpois(length(timepoints) - 1, lambda = diff(simulation$D))
	obs_sim_C_Inc <- cbind(obs_sim_C_Inc, c(0, obs_sim_dC))
	obs_sim_D_Inc <- cbind(obs_sim_D_Inc, c(0, obs_sim_dD))
	obs_sim_C_Cum <- cbind(obs_sim_C_Cum, cumsum(c(0, obs_sim_dC)))
	obs_sim_D_Cum <- cbind(obs_sim_D_Cum, cumsum(c(0, obs_sim_dD)))
}

# Calculate the point-wise quantiles to construct an approximate 95% CI and a 95% prediction interval (cases).
mod_sim <- mod_sim_C_Inc
obs_sim <- obs_sim_C_Inc
x_grid <- timepoints
mod_q <- apply(mod_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
obs_q <- apply(obs_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
simulation <- as.data.frame(ode(init,x_grid,SEIR,parms=trans(coef(fit))))
plot(x_grid, mod_q[2,], ylim = c(0,4), ty="l", lty=2, lwd = 1, col = "red", xlab=NA, ylab="Daily incidence of new cases",frame=FALSE,axes=FALSE)
axis(1, weeks - begin, weeks)
axis(2)
polygon(x = c(x_grid, rev(x_grid)), y = c(obs_q[1,], rev(obs_q[2,])), col = rgb(1, 0, 0, alpha=0.2), border = NA)
lines(x_grid, mod_q[1,], lty = 2, lwd = 1, col = "red")
lines(x_grid, c(0,diff(simulation$C)), lwd = 1, col = "red")
points(data$times[1:57],data$cases[1:57],col="red")

mod_sim <- mod_sim_D_Inc
obs_sim <- obs_sim_D_Inc
mod_q <- apply(mod_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
obs_q <- apply(obs_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
plot(x_grid, mod_q[2,], ylim = c(0,2.5), ty="l", lty=2, lwd = 1, xlab=NA, ylab="Daily incidence of deaths",frame=FALSE,axes=FALSE)
axis(1, weeks - begin, weeks)
axis(2)
polygon(x = c(x_grid, rev(x_grid)), y = c(obs_q[1,], rev(obs_q[2,])), col = grey(0.7, alpha=0.5), border = NA)
lines(x_grid, mod_q[1,], lty = 2, lwd = 1)
lines(x_grid, c(0,diff(simulation$D)), lwd = 1)
points(data$times[1:57],data$deaths[1:57],pch=0)

mod_sim <- mod_sim_C_Cum
obs_sim <- obs_sim_C_Cum
x_grid <- timepoints
mod_q <- apply(mod_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
obs_q <- apply(obs_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
simulation <- as.data.frame(ode(init,x_grid,SEIR,parms=trans(coef(fit))))
plot(x_grid, 1+mod_q[2,], ylim = c(0,40), ty="l", lty=2, lwd = 1, col = "red", xlab=NA, ylab="Cumulative number of cases",frame=FALSE,axes=FALSE)
axis(1, weeks - begin, weeks)
axis(2)
polygon(x = c(x_grid, rev(x_grid)), y = c(1+obs_q[1,], rev(1+obs_q[2,])), col = rgb(1, 0, 0, alpha=0.2), border = NA)
lines(x_grid, 1+mod_q[1,], lty = 2, lwd = 1, col = "red")
lines(x_grid, 1+simulation$C, lwd = 1, col = "red")
points(data$times[1:57],1+cumsum(data$cases)[1:57],col="red")

mod_sim <- mod_sim_D_Cum
obs_sim <- obs_sim_D_Cum
mod_q <- apply(mod_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
obs_q <- apply(obs_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))
plot(x_grid, mod_q[2,], ylim = c(0,20), ty="l", lty=2, lwd = 1, xlab=NA, ylab="Cumulative number of deaths",frame=FALSE,axes=FALSE)
axis(1, weeks - begin, weeks)
axis(2)
polygon(x = c(x_grid, rev(x_grid)), y = c(obs_q[1,], rev(obs_q[2,])), col = grey(0.7, alpha=0.5), border = NA)
lines(x_grid, mod_q[1,], lty = 2, lwd = 1)
lines(x_grid, simulation$D, lwd = 1)
points(data$times[1:57],cumsum(data$deaths)[1:57],pch=0)

# Calculate the effective reproduction number vs time (including approximate pointwise CI)
beta <- function(t, pars) Exponential(t, pars[["beta0"]],pars[["beta1"]],pars[["k"]],pars[["tau1"]])

pars <- trans(coef(fit))
Re <- beta(timepoints, pars) / pars[["gamma"]]    
Re_sim <- numeric()
for(i in 1:n_sim) {
  pars_i <- trans(c(unlist(sim_coef[i,]), fixed))
  Re_i <- beta(timepoints, pars_i) / pars_i[["gamma"]]
  Re_sim <- cbind(Re_sim, Re_i)
}
Re_q <- apply(Re_sim, MAR = 1, FUN = quantile, probs = c(0.025, 0.975))

# Plot the effective reproduction number
par(mfrow=c(1,1))
control <- min(which(Re<1)) - 1.5
control_l <- min(which(Re_q[1,]<1)) - 1.5
control_u <- min(which(Re_q[2,]<1)) - 1.5
plot(timepoints,Re,type="n",col="blue",ylim=c(0,16),xlab=NA,ylab=bquote("Net reproduction number " ~ italic("R"["t"])),axes=FALSE,frame=FALSE)
polygon(x = c(timepoints, rev(timepoints)), y = c(Re_q[1,], rev(Re_q[2,])), col = rgb(0, 0, 1, alpha=0.2), border = NA)
lines(timepoints, Re, col = "blue")
lines(timepoints, Re_q[1,], col = 'blue', lty = 2)
lines(timepoints, Re_q[2,], col = 'blue', lty = 2)
axis(1, weeks - begin, weeks)
axis(2, 0:16)
lines(c(0,control),c(1,1))
lines(c(control,control),c(0,1))
points(control,1,pch=19)

# Plot outbreak probability as a function of R0
x_grid <- seq(1, 16, 0.01)
plot(x_grid, 1 - 1/x_grid, type = 'n', ylim = c(0, 1), xlab=bquote("Basic reproduction number " ~ italic("R")["0"]),ylab="Probability of outbreak from a single case",frame=FALSE,axes=FALSE)
axis(1,seq(1,16,3),seq(1,16,3))
axis(2)
countries <- c("Guinea","Sierra\nLeone","Liberia","Nigeria")
R0 <- c(1.51,2.53,1.59,Re[1])
polygon(x = c(1,Re_q[1, 1],Re_q[1, 1],Re_q[2, 1],Re_q[2, 1],1), y = c(1 - 1/Re_q[1, 1],1 - 1/Re_q[1, 1],0,0,1 - 1/Re_q[2, 1],1 - 1/Re_q[2, 1]), col = grey(0.7, alpha=0.5), border = NA)
lines(x_grid, 1 - 1/x_grid)
lines(c(Re[1],Re[1]),c(0,0.79))
lines(c(1,Re[1]),c(1-1/Re[1],1-1/Re[1]))
lines(c(Re_q[1, 1],Re_q[1, 1]),c(0,1-1/Re_q[1, 1]), lty = 2)
lines(c(Re_q[2, 1],Re_q[2, 1]),c(0,1-1/Re_q[2, 1]), lty = 2)
lines(c(1,Re_q[1, 1]),c(1-1/Re_q[1, 1],1-1/Re_q[1, 1]), lty = 2)
lines(c(1,Re_q[2, 1]),c(1-1/Re_q[2, 1],1-1/Re_q[2, 1]), lty = 2)
points(R0,1-1/R0,pch=19,col="red")
text(R0,1-1/R0,countries,pos=c(4,2,4,1),cex=0.7)
