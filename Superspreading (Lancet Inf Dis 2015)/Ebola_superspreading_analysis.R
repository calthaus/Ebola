###########################################################################
#                                                                         #
#   EBOLA SUPERSPREADING                                                  #
#                                                                         #
#   REFERENCE: ALTHAUS CL. LANCET INFECT DIS. 2015;15(5):507-8.           #
#   DOI: 10.1016/S1473-3099(15)70135-0                                    #
#                                                                         #
#   Copyright (C) 2015 by Christian L. Althaus                            #
#   (christian.althaus@alumni.ethz.ch)                                    #
#                                                                         #
#   This script is free software; you can redistribute it and/or modify   #
#   it under the terms of the GNU General Public License as published by  #
#   the Free Software Foundation; either version 2 of the License, or     #
#   (at your option) any later version.                                   #
#                                                                         #
#   This script is distributed in the hope that it will be useful,        #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#   GNU General Public License for more details.                          #
#                                                                         #
###########################################################################

# Necessary libraries
library(fitdistrplus)

# Number of individuals in the trees
n <- 152
# Number of secondary cases for all individuals
c1 <- c(1,2,2,5,14,1,4,4,1,3,3,8,2,1,1,4,9,9,1,1,17,
        2,1,1,1,4,3,3,4,2,5,1,2,2,1,9,1,3,1,2,1,1,2)
c0 <- c(c1,rep(0,n-length(c1)))
# Fitting a negative binomial distribution to the number of secondary cases
fit.cases <- fitdist(c0,"nbinom")
summary(fit.cases)
plot(fit.cases)

# Range of reported serial intervals
days <- 0:43
# Observed intervals for each day
frequency <- c(0,1,3,1,4,1,6,1,2,2,11,6,0,1,10,3,5,8,4,3,3,1,
               0,2,0,2,0,3,1,1,1,0,0,0,0,0,2,0,1,0,1,1,0,1)
d <- rep(days,frequency)
# Fitting a gamma distribution to the serial interval
fit.serial <- fitdist(d,"gamma")
summary(fit.serial)
plot(fit.serial)

# Set seed for random number generator
set.seed(645)
# Number of simulation runs
runs <- 1e2
# Number of initial cases
seed <- 1
# Initialize plot
plot(NA,xlim=c(0,100),ylim=c(0,100),xlab="Time (days)",
     ylab="Cumulative number of EVD cases",frame=FALSE)
# Set color scheme for different trajectories
cols <- sample(terrain.colors(runs))
# Simulate outbreak trajectories
for(i in 1:runs) {
    cases <- seed
	t <- rep(0,seed)
	times <- t
	while(cases > 0) {
		secondary <- rnbinom(cases,size=fit.cases$estimate[1],mu=fit.cases$estimate[2])
        t.new <- numeric()
		for(j in 1:length(secondary)) {
			t.new <- c(t.new,t[j] + rgamma(secondary[j],shape=fit.serial$estimate[1],
                       rate=fit.serial$estimate[2]))
		}
		cases <- length(t.new)
		t <- t.new
		times <- c(times,t.new)
	}
	lines(sort(times),1:length(times),col=cols[i],lwd=1)
	points(max(times),length(times),col=cols[i],pch=16)
}
