###########################################################################
#                                                                         #
#   RAPID DROP IN THE REPRODUCTION NUMBER DURING THE EBOLA OUTBREAK       #
#   IN THE DEMOCRATIC REPUBLIC OF CONGO                                   #
#                                                                         #
#   REFERENCE: ALTHAUS CL. PEERJ PREPRINTS. 2015;3:E1279.                 #
#   DOI: 10.7287/peerj.preprints.1041v1                                   #
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

# init
rm(list = ls())

# Necessary libraries
library(chron)
library(R0)

# Initialize random number generator
set.seed(679236)

# Read the data
ebola <- read.csv("Ebola_outbreak_DRC_data.csv")
ebola$Date <- chron(as.character(ebola$Date), format=c(dates = "day mon year"))
ebola$Date <- as.Date(ebola$Date)
ebola <- ebola[1:71,]

# Estimating case reproduction number
mGT <- generation.time("gamma", c(15.3,9.3)) # Generation time (serial interval) as reported by WHO Ebola Response Team (2014, NEJM)
TD <- est.R0.TD(ebola$Cases, mGT, t = ebola$Date, begin = 1, end = as.numeric(max(ebola$Date)-min(ebola$Date)), nsim = 1e4, correct = FALSE)

# Plot the case reproduction number during the outbreak
l <- length(ebola$Date)-1
plot(TD)
plot(NA,ty="n",frame=FALSE,axes=FALSE,xlim=c(as.Date("2014-07-26"),as.Date("2014-10-18")),ylim=c(0,8),xlab=NA,ylab=bquote("Case reproduction number " ~ italic("R")),cex.sub = 0.75)
polygon(x = c(ebola$Date[1:l], rev(ebola$Date[1:l])), y = c(TD$conf.int[1:l,1], rev(TD$conf.int[1:l,2])), col = rgb(0, 0, 1, alpha=0.2), border = NA)
lines(ebola$Date,c(TD$R,NA), col = "blue")
axis(1, seq(as.Date("2014-07-26"),as.Date("2014-10-18"),14), c("26 Jul 2014",NA,"23 Aug 2014",NA,"20 Sep 2014",NA,"18 Oct 2014"))
axis(2, 0:8)
control <- as.Date("2014-08-16")
lines(c(0,control),c(1,1))
lines(c(control,control),c(0,1))
points(control,1,pch=19)
