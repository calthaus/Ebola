###################################################################
#   Summary statistics for Monte Carlo simulations                #
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

library(car)
library(MASS)
library(plotrix)
library(userfriendlyscience)

###################################################################
# load data from stochastic simulation files

#    without STI   
out    <- read.table(file="Stoch_noSTI_summaries_1000.csv",header=T)

#    with STI   alpha = 3 months
outsti    <- read.table(file="Stoch_STI_3mo_summaries_1000.csv",header=T)

#    with STI   alpha = 6 months
outsti6    <- read.table(file="Stoch_STI_6mo_summaries_1000.csv",header=T)
###################################################################

#####################################################
#                                                   #
#    Tail Histograms for Stochastic Runs            #
#                                                   #
#####################################################

###################################################################

#    Generates results Figs. 4G, 4H, 4I

#    Sierra Leone Ebola incidence data based on raw case counts from:
#    WHO Ebola data and statistics: Data on new cases per epi week for Sierra Leone
#    18 Nov 2015 http://apps.who.int/gho/data/node.ebola-sitrep.ebola-country-SLE-20151118

thefile 	<- read.csv(file="data_WHO_SierraLeone_updated.csv",header=T)

weekIndexCase 	<-17  # April 23, 2014; week 17 ends April 27, 2014

##################################
##  For No STI
##################################

#week of last case
WeekLastCase<-as.integer(as.character(out[,5]))/7
hist(WeekLastCase,breaks=20,xlim=c(1,300))

pdf("weeklastcase_noSTI_hist.pdf", width = 4, height = 4)

lastweek<-375
xlabels <-thefile[weekIndexCase:(lastweek+weekIndexCase-1),11] #11= the date, 10 is the week number ; epiweek 17 is simulation week 1
xlength <- c(1:lastweek)
maxfreq <-c(0,150)
WeekLastCaseNoZero<-WeekLastCase
WeekLastCaseNoZero[which(out[,3]<50)]<-NA   # remove early die-outs
hist(WeekLastCaseNoZero,breaks=20,xlim=c(1,lastweek),ylim=maxfreq,xaxt="n",xlab=NULL,las=1,main="Week of Last Acute Case",col="grey") # xaxt = xaxis text
axis(1,at=seq(1,lastweek,by=25),labels=FALSE,par(tcl=-0.5))
axis(1,xlength,labels=FALSE,lty=0)
xl<-array(xlabels,lastweek)
xlabsby<-seq(1,lastweek,50)
text(x=xlabsby,par()$usr[3]-8,labels=xl[xlabsby],srt=60,adj=1,xpd=TRUE)

dev.off()

##################################
##  For STI 3 months Conv. Period
##################################

#week of last case
WeekLastCase<-as.integer(as.character(outsti[,10]))/7
hist(WeekLastCase,breaks=20,xlim=c(1,500))

pdf("weeklastcase_STI_3mo_hist.pdf", width = 4, height = 4)

lastweek<-375
xlabels <-thefile[weekIndexCase:(lastweek+weekIndexCase-1),11] #11= the date, 10 is the week number ; epiweek 17 is simulation week 1
xlength <- c(1:lastweek)
maxfreq <-c(0,150)
WeekLastCaseNoZero<-WeekLastCase
WeekLastCaseNoZero[which(outsti[,7]<50)]<-NA   # remove early die-outs
hist(WeekLastCaseNoZero,breaks=20,xlim=c(1,lastweek),ylim=maxfreq,xaxt="n",xlab=NULL,las=1,main="Week of Last Acute Case",col="grey") # xaxt = xaxis text
axis(1,at=seq(1,lastweek,by=25),labels=FALSE,par(tcl=-0.5))
axis(1,xlength,labels=FALSE,lty=0)
xl<-array(xlabels,lastweek)
xlabsby<-seq(1,lastweek,50)
text(x=xlabsby,par()$usr[3]-8,labels=xl[xlabsby],srt=60,adj=1,xpd=TRUE)

dev.off()

##################################
##  For STI 6 months Conv. Period
##################################

#week of last case
WeekLastCase<-as.integer(as.character(outsti6[,10]))/7
hist(WeekLastCase,breaks=20,xlim=c(1,500))

pdf("weeklastcase_STI_6mo_hist.pdf", width = 4, height = 4)

lastweek<-375
xlabels <-thefile[weekIndexCase:(lastweek+weekIndexCase-1),11] #11= the date, 10 is the week number ; epiweek 17 is simulation week 1
xlength <- c(1:lastweek)
maxfreq <-c(0,150)
WeekLastCaseNoZero<-WeekLastCase
WeekLastCaseNoZero[which(outsti6[,7]<50)]<-NA   # remove early die-outs
hist(WeekLastCaseNoZero,breaks=30,xlim=c(1,lastweek),ylim=maxfreq,xaxt="n",xlab=NULL,las=1,main="Week of Last Acute Case",col="grey") # xaxt = xaxis text
axis(1,at=seq(1,lastweek,by=25),labels=FALSE,par(tcl=-0.5))
axis(1,xlength,labels=FALSE,lty=0)
xl<-array(xlabels,lastweek)
xlabsby<-seq(1,lastweek,50)
text(x=xlabsby,par()$usr[3]-8,labels=xl[xlabsby],srt=60,adj=1,xpd=TRUE)

dev.off()

###################################################################

#####################################################
#                                                   #
#    Summary Statistics for Stochastic Runs         #
#                                                   #
#####################################################

#####################################################
#    Cumulative (total) Number of Cases
#####################################################

tot.out<-out["TotalNumCases"]
hist(tot.out[,1])

tot.outsti<-outsti["stiTotalNumCases"]
hist(tot.outsti[,1])

tot.outsti6<-outsti6["stiTotalNumCases"]
hist(tot.outsti6[,1])

means_out<-colMeans(tot.out)
SEM.out<-apply(tot.out,2,std.error)

means.outsti<-colMeans(tot.outsti,na.rm=T)
SEM.outsti<-apply(tot.outsti,2,std.error)

means.outsti6<-colMeans(tot.outsti6,na.rm=T)
SEM.outsti6<-apply(tot.outsti6,2,std.error)


xx<-sapply(c(tot.out),as.numeric)
yy<-sapply(c(tot.outsti),as.numeric)
zz<-sapply(c(tot.outsti6),as.numeric)
wilcox.test(xx,yy,conf.int=TRUE)
wilcox.test(xx,zz,conf.int=TRUE)

#####################################################
#    Peak Size (Incidence) Number of New Cases per day
#####################################################

ps.out<-out["EpidemicPeakNew"]
hist(ps.out[,1])

ps.outsti<-outsti["stiEpidemicPeakNew"] 
hist(ps.outsti[,1])

ps.outsti6<-outsti6["stiEpidemicPeakNew"]
hist(ps.outsti6[,1])

means_out<-colMeans(ps.out)
SEM.out<-apply(ps.out,2,std.error)

means.outsti<-colMeans(ps.outsti,na.rm=T)
SEM.outsti<-apply(ps.outsti,2,std.error)

means.outsti6<-colMeans(ps.outsti6,na.rm=T)
SEM.outsti6<-apply(ps.outsti6,2,std.error)

xxx<-sapply(c(ps.out),as.numeric)
yyy<-sapply(c(ps.outsti),as.numeric)
zzz<-sapply(c(ps.outsti6),as.numeric)
wilcox.test(xxx,yyy,conf.int=TRUE)
wilcox.test(xxx,zzz,conf.int=TRUE)


###################################################################
#    Subset data to only those runs that reached 50 or more total cumulative cases
###################################################################

out<-subset(out,out["TotalNumCases"]>49)

outsti<-subset(outsti,outsti["stiTotalNumCases"]>49)

outsti6<-subset(outsti6,outsti6["stiTotalNumCases"]>49)
###################################################################


#####################################################
#    Date of last Case
#####################################################

lastcase.out<-out["DayLastCase"]
hist(lastcase.out[,1])

lastcase.outsti<-outsti["stiDayLastCase"] 
hist(lastcase.outsti[,1])

lastcase.outsti6<-outsti6["stiDayLastCase"] 
hist(lastcase.outsti6[,1])

means_out<-colMeans(lastcase.out,na.rm=T)
SEM.out<-apply(lastcase.out,2,std.error,na.rm=T)

means.outsti<-colMeans(lastcase.outsti,na.rm=T)
SEM.outsti<-apply(lastcase.outsti,2,std.error,na.rm=T)
sd.outsti<-apply(lastcase.outsti,2,sd,na.rm=T)

means.outsti6<-colMeans(lastcase.outsti6,na.rm=T)
SEM.outsti6<-apply(lastcase.outsti6,2,std.error,na.rm=T)
sd.outsti6<-apply(lastcase.outsti6,2,sd,na.rm=T)

t.test(lastcase.out,lastcase.outsti)
t.test(lastcase.out,lastcase.outsti6)

xxxx<-sapply(c(lastcase.out),as.numeric)
yyyy<-sapply(c(lastcase.outsti),as.numeric)
zzzz<-sapply(c(lastcase.outsti6),as.numeric)
wilcox.test(xxxx,yyyy,conf.int=TRUE)
wilcox.test(xxxx,zzzz,conf.int=TRUE)

sum(lastcase.outsti>730,na.rm=T) # 2 years
sum(lastcase.outsti>0,na.rm=T)

sum(lastcase.outsti6>730,na.rm=T) # 2 years
sum(lastcase.outsti6>0,na.rm=T)


#####################################################
#    Date of Epidemic peak Incidence (max # new cases per day)
#####################################################

pt.out<-out["DayEpidemicPeakNew"]
hist(pt.out[,1])

pt.outsti<-outsti["stiDayEpidemicPeakNew"]
hist(pt.outsti[,1])

pt.outsti6<-outsti6["stiDayEpidemicPeakNew"]
hist(pt.outsti6[,1])

means_out<-colMeans(pt.out,na.rm=T)
SEM.out<-apply(pt.out,2,std.error,na.rm=T)

means.outsti<-colMeans(pt.outsti,na.rm=T)
SEM.outsti<-apply(pt.outsti,2,std.error,na.rm=T)
sd.outsti<-apply(pt.outsti,2,sd,na.rm=T)

means.outsti6<-colMeans(pt.outsti6,na.rm=T)
SEM.outsti6<-apply(pt.outsti6,2,std.error,na.rm=T)
sd.outsti6<-apply(pt.outsti6,2,sd,na.rm=T)

t.test(pt.out,pt.outsti)

xxxxx<-sapply(c(pt.out),as.numeric)
yyyyy<-sapply(c(pt.outsti),as.numeric)
zzzzz<-sapply(c(pt.outsti6),as.numeric)
wilcox.test(xxxxx,yyyyy)
wilcox.test(xxxxx,zzzzz)
###################################################################



###################################################################

#####################################################
#                                                   #
#    Comparing 2-fold reductions in eta             #
#    vs. 2-fold reductions in Conv Period (1/alpha) #
#                                                   #
#####################################################

###################################################################

# Generates results Suppl. Mat. Figure S4 and Table S3

###################################################################
#####################################################
#    Load data from Simulations: eta = 0.0005
#####################################################

# 3 month convalescent period
outsti0005    <- read.table(file="Stoch_STI_3mo_summaries_1000_eta0005.csv",header=T)

# 6 month convalescent period 
outsti6_0005  <- read.table(file="Stoch_STI_6mo_summaries_1000_eta0005.csv",header=T)


# subset these to also exclude those runs with fewer than 50 total cases

outsti0005<-subset(outsti0005,outsti0005["stiTotalNumCases"]>49)

outsti6_0005<-subset(outsti6_0005,outsti6_0005["stiTotalNumCases"]>49)


# dummy dataset for out (No STI) to compare with outsti (All STI Scenarios)

outdummy<-outsti[1:nrow(outdummy),]
outdummy$stiDayLastCase<-out$DayLastCase
outdummy$stiTotalNumCases<-out$TotalNumCases
outdummy$group<-out$group 


# set variable for comparisons "group"
outdummy$group<-"A-noSTI"
outsti0005$group<-"B-STI-3mo-0005"
outsti$group<-"C-STI-3mo-001"
outsti6_0005$group<-"D-STI-6mo-0005"
outsti6$group<-"E-STI-6mo-001"

#####################################################
#    Date of Last Case
#####################################################

# compare differences in means, correcting for unequal sample sizes and variances
bp<-rbind(outdummy,outsti0005)
bp<-rbind(bp,outsti)
bp<-rbind(bp,outsti6_0005)
bp<-rbind(bp,outsti6)
boxplot(bp$stiDayLastCase~bp$group, data=bp)

bp$group <- as.factor(bp$group)
posthocTGH(y=bp$stiDayLastCase, x=bp$group, method="games-howell")
oneway(y=bp$stiDayLastCase, x=bp$group, posthoc="games-howell")

(631- 606)/ 606 #  4.13  25 days STI 3mo: 0.001 vs. 0.0005 eta
(1088-935)/1088 # 14.0% 153 days STI 6mo: 0.001 vs. 0.0005 eta
(1088-631)/1088 # 42.0% 457 days STI 6mo vs. 3mo (eta=0.001)

42/14 # Ratio of effects of reducing ALPHA : reducing ETA = 3.0

# compare differences between variances
test1<-subset(bp,bp$group=="C-STI-3mo-001" | bp$group=="E-STI-6mo-001")
leveneTest(test1$stiDayLastCase~test1$group,data=test1)  #for hov

test2<-subset(bp,bp$group=="D-STI-6mo-0005" | bp$group=="E-STI-6mo-001")
leveneTest(test2$stiDayLastCase~test2$group, data=test2)  #for hov

# violin & boxplot for Fig. S4
g<-ggplot(bp, aes(x=rev(bp$group), y=bp$stiDayLastCase))
vbplot<-g+geom_violin(alpha=0.5, color="grey")+geom_boxplot(fill=c("gray","cadetblue1","cyan","tomato","red"),alpha=0.6)
vbplot<-g+geom_violin(alpha=0.5, color="grey")+geom_boxplot(fill=c("darkorchid","orchid3","red","tomato","gray"),alpha=0.6)+coord_flip()
vbplot + theme(panel.background = element_rect(fill = 'white', colour = 'white'))

