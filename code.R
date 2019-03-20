################################################################################
#AUTHORS
#
#
################################################################################


################################################################################
#SUMMARY (including basic steps of the analysis)
#

#SETUP
################################################################################
#working folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#functions
source("functions.R")

#libraries
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggdendro)
library(viridis)
library(nlme)
library(tidyr)
library(rioja)

#importing data
load("data.RData")


#0. SUBSETTING DATA
###############################################################################
#removing first 60 samples
bale.data <- bale.data[61:nrow(bale.data), ]
#removing ages older than 13.7k
bale.data <- bale.data[bale.data$age <= 13.779,]



#1. CLUSTERING CHARCOAL ACCUMULATION RATE
################################################################################
#subsetting bale.data and removing NAs
bale.charcoal <- na.omit(
  bale.data[, c("age", "charcoal.acc.rate")]
)

#log-transform to reduce skewness
charcoal.log <- log10(bale.charcoal$charcoal.acc.rate + 1)

#standardising for clustering
charcoal.scaled <- data.frame(scale(charcoal.log))

#Computing distance matrix for clustering with CHAR and BA only
charcoal.dist <- dist(charcoal.scaled)

#constrained hierarchical clustering
charcoal.cluster <- chclust(charcoal.dist, method="coniss")

#broken-stick analysis
charcoal.cluster.bstick <- bstick(charcoal.cluster)

#finding number of groups programmatically
charcoal.cluster.bstick$minimum <- charcoal.cluster.bstick$dispersion + charcoal.cluster.bstick$bstick
number.of.groups <- charcoal.cluster.bstick[which.min(charcoal.cluster.bstick$minimum), "nGroups"]

#Separate cluster zones
charcoal.cluster.zones <- cutree(charcoal.cluster, k=number.of.groups)

#Joining grouping vector with bale.charcoal
bale.charcoal$clustering.zone <- charcoal.cluster.zones


#2. CHRONOSTRATIGRAPHICAL DIAGRAM
###############################################################################
plotDiagram(pollen.data=bale.data,
            fire.data=bale.charcoal, 
            fire.cluster=charcoal.cluster, 
            dendrogram=TRUE, 
            title="Garba Guracha (3950 masl)", 
            title.size=16, 
            text.size=6,
            width=12, 
            height=6, 
            filename="Chronostratigraphical_diagram.pdf")

#removing objects we don't need any longer
rm(fire, fire.cluster, fire.log, fire.rotation.period, fire.scaled, fire.cluster.zones, fire.dist, min.age, max.age, zone, number.of.groups, fire.cluster.bstick)
gc()


#3. FEEDBACK BETWEEN FIRE AND ERICA
###############################################################################

#3.1 INTERPOLATING CHARCOAL INTO A REGULAR GRID
#--------------------------------------
#subsetting erica and charcoal and rounding ages
erica <- na.omit(bale.data[, c("age","ericaceae.par")])
erica$age <- round(erica$age, 3)
char <- na.omit(bale.data[, c("age", "charcoal.acc.rate")])
char$age <- round(char$age, 3)

#interpolating charcoal to a regular grid (10 years intervals)
char.interpolated<-interpolateDatasets(
  datasets.list=list(char=char), 
  age.column.name="age", 
  interpolation.time.step=0.01 #ka
)

#matching decimal positions in ages
char.interpolated$age<-round(char.interpolated$age, 2)
char$age <- round(char$age, 2)

#replacing interpolated values by observed ones where possible
#for every age in char.interpolated
for(i in char.interpolated$age){
  
  #getting observed value for the given age
  observed.i <- char[char$age == i,"charcoal.acc.rate"]
  
  #if observed.i is not empty
  if(length(observed.i) > 0){
    
    #replacing interpolated by observed
    char.interpolated[
      char.interpolated$age == i, 
      "charcoal.acc.rate"
      ] <- max(observed.i) #if two or more observed values with same age
    
  }
}

#plot to compare observed and interpolated char curves
plot.1 <- ggplot() + 
  geom_line(data=char, aes(x=age, y=charcoal.acc.rate), col="gray70", size=0.6) +
  geom_line(data=char.interpolated, aes(x=age, y=charcoal.acc.rate), color="red4", size=0.5, alpha=0.5) +
  ylab("CHAR (particles/cm2 yr)") +
  scale_x_continuous(breaks=seq(0, ceiling(max(bale.data$age)), by=1)) +
  xlab("Age (ka BP)")

plot.1.zoom <- ggplot() + 
  geom_line(data=char, aes(x=age, y=charcoal.acc.rate), col="gray70", size=0.6) +
  geom_line(data=char.interpolated, aes(x=age, y=charcoal.acc.rate), color="red4", size=0.5, alpha=0.5) +
  ylab("") +
  scale_x_continuous(breaks=seq(0, ceiling(max(bale.data$age)), by=1)) +
  xlab("Age (ka BP)") +
  coord_cartesian(xlim=c(12, 13.5), ylim=c(0, 20))

plot_grid(plot.1, plot.1.zoom, rel_widths = c(1, 0.6), labels = "auto")

#erica ages should be within char.interpolated ages
erica <- erica[erica$age >= min(char.interpolated$age), ]
erica <- erica[erica$age <= max(char.interpolated$age), ]


#3.2 CORRELATION BETWEEN ERICA AND CHAR
#--------------------------------------
#rounding erica ages
erica$age <- round(erica$age, 2)

#pairing samples of erica and interpolated charcoal
erica.char <- merge(erica, char.interpolated, by="age")

#modeling erica as a function of char
#using generalised least squares with gls "fits a linear model using generalized least squares. The errors are allowed to be correlated and/or have unequal variances."
erica.char.gls <- gls(ericaceae.par ~ charcoal.acc.rate, data=erica.char)
summary(erica.char.gls)

#pseudo R2 (gls doesn't provide R2)
cor(erica.char$ericaceae.par, predict(erica.char.gls))^2
# [1] 0.1644662

#adding predicted and residual values to erica.char
erica.char$ericaceae.par.predicted <- predict(erica.char.gls)

#plotting erica vs char
ggplot(data=erica.char, aes(x=charcoal.acc.rate, y=ericaceae.par)) + 
  geom_point(shape=21, fill="gray50", color="black", size=3, alpha=0.5) +
  geom_line(aes(x=charcoal.acc.rate, y=ericaceae.par.predicted), size=2, color="red4", alpha=0.6) +
  xlab("CHAR (particles/cm2 yr)") + 
  ylab("Erica PAR (pollen grains/cm2yr)") +
  ggtitle("Erica vs. CHAR") +
  theme(text=element_text(size=12), plot.title=element_text(size = 16))
ggsave(filename="output/erica_vs_char.pdf", width=6, height=4)



#3.3 GENERATING LAGGED DATA
#--------------------------------------
#generates lagged charcoal data after every erica sample, up to 25 lags (250 years)
lags<-1:101
lag.data.forward <- forwardLags(lags=lags, reference.data=erica, data.to.lag=char.interpolated)
lag.data.backward <- backwardLags(lags=lags, reference.data=erica, data.to.lag=char.interpolated)

#scatterplots of lagged datasets
plotLags(forward.lagged.data=lag.data.forward, backward.lagged.data=lag.data.backward, lags=lags, filename="output/lagged_data_scatterplots.pdf")


#3.4 ANALYSING LAGGED DATA
#--------------------------------------

#fitting a GLS model per lag
backward.results <- modelLagData(model.formula="ericaceae.par ~ charcoal.acc.rate", lagged.data=lag.data.backward)
forward.results <- modelLagData(model.formula="charcoal.acc.rate ~ ericaceae.par", lagged.data=lag.data.forward)

#computing GLS models with randomized responses (this can take a while depending on the number of iterations)
backward.results.random <- modelRandomLagData(lagged.data=lag.data.backward, model.formula="ericaceae.par ~ charcoal.acc.rate", iterations=1000)
forward.results.random <- modelRandomLagData(lagged.data=lag.data.forward, model.formula="charcoal.acc.rate ~ ericaceae.par", iterations=1000)

#plotting model results
plotModelOutput(forward.results=forward.results, 
              forward.results.random=forward.results.random, 
              backward.results=backward.results, 
              backward.results.random=backward.results.random, 
              filename="output/Figure_2.pdf", 
              width=9, 
              height=6, 
              title.size=16, 
              text.size=10)



