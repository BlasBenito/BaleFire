################################################################################
#AUTHORS
#Blas M. Benito.
#Graciela Gil-Romera
#Carole Adolf
################################################################################

#SETUP
################################################################################
#working folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#functions
source("functions.R")

#libraries
library(rioja)
library(dplyr)
library(ggplot2)
library(cowplot) #neat publishing theme for ggplot
library(ggdendro) #to plot dendrograms with ggplot
library(viridis)
library(nlme)
library(tidyr)

#importing data
load("data.RData")


#0. SUBSETTING DATA
###############################################################################
#removing first 60 samples
bale.data <- bale.data[61:nrow(bale.data), ]
#removing ages older than 13.7k
bale.data <- bale.data[bale.data$age <= 13.779,]



#1. TRANSFER FUNCTIONS TO COMPUTE BURNED AREA
################################################################################
#subsetting bale.data and removing NAs
bale.charcoal <- na.omit(bale.data[, c("age", "charcoal.acc.rate")])

#computing fire measures
fire <- data.frame(
  fires.per.year=FN(bale.charcoal$charcoal.acc.rate), 
  fire.radiative.power=FRP(bale.charcoal$charcoal.acc.rate), 
  burned.area=BA(bale.charcoal$charcoal.acc.rate)
  )

#joining with bale charcoal to plot it later
bale.charcoal <- data.frame(bale.charcoal, fire)

#adding charcoal to fire
fire$charcoal.acc.rate <- bale.charcoal$charcoal.acc.rate

#Transformation to correct for skewness
fire.log <- log10(fire + 1)

#Standardising charcoal and fire variables for clustering
fire.scaled <- data.frame(scale(fire.log))

#Computing distance matrix for clustering only with charcoal and BA
fire.dist <- dist(fire.scaled[, c("charcoal.acc.rate", "burned.area")])

#Constrained hierarchical clustering
fire.cluster <- chclust(fire.dist, method="coniss")

#computing optimum number of groups by minimizing the sum of dispersion and bstick
fire.cluster.bstick <- bstick(fire.cluster)
fire.cluster.bstick$minimum <- fire.cluster.bstick$dispersion + fire.cluster.bstick$bstick
number.of.groups <- fire.cluster.bstick[which.min(fire.cluster.bstick$minimum), "nGroups"]

#Separate cluster zones
fire.cluster.zones <- cutree(fire.cluster, k=number.of.groups)

#Joining grouping vector with bale.charcoal
bale.charcoal$clustering.zone <- fire.cluster.zones

#Computing fire rotation period for each zone
fire.rotation.period <- vector()
min.age <- vector()
max.age <- vector()
for(zone in unique(fire.cluster.zones)){
  fire.rotation.period <- c(fire.rotation.period, 1000/mean(bale.charcoal[bale.charcoal$clustering.zone==zone, "burned.area"]))
  min.age <- c(min.age, min(bale.charcoal[bale.charcoal$clustering.zone==zone, "age"]))
  max.age <- c(max.age, max(bale.charcoal[bale.charcoal$clustering.zone==zone, "age"]))
}

#to dataframe
fire.rotation.period <- data.frame(zone=unique(fire.cluster.zones), fire.rotation.period=fire.rotation.period, min.age=min.age, max.age=max.age)

#saving to csv
write.table(fire.rotation.period, file="output/fire_rotation_period.csv", col.names = TRUE, row.names = FALSE, sep=";")


#2. CHRONOSTRATIGRAPHICAL DIAGRAM
###############################################################################
plotDiagram(pollen.data=bale.data, fire.data=bale.charcoal, fire.cluster=fire.cluster, dendrogram=FALSE, title="Garba Guracha (3950 masl)", title.size=16, text.size=10, width=12, height=9, filename="output/Figure_1.pdf")

plotDiagram(pollen.data=bale.data, fire.data=bale.charcoal, fire.cluster=fire.cluster, dendrogram=TRUE, title="Garba Guracha (3950 masl)", title.size=16, text.size=10, width=12, height=9, filename="output/Figure_1_with_dendrogram.pdf")

#removing objects we don't need any longer
rm(fire, fire.cluster, fire.log, fire.rotation.period, fire.scaled, fire.cluster.zones, fire.dist, min.age, max.age, zone, number.of.groups, fire.cluster.bstick)
gc()


#3. FEEDBACK BETWEEN FIRE AND ERICA
###############################################################################

#3.1 CORRELATION BETWEEN ERICA AND CHAR
#--------------------------------------
#removing NAs
erica.char<-na.omit(bale.data[, c("ericaceae.par", "charcoal.acc.rate")])
erica.char <- data.frame(scale(erica.char))

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


#3.2 PREPARING LAGGED DATA
#--------------------------------------
#subsetting erica and charcoal and rounding ages
erica <- na.omit(bale.data[, c("age","ericaceae.par")])
erica$age <- round(erica$age, 3)

char <- na.omit(bale.data[, c("age", "charcoal.acc.rate")])
char$age <- round(char$age, 3)

#interpolating charcoal to a regular grid (10 years intervals)
char.interpolated<-interpolateDatasets(datasets.list=list(char=char), age.column.name="age", interpolation.time.step=1/100)
char.interpolated$age<-round(char.interpolated$age, 3)

#plot to compare observed and interpolated char curves
ggplot() + 
  geom_line(data=char, aes(x=age, y=charcoal.acc.rate), col="gray50", size=1) +
  geom_line(data=char.interpolated, aes(x=age, y=charcoal.acc.rate), color="red4", size=0.5, alpha=0.5) +
  ylab("CHAR (particles/cm2 yr)") +
  scale_x_continuous(breaks=seq(0, ceiling(max(bale.data$age)), by=1)) +
  xlab("Age (ka BP)") +
  ggtitle("CHAR: observed (gray) vs. interpolated (red)")
ggsave(filename="output/observed_vs_interpolated_charcoal.pdf", width=12, height=4)

#erica ages should be within char.interpolated ages
erica <- erica[erica$age >= min(char.interpolated$age), ]
erica <- erica[erica$age <= max(char.interpolated$age), ]

#generates lagged charcoal data after every erica sample, up to 25 lags (250 years)
lags<-1:101
lag.data.forward <- forwardLags(lags=lags, reference.data=erica, data.to.lag=char.interpolated)
lag.data.backward <- backwardLags(lags=lags, reference.data=erica, data.to.lag=char.interpolated)

#scatterplots of lagged datasets
plotLags(forward.lagged.data=lag.data.forward, backward.lagged.data=lag.data.backward, lags=lags, filename="output/lagged_data_scatterplots.pdf")


#3.3 ANALYSING LAGGED DATA
#--------------------------------------

#fitting a GLS model per lag
backward.results <- modelLagData(model.formula="ericaceae.par ~ charcoal.acc.rate", lagged.data=lag.data.backward)
forward.results <- modelLagData(model.formula="charcoal.acc.rate ~ ericaceae.par", lagged.data=lag.data.forward)

#computing GLS models with randomized responses (this can take a while depending on the number of iterations)
backward.results.random <- modelRandomLagData(lagged.data=lag.data.backward, model.formula="ericaceae.par ~ charcoal.acc.rate", iterations=100)
forward.results.random <- modelRandomLagData(lagged.data=lag.data.forward, model.formula="charcoal.acc.rate ~ ericaceae.par", iterations=100)

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



