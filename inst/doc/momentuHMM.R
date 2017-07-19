### R code from vignette source 'momentuHMM.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: load-library (eval = FALSE)
###################################################
## library(momentuHMM)
## ### Load raw data
## rawHaggis<-read.csv("rawHaggises.csv")
## ### Process data
## processedHaggis<-prepData(data=rawHaggis,covNames=c("slope","temp"))
## 
## ### Fit HMM	
## # initial step distribution natural scale parameters
## stepPar0 <- c(1,5,0.5,3) # (mu_1,mu_2,sd_1,sd_2)
## # initial angle distribution natural scale parameters 
## anglePar0 <- c(0,0,1,8) # (mean_1,mean_2,concentration_1,concentration_2)       	
## fitHaggis <- fitHMM(data = processedHaggis, nbStates = 2,
##                     dist = list(step = "gamma", angle = "vm"),
##                     Par0 = list(step = stepPar0, angle = anglePar0),
##                     formula = ~ slope + I(slope^2),
##                     estAngleMean = list(angle=TRUE))


###################################################
### code chunk number 2: fit-haggis (eval = FALSE)
###################################################
## stepDM <- list(mean = ~1, sd = ~1)
## 
## ### Fit HMM	using user-specified DM
## fitHaggisDM <- fitHMM(data = processedHaggis, nbStates = 2,
##                       dist = list(step = "gamma", angle = "vm"),
##                       DM = list(step = stepDM),
##                       Par0 = list(step = log(stepPar0), angle = anglePar0),
##                       formula = ~ slope + I(slope^2),
##                       estAngleMean = list(angle=TRUE))


###################################################
### code chunk number 3: fit-haggis-2 (eval = FALSE)
###################################################
## stepDMp <- matrix(c(1,0,0,0,
##                    0,1,0,0,
##                    0,0,1,0,
##                    0,0,0,1),4,4,byrow=TRUE)
## rownames(stepDMp) <- c("mean_1","mean_2","sd_1","sd_2")
## colnames(stepDMp) <- c("mean_1:(Intercept)","mean_2:(Intercept)",
##                       "sd_1:(Intercept)","sd_2:(Intercept)")
## 
## ### Fit HMM	using user-specified DM
## fitHaggisDMp <- fitHMM(data = processedHaggis, nbStates = 2,
##                        dist = list(step = "gamma", angle = "vm"),
##                        DM = list(step = stepDMp),
##                        Par0 = list(step = log(stepPar0), angle = anglePar0),
##                        formula = ~ slope + I(slope^2),
##                        estAngleMean = list(angle=TRUE))


###################################################
### code chunk number 4: load-vignette_inputs
###################################################
library(momentuHMM)
#example_wd <- ("~/Dropbox/vignette examples/")
example_wd <- ("~/Documents/Dropbox/current projects/moveHMM extension/momentuHMM/vignette examples/")

#source("vignette_examples.R")
load("vignette_inputs.RData")


###################################################
### code chunk number 5: load-data
###################################################
URL <- paste0("https://www.datarepository.movebank.org/bitstream/handle/",
              "10255/move.373/Elliptical%20Time-Density%20Model%20%28Wall%",
              "20et%20al.%202014%29%20African%20Elephant%20Dataset%20%",
              "28Source-Save%20the%20Elephants%29.csv")
rawData <- read.csv(url(URL))


###################################################
### code chunk number 6: subset-data
###################################################
# select and rename relevant columns
rawData <- rawData[,c(11,3,4,5,6)]
colnames(rawData) <- c("ID","time","lon","lat","temp")

# only keep first track
rawData <- subset(rawData,ID==unique(ID)[1])


###################################################
### code chunk number 7: head-data
###################################################
head(rawData)


###################################################
### code chunk number 8: project-data (eval = FALSE)
###################################################
## # convert times from factors to POSIX
## rawData$time <- as.POSIXct(rawData$time,tz="GMT")
## 
## # project to UTM coordinates using package rgdal
## library(rgdal)
## llcoord <- SpatialPoints(rawData[,3:4], 
##                          proj4string=CRS("+proj=longlat +datum=WGS84"))
## utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84"))
## 
## # add UTM locations to data frame
## rawData$x <- attr(utmcoord,"coords")[,1]
## rawData$y <- attr(utmcoord,"coords")[,2]


###################################################
### code chunk number 9: crawl (eval = FALSE)
###################################################
## # initial parameters for crawl fit
## inits <- list(a = c(rawData$x[1],0,rawData$y[1],0),
##               P = diag(c(5000^2, 10*3600^2, 5000^2, 10*3600^2)))
## 
## # fit crawl model
## crwOut <- crawlWrap(obsData=rawData, timeStep="hour", initial.state=inits,
##                     theta=c(4,-10), fixPar=c(NA,NA))


###################################################
### code chunk number 10: prepData (eval = FALSE)
###################################################
## # create momentuHMMData object from crwData object
## elephantData <- prepData(data=crwOut, covNames="temp")
## 
## # add cosinor covariate based on hour of day
## elephantData$hour <- as.integer(strftime(elephantData$time, format = "%H", tz="GMT"))


###################################################
### code chunk number 11: acf-before (eval = FALSE)
###################################################
## acf(elephantData$step[!is.na(elephantData$step)],lag.max=300)


###################################################
### code chunk number 12: fitm1 (eval = FALSE)
###################################################
## # label states
## stateNames <- c("encamped","exploratory")
## # distributions for observation processes
## dist = list(step = "gamma", angle = "wrpcauchy")
## 
## # initial parameters
## Par0_m1 <- list(step=c(100,500,100,200),angle=c(0.3,0.7))
## 
## # fit model
## m1 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m1, 
##              estAngleMean = list(angle=FALSE), stateNames = stateNames)


###################################################
### code chunk number 13: fitm2 (eval = FALSE)
###################################################
## # formula for transition probabilities
## formula <- ~ temp * cosinor(hour, period = 24)
## 
## # initial parameters (obtained from nested model m1)
## Par0_m2 <- getPar0(model=m1, formula=formula)
## 
## # fit model
## m2 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par, 
##              beta0=Par0_m2$beta, stateNames = stateNames, formula=formula)


###################################################
### code chunk number 14: fitm3 (eval = FALSE)
###################################################
## # formulas for parameters of state-dependent observation distributions
## DM <- list(step = list(mean = ~ temp * cosinor(hour, period = 24),
##                        sd = ~ temp * cosinor(hour, period = 24)),
##            angle = list(concentration = ~ temp))
## 
## # initial parameters (obtained from nested model m2)
## Par0_m3 <- getPar0(model=m2, formula=formula, DM=DM)
## 
## # fit model
## m3 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par, 
##              beta0 = Par0_m3$beta, DM = DM, stateNames = stateNames,
##              formula = formula)


###################################################
### code chunk number 15: viterbi (eval = FALSE)
###################################################
## # decode most likely state sequence
## states <- viterbi(m3)
## # derive percentage of time spent in each state
## table(states)/nrow(elephantData)


###################################################
### code chunk number 16: plotm3 (eval = FALSE)
###################################################
## plot(m3, plotCI = TRUE, covs = data.frame(hour=12))


###################################################
### code chunk number 17: pseudoRes (eval = FALSE)
###################################################
## # compute pseudo-residuals for the steps and the angles
## pr <- pseudoRes(m3)
## 
## # plot the ACF of step pseudo-residuals
## acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)


###################################################
### code chunk number 18: fit-nfs (eval = FALSE)
###################################################
## nbStates <- 3
## stateNames <- c("resting", "foraging", "transit")
## dist <- list(step = "gamma", angle = "wrpcauchy", dive = "pois")
## Par0 <- getParDM(nbStates = nbStates, dist = dist,
##                  Par = Par, DM = DM, cons = cons,
##                  estAngleMean = list(angle = FALSE))
## Fixpar <- list(dive = c(-100, NA, NA))
## nfsFits <- MIfitHMM(crwOut, nSims = 100, nbStates = nbStates, dist = dist,
##                     Par0 = Par0, DM = DM, cons = cons,
##                     estAngleMean = list(angle = FALSE), 
##                     fixPar = fixPar, retryFits = 30,
##                     stateNames=stateNames)
## plot(nfsFits)


###################################################
### code chunk number 19: fit-turtle (eval = FALSE)
###################################################
## miTurtleData <- MIfitHMM(crwOut, nSims = 100, fit=FALSE,
##                  spatialCovs = list(w = speedBrick, d = dirBrick, r = dirBrick),
##                  angleCovs = "d")


###################################################
### code chunk number 20: fit-turtle-2 (eval = FALSE)
###################################################
## nbStates<-2
## dist <- list(step = "gamma", angle = "wrpcauchy")
## DM <- list(step = list(mean = ~state2(w:angle_osc), sd = ~1),
##            angle = list(mean = ~state2(d), concentration= ~1))
## turtleFits <- MIfitHMM(miTurtleData$miData, nbStates = nbStates, dist = dist, 
##                        Par0 = Par0, DM = DM, 
##                        estAngleMean = list(angle = TRUE),
##                        circularAngleMean = list(angle = TRUE))
## plot(turtleFits, plotCI = TRUE, covs = data.frame(angle_osc = cos(0)))


###################################################
### code chunk number 21: fit-grey-seal (eval = FALSE)
###################################################
## crwSim <- MIfitHMM(crwOut, nSims = 100, fit=FALSE,
##                   center = centers)


###################################################
### code chunk number 22: spec-grey-seal (eval = FALSE)
###################################################
## dist <- list(step = "weibull", angle = "wrpcauchy")
## distFormula <- ~state1(I(Abertay.dist>2500)) + state2(I(Farne.dist>2500)) 
##                   + state3(I(Dogger.dist>15000))
## angleFormula <- ~state1(Abertay.angle) + state2(Farne.angle) 
##                   + state3(Dogger.angle)
## stepDM <- list(shape = distFormula, scale = distFormula)
## angleDM <- list(mean = angleFormula, concentration = distFormula)
## DM <- list(step = stepDM, angle = angleDM)


###################################################
### code chunk number 23: fit-grey-seal-2 (eval = FALSE)
###################################################
## greySealFits <- MIfitHMM(miDat, nSims = 400,
##                          nbStates = 5, dist = dist,
##                          Par0 = Par0, beta0 = beta0, fixPar = fixPar,
##                          formula = distFormula,
##                          estAngleMean = list(angle=TRUE), 
##                          circularAngleMean = list(angle=TRUE),
##                          DM = DM, knownStates = knownStates)
## plot(greySealFits, plotCI = TRUE)


###################################################
### code chunk number 24: sim-grey-seal (eval = FALSE)
###################################################
## greySealSim<-simData(model = greySealFits, centers = centers,
##                      initialPosition = centers[1,],
##                      obsPerAnimal = 1515)


###################################################
### code chunk number 25: head-ses
###################################################
head(tracks)


###################################################
### code chunk number 26: prep-ses
###################################################
center <- matrix(c(70,-49),nrow=1,dimnames=list("colony"))
data <- prepData(data=tracks, type="LL", centers=center)


###################################################
### code chunk number 27: m1-ses (eval = FALSE)
###################################################
## stateNames <- c("outbound","search","forage","inbound")
## 
## # initial parameters
## stepPar0 <- c(25,5,1,25,10,5,3,10)
## anglePar0 <- c(15,5,2,15)
## 
## # constrain transition probabilities
## fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100),
##                   nrow=1)
## 
## m1 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), 
##              Par0=list(step=stepPar0, angle=anglePar0),
##              fixPar=list(beta=fixbeta), stateNames = stateNames)


###################################################
### code chunk number 28: formula-ses (eval = FALSE)
###################################################
## # time spent since left colony
## time <- NULL
## for(id in unique(data$ID)) {
##     nbSubObs <- length(which(data$ID==id))
##     
##     # approximately in months for interval = 9.6h
##     time <- c(time, (1:nbSubObs)/75)
## }
## 
## data$time <- time
## 
## # compute time since departure and include in formula below
## formula <- ~ colony.dist + time


###################################################
### code chunk number 29: fixpar-ses (eval = FALSE)
###################################################
## fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100,
##                     NA,   0,   0,   0, 0, 0,   0, 0,   0,   0,   0,   0,
##                      0,   0,   0,   0, 0,NA,   0, 0,   0,   0,   0,   0),
##                   nrow=3,byrow=TRUE)


###################################################
### code chunk number 30: crw-ses-1 (eval = FALSE)
###################################################
## angleFormula <- ~ state1(colony.angle) + state4(colony.angle)


###################################################
### code chunk number 31: crw-ses-2 (eval = FALSE)
###################################################
## fixPar <- list(angle=c(-100,100,NA,NA,NA,NA),beta=fixbeta)


###################################################
### code chunk number 32: fit-ses (eval = FALSE)
###################################################
## Par0 <- getPar0(model=m1, nbStates=4, 
##                 DM=list(angle=list(mean=angleFormula, concentration=~1)), 
##                 estAngleMean=list(angle=TRUE), 
##                 circularAngleMean=list(angle=TRUE), formula=formula)
## 
## m2 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), 
##              Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
##              beta0=Par0$beta, fixPar=fixPar, formula=formula, 
##              DM=list(angle=list(mean=angleFormula, concentration=~1)), 
##              estAngleMean=list(angle=TRUE), circularAngleMean=list(angle=TRUE), 
##              stateNames = stateNames)


###################################################
### code chunk number 33: crw-ses-2 (eval = FALSE)
###################################################
## formula <- ~ betaCol1(colony.dist) + betaCol6(time)
## 
## fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100,
##                     rep(NA,12),
##                     rep(NA,12)),
##                   nrow=3,byrow=TRUE)
## 
## fixPar <- list(angle=c(-100,100,NA,NA,NA,NA),beta=fixbeta)


###################################################
### code chunk number 34: mod3-ses (eval = FALSE)
###################################################
## distFormula <- ~ state1(colony.dist) + state4(colony.dist)
## stepDM <- list(mean=distFormula, sd=distFormula)
## angleDM <- list(mean=angleFormula, concentration=distFormula)


###################################################
### code chunk number 35: fit-ses-2 (eval = FALSE)
###################################################
## # remove fixed angle parameters
## fixPar <- list(beta=fixbeta)
## 
## # get starting parameters from m2
## Par0 <- getPar0(model=m2, nbStates=4, 
##                 DM = list(step=stepDM, angle=angleDM), 
##                 estAngleMean=list(angle=TRUE), 
##                 circularAngleMean=list(angle=TRUE), 
##                 formula=formula)
## 
## # the bias is estimated rather than fixed
## Par0$Par$angle[c("mean_1:colony.angle","mean_4:colony.angle")] <- 0
## 
## m3 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), 
##              Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
##              beta0=Par0$beta, fixPar=fixPar, formula=formula, 
##              DM = list(step=stepDM, angle=angleDM), 
##              estAngleMean=list(angle=TRUE), 
##              circularAngleMean=list(angle=TRUE), 
##              stateNames = stateNames)


