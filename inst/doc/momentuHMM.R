## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----load-library, echo=TRUE, eval=FALSE---------------------------------
#  library(momentuHMM)
#  ### Load raw data
#  rawHaggis<-read.csv("rawHaggises.csv")
#  ### Process data
#  processedHaggis<-prepData(data=rawHaggis,covNames=c("slope","temp"))
#  
#  ### Fit HMM	
#  # initial step distribution natural scale parameters
#  stepPar0 <- c(1,5,0.5,3) # (mu_1,mu_2,sd_1,sd_2)
#  # initial angle distribution natural scale parameters
#  anglePar0 <- c(0,0,1,8) # (mean_1,mean_2,concentration_1,concentration_2)       	
#  fitHaggis <- fitHMM(data = processedHaggis, nbStates = 2,
#                      dist = list(step = "gamma", angle = "vm"),
#                      Par0 = list(step = stepPar0, angle = anglePar0),
#                      formula = ~ slope + I(slope^2),
#                      estAngleMean = list(angle=TRUE))

## ----fit-haggis, echo=TRUE, eval=FALSE-----------------------------------
#  stepDM <- list(mean = ~1, sd = ~1)
#  
#  ### Fit HMM	using user-specified DM
#  fitHaggisDM <- fitHMM(data = processedHaggis, nbStates = 2,
#                        dist = list(step = "gamma", angle = "vm"),
#                        DM = list(step = stepDM),
#                        Par0 = list(step = log(stepPar0), angle = anglePar0),
#                        formula = ~ slope + I(slope^2),
#                        estAngleMean = list(angle=TRUE))

## ----fit-haggis-2, echo=TRUE, eval=FALSE---------------------------------
#  stepDMp <- matrix(c(1,0,0,0,
#                     0,1,0,0,
#                     0,0,1,0,
#                     0,0,0,1),4,4,byrow=TRUE)
#  rownames(stepDMp) <- c("mean_1","mean_2","sd_1","sd_2")
#  colnames(stepDMp) <- c("mean_1:(Intercept)","mean_2:(Intercept)",
#                        "sd_1:(Intercept)","sd_2:(Intercept)")
#  
#  ### Fit HMM	using user-specified DM
#  fitHaggisDMp <- fitHMM(data = processedHaggis, nbStates = 2,
#                         dist = list(step = "gamma", angle = "vm"),
#                         DM = list(step = stepDMp),
#                         Par0 = list(step = log(stepPar0), angle = anglePar0),
#                         formula = ~ slope + I(slope^2),
#                         estAngleMean = list(angle=TRUE))

## ----load-vignette_inputs, results='hide', echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE----
library(momentuHMM)
#example_wd <- ("~/Dropbox/vignette examples/")
example_wd <- ("~/Documents/Dropbox/current projects/moveHMM extension/momentuHMM/vignette examples/")

#source("vignette_examples.R")
load("vignette_inputs.RData")

## ----load-data, echo=TRUE, eval=TRUE-------------------------------------
URL <- paste0("https://www.datarepository.movebank.org/bitstream/handle/",
              "10255/move.373/Elliptical%20Time-Density%20Model%20%28Wall%",
              "20et%20al.%202014%29%20African%20Elephant%20Dataset%20%",
              "28Source-Save%20the%20Elephants%29.csv")
rawData <- read.csv(url(URL))

## ----subset-data, echo=TRUE----------------------------------------------
# select and rename relevant columns
rawData <- rawData[,c(11,3,4,5,6)]
colnames(rawData) <- c("ID","time","lon","lat","temp")

# only keep first track
rawData <- subset(rawData,ID==unique(ID)[1])

## ----head-data, echo=TRUE------------------------------------------------
head(rawData)

## ----project-data, echo=TRUE, eval=FALSE---------------------------------
#  # convert times from factors to POSIX
#  rawData$time <- as.POSIXct(rawData$time,tz="GMT")
#  
#  # project to UTM coordinates using package rgdal
#  library(rgdal)
#  llcoord <- SpatialPoints(rawData[,3:4],
#                           proj4string=CRS("+proj=longlat +datum=WGS84"))
#  utmcoord <- spTransform(llcoord,CRS("+proj=utm +zone=30 ellps=WGS84"))
#  
#  # add UTM locations to data frame
#  rawData$x <- attr(utmcoord,"coords")[,1]
#  rawData$y <- attr(utmcoord,"coords")[,2]

## ----crawl, echo=TRUE, eval=FALSE----------------------------------------
#  # initial parameters for crawl fit
#  inits <- list(a = c(rawData$x[1],0,rawData$y[1],0),
#                P = diag(c(5000^2, 10*3600^2, 5000^2, 10*3600^2)))
#  
#  # fit crawl model
#  crwOut <- crawlWrap(obsData=rawData, timeStep="hour", initial.state=inits,
#                      theta=c(4,-10), fixPar=c(NA,NA))

## ----prepData, echo=TRUE, eval=FALSE-------------------------------------
#  # create momentuHMMData object from crwData object
#  elephantData <- prepData(data=crwOut, covNames="temp")
#  
#  # add cosinor covariate based on hour of day
#  elephantData$hour <- as.integer(strftime(elephantData$time, format = "%H", tz="GMT"))

## ----acf-before, echo=TRUE, eval=FALSE-----------------------------------
#  acf(elephantData$step[!is.na(elephantData$step)],lag.max=300)

## ----fitm1, echo=TRUE, eval=FALSE----------------------------------------
#  # label states
#  stateNames <- c("encamped","exploratory")
#  # distributions for observation processes
#  dist = list(step = "gamma", angle = "wrpcauchy")
#  
#  # initial parameters
#  Par0_m1 <- list(step=c(100,500,100,200),angle=c(0.3,0.7))
#  
#  # fit model
#  m1 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m1,
#               estAngleMean = list(angle=FALSE), stateNames = stateNames)

## ----fitm2, echo=TRUE, eval=FALSE----------------------------------------
#  # formula for transition probabilities
#  formula <- ~ temp * cosinor(hour, period = 24)
#  
#  # initial parameters (obtained from nested model m1)
#  Par0_m2 <- getPar0(model=m1, formula=formula)
#  
#  # fit model
#  m2 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
#               beta0=Par0_m2$beta, stateNames = stateNames, formula=formula)

## ----fitm3, echo=TRUE, eval=FALSE----------------------------------------
#  # formulas for parameters of state-dependent observation distributions
#  DM <- list(step = list(mean = ~ temp * cosinor(hour, period = 24),
#                         sd = ~ temp * cosinor(hour, period = 24)),
#             angle = list(concentration = ~ temp))
#  
#  # initial parameters (obtained from nested model m2)
#  Par0_m3 <- getPar0(model=m2, formula=formula, DM=DM)
#  
#  # fit model
#  m3 <- fitHMM(data = elephantData, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
#               beta0 = Par0_m3$beta, DM = DM, stateNames = stateNames,
#               formula = formula)

## ----viterbi, echo=TRUE, eval=FALSE--------------------------------------
#  # decode most likely state sequence
#  states <- viterbi(m3)
#  # derive percentage of time spent in each state
#  table(states)/nrow(elephantData)

## ----plotm3, echo=TRUE, eval=FALSE---------------------------------------
#  plot(m3, plotCI = TRUE, covs = data.frame(hour=12))

## ----pseudoRes, echo=TRUE, eval=FALSE------------------------------------
#  # compute pseudo-residuals for the steps and the angles
#  pr <- pseudoRes(m3)
#  
#  # plot the ACF of step pseudo-residuals
#  acf(pr$stepRes[!is.na(pr$stepRes)],lag.max = 300)

## ----fit-nfs, echo=TRUE, eval=FALSE--------------------------------------
#  nbStates <- 3
#  stateNames <- c("resting", "foraging", "transit")
#  dist <- list(step = "gamma", angle = "wrpcauchy", dive = "pois")
#  Par0 <- getParDM(nbStates = nbStates, dist = dist,
#                   Par = Par, DM = DM, workBounds = workBounds,
#                   estAngleMean = list(angle = FALSE))
#  Fixpar <- list(dive = c(-100, NA, NA))
#  nfsFits <- MIfitHMM(crwOut, nSims = 100, nbStates = nbStates, dist = dist,
#                      Par0 = Par0, DM = DM, workBounds = workBounds,
#                      estAngleMean = list(angle = FALSE),
#                      fixPar = fixPar, retryFits = 30,
#                      stateNames=stateNames)
#  plot(nfsFits)

## ----fit-turtle, echo=TRUE, eval=FALSE-----------------------------------
#  miTurtleData <- MIfitHMM(crwOut, nSims = 100, fit=FALSE,
#                   spatialCovs = list(w = speedBrick, d = dirBrick, r = dirBrick),
#                   angleCovs = "d")

## ----fit-turtle-2, echo=TRUE, eval=FALSE---------------------------------
#  nbStates<-2
#  dist <- list(step = "gamma", angle = "wrpcauchy")
#  DM <- list(step = list(mean = ~state2(w:angle_osc), sd = ~1),
#             angle = list(mean = ~state2(d), concentration= ~1))
#  turtleFits <- MIfitHMM(miTurtleData$miData, nbStates = nbStates, dist = dist,
#                         Par0 = Par0, DM = DM,
#                         estAngleMean = list(angle = TRUE),
#                         circularAngleMean = list(angle = TRUE))
#  plot(turtleFits, plotCI = TRUE, covs = data.frame(angle_osc = cos(0)))

## ----fit-turtle-3, echo=TRUE, eval=FALSE---------------------------------
#  DM$angle = list(mean = ~state2(angleFormula(d, strength = w)),
#                  concentration= ~1))

## ----fit-turtle-4, echo=TRUE, eval=FALSE---------------------------------
#  dist$angle = "vmConsensus"
#  DM$angle = list(mean = ~state2(angleFormula(d, strength = w)),
#                  kappa = ~1)

## ----fit-grey-seal, echo=TRUE, eval=FALSE--------------------------------
#  crwSim <- MIfitHMM(crwOut, nSims = 100, fit=FALSE,
#                    center = centers)

## ----spec-grey-seal, echo=TRUE, eval=FALSE-------------------------------
#  dist <- list(step = "weibull", angle = "wrpcauchy")
#  distFormula <- ~state1(I(Abertay.dist>2500)) + state2(I(Farne.dist>2500))
#                    + state3(I(Dogger.dist>15000))
#  angleFormula <- ~state1(Abertay.angle) + state2(Farne.angle)
#                    + state3(Dogger.angle)
#  stepDM <- list(shape = distFormula, scale = distFormula)
#  angleDM <- list(mean = angleFormula, concentration = distFormula)
#  DM <- list(step = stepDM, angle = angleDM)

## ----fit-grey-seal-2, echo=TRUE, eval=FALSE------------------------------
#  greySealFits <- MIfitHMM(miDat, nSims = 400,
#                           nbStates = 5, dist = dist,
#                           Par0 = Par0, beta0 = beta0, fixPar = fixPar,
#                           formula = distFormula,
#                           estAngleMean = list(angle=TRUE),
#                           circularAngleMean = list(angle=TRUE),
#                           DM = DM, knownStates = knownStates)
#  plot(greySealFits, plotCI = TRUE)

## ----sim-grey-seal, echo=TRUE, eval=FALSE--------------------------------
#  greySealSim<-simData(model = greySealFits, centers = centers,
#                       initialPosition = centers[1,],
#                       obsPerAnimal = 1515)

## ----head-ses, echo=TRUE, eval=TRUE--------------------------------------
head(tracks)

## ----prep-ses, echo=TRUE-------------------------------------------------
center <- matrix(c(70,-49),nrow=1,dimnames=list("colony"))
data <- prepData(data=tracks, type="LL", centers=center)

## ----m1-ses, echo=TRUE, eval=FALSE---------------------------------------
#  stateNames <- c("outbound","search","forage","inbound")
#  
#  # initial parameters
#  stepPar0 <- c(25,5,1,25,10,5,3,10)
#  anglePar0 <- c(15,5,2,15)
#  
#  # constrain transition probabilities
#  fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100),
#                    nrow=1)
#  
#  m1 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"),
#               Par0=list(step=stepPar0, angle=anglePar0),
#               fixPar=list(beta=fixbeta), stateNames = stateNames)

## ----formula-ses, echo=TRUE, eval=FALSE----------------------------------
#  # time spent since left colony
#  time <- NULL
#  for(id in unique(data$ID)) {
#      nbSubObs <- length(which(data$ID==id))
#  
#      # approximately in months for interval = 9.6h
#      time <- c(time, (1:nbSubObs)/75)
#  }
#  
#  data$time <- time
#  
#  # compute time since departure and include in formula below
#  formula <- ~ colony.dist + time

## ----fixpar-ses, echo=TRUE, eval=FALSE-----------------------------------
#  fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100,
#                      NA,   0,   0,   0, 0, 0,   0, 0,   0,   0,   0,   0,
#                       0,   0,   0,   0, 0,NA,   0, 0,   0,   0,   0,   0),
#                    nrow=3,byrow=TRUE)

## ----crw-ses-1, echo=TRUE, eval=FALSE------------------------------------
#  angleFormula <- ~ state1(colony.angle) + state4(colony.angle)

## ----crw-ses-2, echo=TRUE, eval=FALSE------------------------------------
#  fixPar <- list(angle=c(-100,100,NA,NA,NA,NA),beta=fixbeta)

## ----fit-ses, echo=TRUE, eval=FALSE--------------------------------------
#  Par0 <- getPar0(model=m1, nbStates=4,
#                  DM=list(angle=list(mean=angleFormula, concentration=~1)),
#                  estAngleMean=list(angle=TRUE),
#                  circularAngleMean=list(angle=TRUE), formula=formula)
#  
#  m2 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"),
#               Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
#               beta0=Par0$beta, fixPar=fixPar, formula=formula,
#               DM=list(angle=list(mean=angleFormula, concentration=~1)),
#               estAngleMean=list(angle=TRUE), circularAngleMean=list(angle=TRUE),
#               stateNames = stateNames)

## ----crw-ses-3, echo=TRUE, eval=FALSE------------------------------------
#  formula <- ~ betaCol1(colony.dist) + betaCol6(time)
#  
#  fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100,
#                      rep(NA,12),
#                      rep(NA,12)),
#                    nrow=3,byrow=TRUE)
#  
#  fixPar <- list(angle=c(-100,100,NA,NA,NA,NA),beta=fixbeta)

## ----mod3-ses, echo=TRUE, eval=FALSE-------------------------------------
#  distFormula <- ~ state1(colony.dist) + state4(colony.dist)
#  stepDM <- list(mean=distFormula, sd=distFormula)
#  angleDM <- list(mean=angleFormula, concentration=distFormula)

## ----fit-ses-2, echo=TRUE, eval=FALSE------------------------------------
#  # remove fixed angle parameters
#  fixPar <- list(beta=fixbeta)
#  
#  # get starting parameters from m2
#  Par0 <- getPar0(model=m2, nbStates=4,
#                  DM = list(step=stepDM, angle=angleDM),
#                  estAngleMean=list(angle=TRUE),
#                  circularAngleMean=list(angle=TRUE),
#                  formula=formula)
#  
#  # the bias is estimated rather than fixed
#  Par0$Par$angle[c("mean_1:(colony.angle)","mean_4:(colony.angle)")] <- 0
#  
#  m3 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"),
#               Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
#               beta0=Par0$beta, fixPar=fixPar, formula=formula,
#               DM = list(step=stepDM, angle=angleDM),
#               estAngleMean=list(angle=TRUE),
#               circularAngleMean=list(angle=TRUE),
#               stateNames = stateNames)

## ----groupModel-centroid, echo=TRUE, eval=FALSE--------------------------
#  dist <- list(step="gamma", angle="vm")
#  nbObs <- 250
#  
#  Parc <- list(step = c(15,10),
#               angle = c(0.15,log(1)))
#  
#  DMc <- list(angle=list(mean = ~center1.angle,
#                         concentration=~1))
#  
#  centroidData <- simData(nbStates=1, dist=dist, Par=Parc, DM=DMc,
#                          circularAngleMean = list(angle = TRUE),
#                          centers = matrix(0,1,2),
#                          obsPerAnimal = nbObs)

## ----groupModel-groupData, echo=TRUE, eval=FALSE-------------------------
#  nbAnimals <- 20
#  nbStates <- 2
#  stateNames <- c("group","solitary")
#  
#  Par <- list(step = c(30,50,15,25),
#              angle = c(1.e+7,log(2.5),log(5)))
#  
#  beta <- matrix(c(-2.944439,-1.734601),1,nbStates)
#  
#  DM <- list(angle=list(mean = ~state1(centroid.angle),
#                        concentration = ~1))
#  
#  # calculate stationary distribution
#  gamma <- diag(nbStates)
#  gamma[!gamma] <- exp(beta)
#  gamma <- t(gamma)
#  gamma <- gamma/apply(gamma,1,sum)
#  delta <- solve(diag(nbStates) - t(gamma) + 1, rep(1, nbStates))
#  
#  # draw random initial locations for each individual
#  initialPositions <- vector("list")
#  for (i in 1:nbAnimals) {
#    initialPositions[[i]] <- runif(2, -10, 10)
#  }
#  
#  # create centroid data frame
#  cD <- data.frame(x = centroidData$x, y = centroidData$y)
#  
#  groupData <- simData(nbAnimals=nbAnimals, nbStates=nbStates, dist=dist,
#                       Par = Par, beta = beta, delta = delta, DM = DM,
#                       circularAngleMean = list(angle = TRUE),
#                       centroids = list(centroid = cD),
#                       obsPerAnimal = nbObs,
#                       initialPosition = initialPositions,
#                       states = TRUE, stateNames = stateNames)

## ----groupModel-fit, echo=TRUE, eval=FALSE-------------------------------
#  Par0 <- list(step = c(30,50,15,25),
#               angle = c(0,log(2.5),log(5)))
#  
#  groupFit <- fitHMM(groupData, nbStates=nbStates, dist=dist, Par=Par0,
#                     DM = DM, stationary = TRUE,
#                     estAngleMean = list(angle = TRUE),
#                     circularAngleMean = list(angle = TRUE),
#                     stateNames = stateNames)

## ----spec-hsStep, echo=TRUE, eval=TRUE-----------------------------------
nbStates <- 3
stateNames <- c("resting", "foraging", "transit")
dist <- list(step = "weibull", angle = "wrpcauchy", dive = "beta")
stepDM<-matrix(c(1,0,0,0,0,0,0,0,0,
                 0,1,0,0,0,0,0,0,0,
                 0,0,1,0,0,0,0,0,0,
                 0,0,0,1,0,0,0,0,0,
                 0,0,0,1,1,0,0,0,0,
                 0,0,0,1,1,1,0,0,0,
                 0,0,0,0,0,0,1,0,0,
                 0,0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,0,0,1),nrow=3*nbStates,byrow=TRUE,
        dimnames=list(c(paste0("shape_",1:nbStates),
                        paste0("scale_",1:nbStates),
                        paste0("zeromass_",1:nbStates)),
                      c(paste0("shape_",1:nbStates,":(Intercept)"),
                        "scale:(Intercept)","scale_2","scale_3",
                        paste0("zeromass_",1:nbStates,":(Intercept)"))))
stepworkBounds<-matrix(c(rep(-Inf,4),0,0,rep(-Inf,3),
                   rep(Inf,ncol(stepDM))),ncol(stepDM),2)
stepBounds<-matrix(c(0,5,
                     0,5,
                     0,5,
                     0,14400,
                     0,14400,
                     0,14400,
                     0,1,
                     0,1,
                     0,1),nrow=3*nbStates,byrow=TRUE,
                    dimnames=list(c(paste0("shape_",1:nbStates),
                                    paste0("scale_",1:nbStates),
                                    paste0("zeromass_",1:nbStates))))

## ----spec-hsAngle, echo=TRUE, eval=TRUE----------------------------------
angleDM <- matrix(c(1,0,0,
                    0,1,1,
                    0,1,0),nrow=nbStates,byrow=TRUE,
                  dimnames=list(paste0("concentration_",1:nbStates),
                                c("concentration_1:(Intercept)",
                                  "concentration_23:(Intercept)",
                                  "concentration_2")))
angleBounds <- matrix(c(0,0.95,
                        0,0.95,
                        0,0.95),nrow=nbStates,byrow=TRUE,
                      dimnames=list(paste0("concentration_",1:nbStates)))
transitcons <- boot::logit((0.75 - angleBounds[3,1])
                           /(angleBounds[3,2] - angleBounds[3,1]))
angleworkBounds <- matrix(c(-Inf,transitcons,-Inf,
                            rep(Inf,2),0),ncol(angleDM),2)

## ----spec-hsDive, echo=TRUE, eval=FALSE----------------------------------
#  omegaDM <- matrix(c(1,0,0,0,0,0,
#                      0,0,1,1,0,0,
#                      0,0,1,1,0,0,
#                      1,1,0,0,0,0,
#                      0,0,1,0,0,0,
#                      0,0,1,0,0,0,
#                      0,0,0,0,1,0,
#                      0,0,0,0,0,1,
#                      0,0,0,0,0,1),nrow=nbStates*3,byrow=TRUE,
#                  dimnames=list(c(paste0("shape1_",1:nbStates),
#                                  paste0("shape2_",1:nbStates),
#                                  paste0("zeromass_",1:nbStates)),
#                                c("shape_1:(Intercept)","shape2_1",
#                                  "shape_2:(Intercept)","shape1_2",
#                                  "zeromass_1:(Intercept)",
#                                  "zeromass_23:(Intercept)")))
#  omegaworkBounds <- matrix(c(-Inf,0,-Inf,0,-Inf,-Inf,
#                              rep(Inf,ncol(omegaDM))),ncol(omegaDM),2)
#  omegaBounds <- matrix(c(1,10,
#                          1,10,
#                          1,10,
#                          1,10,
#                          1,10,
#                          1,10,
#                          0,1,
#                          0,1,
#                          0,1),nrow=nbStates*3,byrow=TRUE,
#                      dimnames=list(c(paste0("shape1_",1:nbStates),
#                                      paste0("shape2_",1:nbStates),
#                                      paste0("zeromass_",1:nbStates))))

## ----spec-hsFixPar, echo=TRUE, eval=FALSE--------------------------------
#  fixPar <- list(step=c(rep(NA,nbStates*2),NA,NA,boot::logit(1.e-100)),
#               omega=c(rep(NA,4),NA,boot::logit(1.e-100)))

## ----fit-hs, echo=TRUE, eval=FALSE---------------------------------------
#  DM <- list(step = stepDM, angle = angleDM, omega = omegaDM)
#  userBounds <- list(step = stepBounds,
#                     angle = angleBounds,
#                     omega = omegaBounds)
#  workBounds <- list(step = stepworkBounds,
#                     angle = angleworkBounds,
#                     omega = omegaworkBounds)
#  hsFits <- MIfitHMM(crwOut, nSims = 30,
#                     nbStates = nbStates, dist = dist, Par0 = Par0,
#                     DM = DM, workBounds = workBounds,
#                     userBounds = userBounds, workBounds = workBounds,
#                     fixPar = fixPar, stateNames = stateNames)

## ----prep-fulmar-1, echo=TRUE, eval=FALSE--------------------------------
#  library(sp)
#  
#  # load data provided by Pirotta et al
#  fulmarURL <- "https://datadryad.org/bitstream/handle/10255/dryad.174482/"
#  fileName <- "Fulmar_trackingData.csv?sequence=1"
#  raw_data <- read.csv(url(paste0(fulmarURL,fileName)),
#                       stringsAsFactors = FALSE)
#  
#  raw_data$ID <- raw_data$tripID
#  raw_data$Date <- as.POSIXct(raw_data$Date,tz="UTC",
#                              format="%d/%m/%Y %H:%M")
#  
#  # project data
#  oldProj <- CRS("+proj=longlat +datum=WGS84")
#  newProj <- CRS("+init=epsg:27700")
#  coordinates(raw_data) <- c("Longitude","Latitude")
#  proj4string(raw_data) <- oldProj
#  raw_data <- as.data.frame(spTransform(raw_data, newProj))
#  
#  coordinates(raw_data) <- c("Boat_Longitude","Boat_Latitude")
#  proj4string(raw_data) <- oldProj
#  raw_data <- as.data.frame(spTransform(raw_data, newProj))

## ----prep-fulmar-2, echo=TRUE, eval=FALSE--------------------------------
#  # use prepData to calculate colony distance covariate ('sea.angle')
#  colony <- data.frame(x = -3.1, y = 59.12)
#  coordinates(colony) <- c("x", "y")
#  proj4string(colony) <- oldProj
#  colony <- as.matrix(as.data.frame(spTransform(colony, newProj)))
#  rownames(colony) <- "colony"
#  colony_dist <- prepData(raw_data, coordNames = c("Longitude","Latitude"),
#                          centers = colony)
#  
#  # calculate "sea" mean angle covariate
#  sea.angle <- NULL
#  for(id in unique(colony_dist$ID)) {
#    idat <- subset(colony_dist,ID==id)
#    nbSubObs <- length(which(colony_dist$ID==id))
#    max_dist <- as.numeric(idat[which.max(idat$colony.dist),c("x","y")])
#    max_angle <- momentuHMM:::distAngle(colony,colony,max_dist)[2]
#    sea.angle <- c(sea.angle, rep(max_angle,nbSubObs))
#  }
#  raw_data$sea.angle <- sea.angle

## ----prep-fulmar-3, echo=TRUE, eval=FALSE--------------------------------
#  # calculate time since left colony covariate ('time')
#  time <- aInd <- NULL
#  for(id in unique(raw_data$ID)) {
#    idInd <- which(raw_data$ID==id)
#    aInd <- c(aInd,idInd[1])
#    nbSubObs <- length(idInd)
#    time <- c(time, (1:nbSubObs)/nbSubObs)
#  }
#  raw_data$time <- time

## ----prep-fulmar-4, echo=TRUE, eval=FALSE--------------------------------
#  # get boat data into centroids argument format
#  boat_data <- list(boat=data.frame(Date = raw_data$Date,
#                                    x = raw_data$Boat_Longitude,
#                                    y = raw_data$Boat_Latitude))
#  
#  # format and merge all data and covariates for analysis
#  fulmar_data <- prepData(raw_data, coordNames = c("Longitude","Latitude"),
#                          centers = colony,
#                          centroids = boat_data,
#                          covNames = "time",
#                          angleCovs = "sea.angle")
#  
#  # momentuHMM doesn't like data streams and covariates to have same name,
#  # so create identical data column with different name
#  fulmar_data$d <- fulmar_data$boat.dist
#  
#  # standarize boat.dist covariate
#  fulmar_data$boat.dist <- scale(fulmar_data$boat.dist)

## ----spec-fulmar-1, echo=TRUE, eval=FALSE--------------------------------
#  nbStates <- 6
#  
#  stateNames <- c("seaARS", "seaTr",
#                  "boatARS", "boatTr",
#                  "colonyARS", "colonyTr")
#  
#  dist <- list(step = "weibull",
#               angle = "wrpcauchy",
#               d = "lnorm")
#  

## ----spec-fulmar-2, echo=TRUE, eval=FALSE--------------------------------
#  # specify data stream probability distribution parameter constraints
#  stepDM <- matrix(c(1,0,0,0,
#                     0,1,0,0,
#                     1,0,0,0,
#                     0,1,0,0,
#                     1,0,0,0,
#                     0,1,0,0,
#                     0,0,1,0,
#                     0,0,1,1,
#                     0,0,1,0,
#                     0,0,1,1,
#                     0,0,1,0,
#                     0,0,1,1),2*nbStates,4,byrow=TRUE,
#                   dimnames=list(c(paste0("shape_",1:nbStates),
#                                   paste0("scale_",1:nbStates)),
#                                 c("shape:ARS","shape:Tr",
#                                   "scale:(Intercept)","scale:Tr")))
#  
#  # constrain scale parameters such that Tr > ARS
#  stepworkBounds <- matrix(c(-Inf,Inf,
#                             -Inf,Inf,
#                             -Inf,Inf,
#                             0,Inf),ncol(stepDM),2,byrow=TRUE)
#  
#  
#  nbTrips <- length(unique(fulmar_data$ID))
#  
#  angleDM <- matrix(c("sea.angle",0,0,0,0,rep(0,2*nbTrips),
#                      "sea.angle",0,0,0,0,rep(0,2*nbTrips),
#                      0,"boat.angle",0,0,0,rep(0,2*nbTrips),
#                      0,"boat.angle",0,0,0,rep(0,2*nbTrips),
#                      0,0,"colony.angle",0,0,rep(0,2*nbTrips),
#                      0,0,"colony.angle",0,0,rep(0,2*nbTrips),
#                      0,0,0,1,0,paste0("ID",1:nbTrips),rep(0,nbTrips),
#                      0,0,0,1,1,paste0("ID",1:nbTrips),paste0("ID",1:nbTrips),
#                      0,0,0,1,0,rep(0,2*nbTrips),
#                      0,0,0,1,1,rep(0,2*nbTrips),
#                      0,0,0,1,0,rep(0,2*nbTrips),
#                      0,0,0,1,1,rep(0,2*nbTrips)),2*nbStates,3+2+2*nbTrips,byrow=TRUE,
#                    dimnames=list(c(paste0("mean_",1:nbStates),
#                                    paste0("concentration_",1:nbStates)),
#                                  c("mean:sea","mean:boat","mean:colony",
#                                    "concentration:(Intercept)","concentration:Tr",
#                                    paste0("concentration:ID",1:nbTrips,":(Intercept)"),
#                                    paste0("concentration:ID",1:nbTrips,":Tr"))))
#  
#  # constrain concentration parameters such that Tr > ARS
#  angleworkBounds <- matrix(c(-Inf,Inf,
#                              -Inf,Inf,
#                              -Inf,Inf,
#                              -Inf,Inf,
#                              0,Inf,
#                              rep(c(-Inf,Inf),nbTrips),
#                              rep(c(0,Inf),nbTrips)),ncol(angleDM),2,byrow=TRUE)
#  
#  dDM <- matrix(c(1,1,0,0,
#                  1,1,0,0,
#                  1,0,0,0,
#                  1,0,0,0,
#                  1,1,0,0,
#                  1,1,0,0,
#                  0,0,1,1,
#                  0,0,1,1,
#                  0,0,1,0,
#                  0,0,1,0,
#                  0,0,1,1,
#                  0,0,1,1),2*nbStates,4,byrow=TRUE,
#                dimnames=list(c(paste0("location_",1:nbStates),
#                                paste0("scale_",1:nbStates)),
#                              c("location:(Intercept)","location:noboat",
#                                "scale:(Intercept)","scale:noboat")))
#  
#  # constrain location and scale parameters such that sea and colony > boat
#  dworkBounds <- matrix(c(-Inf,Inf,
#                          0,Inf,
#                          -Inf,Inf,
#                          0,Inf),ncol(dDM),2,byrow=TRUE)
#  
#  DM<-list(step = stepDM, angle = angleDM, d = dDM)
#  
#  workBounds <- list(step = stepworkBounds,
#                     angle = angleworkBounds,
#                     d = dworkBounds)

## ----spec-fulmar-3, echo=TRUE, eval=FALSE--------------------------------
#  # state transition formula similar to Pirotta et al
#  formula <- ~ toState3(boat.dist) + toState4(boat.dist) +
#               toState5(time) + toState6(time)
#  
#  # specify knownStates
#  # Priotta et al assumed all animals start in state 2 ('seaTr')
#  knownStates <- rep(NA,nrow(fulmar_data))
#  knownStates[aInd] <- 2
#  
#  # fix delta_2 = 1 because assuming initial state is known for each track
#  fixPar <- list(delta=c(100,rep(0,nbStates-2)))
#  fixPar$delta <- exp(c(0,fixPar$delta))/sum(exp(c(0,fixPar$delta)))
#  # Constrain model to BRW (instead of BCRW)
#  fixPar$angle <- c(rep(1.e+7, 3), rep(NA, 2+2*nbTrips))

## ----fit-fulmar-1, echo=TRUE, eval=FALSE---------------------------------
#  prior <- function(par) {sum(dnorm(par[27+1:90],0,100,log=TRUE))}
#  
#  m2 <- fitHMM(fulmar_data, nbStates, dist,
#               Par0 = Par0$Par, beta0 = Par0$beta0,
#               formula = formula,
#               estAngleMean = list(angle = TRUE),
#               circularAngleMean = list(angle = TRUE),
#               DM = DM, workBounds = workBounds,
#               fixPar = fixPar, knownStates = knownStates,
#               stateNames = stateNames,
#               prior = prior)

## ----results-fulmar-1a, echo=TRUE, eval=FALSE----------------------------
#  timeInStates(m2)

## ----results-fulmar-1b, echo=FALSE---------------------------------------
timeIn1

## ----results-fulmar-2a, echo=TRUE, eval=FALSE----------------------------
#  timeInStates(m2, by = "birdID")

## ----results-fulmar-2b, echo=FALSE---------------------------------------
timeIn2

## ----results-fulmar-3, echo=TRUE, eval=FALSE-----------------------------
#  plotSat(m2, zoom = 7, shape = c(17,1,17,1,17,1), size = 2,
#          col = rep(c("#E69F00", "#56B4E9", "#009E73"), each = 2),
#          stateNames = c("sea ARS", "sea Transit",
#                         "boat ARS", "boat Transit",
#                         "colony ARS", "colony Transit"),
#          projargs = newProj, ask = FALSE)

