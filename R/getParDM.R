#' Get starting values on working scale based on design matrix and other parameter constraints
#' 
#' Convert starting values on the natural scale of data stream probability distributions to
#' a feasible set of working scale parameters based on a design matrix and other parameter constraints.
#' 
#' If design matrix includes non-factor covariates, then natural scale parameters are assumed to correspond to the 
#' mean value(s) for the covariate(s) (if \code{nrow(data)>1}) and \code{getParDM} simply returns one possible solution to the 
#' system of linear equations defined by \code{Par}, \code{DM}, and any other constraints using singular value decomposition. 
#' This can be helpful for exploring relationships between the natural and working scale parameters when covariates are included, but \code{getParDM}
#' will not necessarily return ``good'' starting values (i.e., \code{Par0}) for \code{\link{fitHMM}} or \code{\link{MIfitHMM}}. 
#' 
#' @param data Optional \code{\link{momentuHMMData}} object, \code{\link{momentuHierHMMData}} object, or a data frame containing the covariate values. 
#' \code{data} must be specified if covariates are included in \code{DM}.
#' 
#' If a data frame is provided, then either \code{nbStates} and \code{dist} must be specified (for a regular HMM) or \code{hierStates} and \code{hierDist}
#' must be specified (for a hierarchical HMM).
#' 
#' @param ... further arguments passed to or from other methods
#' @export
getParDM <- function(data, ...) {
  UseMethod("getParDM")
}

#' @rdname getParDM
#' @method getParDM default
#' @param nbStates Number of states of the HMM.
#' @param dist A named list indicating the probability distributions of the data streams. Currently
#' supported distributions are 'bern', 'beta', 'exp', 'gamma', 'lnorm', 'norm', 'mvnorm2' (bivariate normal distribution), 'mvnorm3' (trivariate normal distribution),
#' 'pois', 'rw_norm' (normal random walk), 'rw_mvnorm2' (bivariate normal random walk), 'rw_mvnorm3' (trivariate normal random walk), 'vm', 'vmConsensus', 'weibull', and 'wrpcauchy'. For example,
#' \code{dist=list(step='gamma', angle='vm', dives='pois')} indicates 3 data streams ('step', 'angle', and 'dives')
#' and their respective probability distributions ('gamma', 'vm', and 'pois').
#' @param Par A named list containing vectors of state-dependent probability distribution parameters for 
#' each data stream specified in \code{dist}. The parameters should be on the natural scale,
#' in the order expected by the pdfs of \code{dist}, and any zero-mass parameters should be the last.
#' @param zeroInflation A named list of logicals indicating whether the probability distributions of the data streams should be zero-inflated. If \code{zeroInflation} is \code{TRUE} 
#' for a given data stream, then values for the zero-mass parameters should be
#' included in the corresponding element of \code{Par}. Ignored if \code{data} is a \code{\link{momentuHMMData}} or \code{\link{momentuHierHMMData}} object.
#' @param oneInflation Named list of logicals indicating whether the probability distributions of the data streams are one-inflated. If \code{oneInflation} is \code{TRUE} 
#' for a given data stream, then values for the one-mass parameters should be
#' included in the corresponding element of \code{Par}. Ignored if \code{data} is a \code{\link{momentuHMMData}} or \code{\link{momentuHierHMMData}} object.
#' @param estAngleMean An optional named list indicating whether or not to estimate the angle mean for data streams with angular 
#' distributions ('vm' and 'wrpcauchy'). Any \code{estAngleMean} elements corresponding to data streams that do not have angular distributions are ignored.
#' @param circularAngleMean An optional named list indicating whether to use circular-linear or circular-circular 
#' regression on the mean of circular distributions ('vm' and 'wrpcauchy') for turning angles.  See \code{\link{fitHMM}}. \code{circularAngleMean} elements corresponding to angular data 
#' streams are ignored unless the corresponding element of \code{estAngleMean} is \code{TRUE}. Any \code{circularAngleMean} elements 
#' corresponding to data streams that do not have angular distributions are ignored.
#' @param DM A named list indicating the design matrices to be used for the probability distribution parameters of each data 
#' stream. Each element of \code{DM} can either be a named list of linear regression formulas or a matrix.  For example, for a 2-state 
#' model using the gamma distribution for a data stream named 'step', \code{DM=list(step=list(mean=~cov1, sd=~1))} specifies the mean 
#' parameters as a function of the covariate 'cov1' for each state.  This model could equivalently be specified as a 4x6 matrix using 
#' character strings for the covariate: 
#' \code{DM=list(step=matrix(c(1,0,0,0,'cov1',0,0,0,0,1,0,0,0,'cov1',0,0,0,0,1,0,0,0,0,1),4,6))}
#' where the 4 rows correspond to the state-dependent paramaters (mean_1,mean_2,sd_1,sd_2) and the 6 columns correspond to the regression 
#' coefficients. 
#' @param userBounds An optional named list of 2-column matrices specifying bounds on the natural (i.e, real) scale of the probability 
#' distribution parameters for each data stream. For example, for a 2-state model using the wrapped Cauchy ('wrpcauchy') distribution for 
#' a data stream named 'angle' with \code{estAngleMean$angle=TRUE)}, \code{userBounds=list(angle=matrix(c(-pi,-pi,-1,-1,pi,pi,1,1),4,2))} 
#' specifies (-1,1) bounds for the concentration parameters instead of the default [0,1) bounds.
#' @param workBounds An optional named list of 2-column matrices specifying bounds on the working scale of the probability distribution, transition probability, and initial distribution parameters. For each matrix, the first column pertains to the lower bound and the second column the upper bound.
#' For data streams, each element of \code{workBounds} should be a k x 2 matrix with the same name of the corresponding element of 
#' \code{Par0}, where k is the number of parameters. For transition probability parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``beta'', where k=\code{length(beta0)}. For initial distribution parameters, the corresponding element of \code{workBounds} must be a k x 2 matrix named ``delta'', where k=\code{length(delta0)}.
#'
#' @return A list of parameter values that can be used as starting values (\code{Par0}) in \code{\link{fitHMM}} or \code{\link{MIfitHMM}}
#' 
#' @seealso \code{\link{getPar}}, \code{\link{getPar0}}, \code{\link{fitHMM}}, \code{\link{MIfitHMM}}
#'
#' @examples
#' # data is a momentuHMMData object, automatically loaded with the package
#' data <- example$m$data
#' stepDist <- "gamma"
#' angleDist <- "vm"
#' nbStates <- 2
#' stepPar0 <- c(15,50,10,20) # natural scale mean_1, mean_2, sd_1, sd_2
#' anglePar0 <- c(0.7,1.5) # natural scale conentration_1, concentration_2
#' 
#' # get working parameters for 'DM' that constrains step length mean_1 < mean_2
#' stepDM <- matrix(c(1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,
#'           dimnames=list(NULL,c("mean:(Intercept)","mean_2",
#'                                "sd_1:(Intercept)","sd_2:(Intercept)")))
#' stepworkBounds <- matrix(c(-Inf,Inf),4,2,byrow=TRUE,
#'                          dimnames=list(colnames(stepDM),c("lower","upper")))
#' stepworkBounds["mean_2","lower"] <- 0 #coefficient for 'mean_2' constrained to be positive
#' wPar0 <- getParDM(nbStates=2,dist=list(step=stepDist),
#'                       Par=list(step=stepPar0),
#'                       DM=list(step=stepDM),workBounds=list(step=stepworkBounds))
#'
#' \dontrun{
#' # Fit HMM using wPar0 as initial values for the step data stream
#' mPar <- fitHMM(data,nbStates=2,dist=list(step=stepDist,angle=angleDist),
#'                Par0=list(step=wPar0$step,angle=anglePar0),
#'                DM=list(step=stepDM),workBounds=list(step=stepworkBounds))
#' }
#' 
#' # get working parameters for 'DM' using 'cov1' and 'cov2' covariates
#' stepDM2 <- list(mean=~cov1,sd=~cov2)
#' wPar20 <- getParDM(data,nbStates=2,dist=list(step=stepDist),
#'                       Par=list(step=stepPar0),
#'                       DM=list(step=stepDM2))
#'
#' \dontrun{
#' # Fit HMM using wPar20 as initial values for the step data stream
#' mPar2 <- fitHMM(data,nbStates=2,dist=list(step=stepDist,angle=angleDist),
#'                Par0=list(step=wPar20$step,angle=anglePar0),
#'                DM=list(step=stepDM2))
#' }
#'
#' @importFrom CircStats circ.mean
# #' @importFrom nleqslv nleqslv
#' @importFrom stats plogis
#' @export
getParDM.default<-function(data=data.frame(),nbStates,dist,
                 Par,
                 zeroInflation=NULL,
                 oneInflation=NULL,
                 estAngleMean=NULL,
                 circularAngleMean=NULL,
                 DM=NULL,userBounds=NULL,workBounds=NULL, ...){
  
  hierArgs <- list(...)
  argNames <- names(hierArgs)[which(names(hierArgs) %in% c("hierStates","hierDist"))]

  ## check that the data is a momentuHMMData object or valid data frame
  if(!is.momentuHMMData(data)){ 
    if(missing(nbStates) & missing(dist)){
      if(all(c("hierStates","hierDist") %in% argNames)){
        return(getParDM.hierarchical(data,hierStates=hierArgs$hierStates,hierDist=hierArgs$hierDist,Par,zeroInflation,oneInflation,estAngleMean,circularAngleMean,DM,userBounds,workBounds))
      }
    }
    if(!is.data.frame(data)) stop('data must be a data.frame')
  }
  if(!missing(nbStates) | !missing(dist)){
    if(any(c("hierStates","hierDist") %in% argNames))
      stop("Either nbStates and dist must be specified (for a regular HMM) or hierStates and hierDist must be specified (for a hierarchical HMM)")
  }
  
  if(nbStates<1) stop('nbStates must be >0')
  
  chkDots(...)
  
  if(!is.list(dist) | is.null(names(dist))) stop("'dist' must be a named list")
  if(!is.list(Par) | is.null(names(Par))) stop("'Par' must be a named list")
  distnames<-names(dist)
  if(!all(distnames %in% names(Par))) stop(distnames[which(!(distnames %in% names(Par)))]," is missing in 'Par'")
  Par <- Par[distnames]
  
  if(is.momentuHMMData(data)){
    if(any(is.na(match(distnames,names(data))))){
      tmpdistnames <- distnames
      for(i in which(is.na(match(distnames,names(data))))){
        if(dist[[distnames[i]]] %in% mvndists){
          if(dist[[distnames[i]]] %in% c("mvnorm2","rw_mvnorm2")){
            tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"))
          } else if(dist[[distnames[i]]] %in% c("mvnorm3","rw_mvnorm3")){
            tmpdistnames <- c(tmpdistnames[-i],paste0(distnames[i],".x"),paste0(distnames[i],".y"),paste0(distnames[i],".z"))          
          }
        }
      }
      if(any(is.na(match(tmpdistnames,names(data))))) stop(paste0(tmpdistnames[is.na(match(tmpdistnames,names(data)))],collapse=", ")," not found in data")
    }
    
    zeroInflation <- vector('list',length(distnames))
    names(zeroInflation) <- distnames
    for(i in distnames){
      if(dist[[i]] %in% zeroInflationdists){
        if(length(which(data[[i]]==0))>0) {
          zeroInflation[[i]]<-TRUE
        }
        else 
          zeroInflation[[i]]<-FALSE
      }
      else zeroInflation[[i]]<-FALSE
    }
    oneInflation <- vector('list',length(distnames))
    names(oneInflation) <- distnames
    for(i in distnames){
      if(dist[[i]] %in% oneInflationdists){
        if(length(which(data[[i]]==1))>0) {
          oneInflation[[i]]<-TRUE
        }
        else 
          oneInflation[[i]]<-FALSE
      }
      else oneInflation[[i]]<-FALSE
    }
  } else {
    if(is.null(zeroInflation)){
      zeroInflation <- vector('list',length(distnames))
      names(zeroInflation) <- distnames
      for(i in distnames){
        zeroInflation[[i]]<-FALSE
      }
    } else {
      if(!is.list(zeroInflation) | is.null(names(zeroInflation))) stop("'zeroInflation' must be a named list")
      for(i in distnames){
        if(is.null(zeroInflation[[i]])) zeroInflation[[i]] <- FALSE
      }
    }
    if(is.null(oneInflation)){
      oneInflation <- vector('list',length(distnames))
      names(oneInflation) <- distnames
      for(i in distnames){
        oneInflation[[i]]<-FALSE
      }
    } else {
      if(!is.list(oneInflation) | is.null(names(oneInflation))) stop("'oneInflation' must be a named list")
      for(i in distnames){
        if(is.null(oneInflation[[i]])) oneInflation[[i]] <- FALSE
      }
    }
    
    if(!all(unlist(lapply(zeroInflation,is.logical)))) stop("zeroInflation must be a list of logical objects")
    if(!all(unlist(lapply(oneInflation,is.logical)))) stop("oneInflation must be a list of logical objects")
    for(i in distnames){
      if(!(dist[[i]] %in% zeroInflationdists) & zeroInflation[[i]])
        stop(dist[[i]]," distribution cannot be zero inflated")
      if(!(dist[[i]] %in% oneInflationdists) & oneInflation[[i]])
        stop(dist[[i]]," distribution cannot be one inflated")
    }
  }
  
  tempCovs <- data[1,,drop=FALSE]
  if(length(data)){
    for(j in names(data)[which(unlist(lapply(data,function(x) any(class(x) %in% meansList))))]){
      if(inherits(data[[j]],"angle")) tempCovs[[j]] <- CircStats::circ.mean(data[[j]][!is.na(data[[j]])])
      else tempCovs[[j]]<-mean(data[[j]],na.rm=TRUE)
    }
  }
  inputs <- checkInputs(nbStates,dist,Par,estAngleMean,circularAngleMean,zeroInflation,oneInflation,DM,userBounds,stateNames=NULL,checkInflation = TRUE)
  
  DMinputs<-getDM(tempCovs,inputs$DM,inputs$dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,zeroInflation,oneInflation,inputs$circularAngleMean,FALSE)
  fullDM <- DMinputs$fullDM
  if(length(data))
    DMind <- getDM(data,inputs$DM,inputs$dist,nbStates,inputs$p$parNames,inputs$p$bounds,Par,zeroInflation,oneInflation,inputs$circularAngleMean,FALSE)$DMind
  else DMind <- DMinputs$DMind
  
  nc <- meanind <- vector('list',length(distnames))
  names(nc) <- names(meanind) <- distnames
  for(i in distnames){
    nc[[i]] <- apply(fullDM[[i]],1:2,function(x) !all(unlist(x)==0))
    # deal with factors
    if(length(tempCovs)){
      for(j in names(which(unlist(lapply(tempCovs,function(x) inherits(x,"factor")))))){
        if(any(grepl(j,inputs$DM[[i]]))){
          tmpCov <- tempCovs
          for(jj in levels(tempCovs[[j]])){
            tmpCov[[j]] <- factor(jj,levels=levels(tempCovs[[j]]))
            tmpgDM<-getDM(tmpCov,inputs$DM[i],inputs$dist[i],nbStates,inputs$p$parNames[i],inputs$p$bounds[i],Par[i],zeroInflation[i],oneInflation[i],inputs$circularAngleMean[i],FALSE)$fullDM[[i]]
            nc[[i]] <- nc[[i]] + apply(tmpgDM,1:2,function(x) !all(unlist(x)==0))
          }
        }
      }
    }
    nc[[i]] <- nc[[i]]>0
    if(!isFALSE(inputs$circularAngleMean[[i]])) {
      meanind[[i]] <- which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
      # deal with angular covariates that are exactly zero
      if(length(meanind[[i]])){
        angInd <- which(is.na(match(gsub("cos","",gsub("sin","",colnames(nc[[i]]))),colnames(nc[[i]]),nomatch=NA)))
        sinInd <- colnames(nc[[i]])[which(grepl("sin",colnames(nc[[i]])[angInd]))]
        nc[[i]][meanind[[i]],sinInd]<-ifelse(nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],sinInd],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)])
        nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)]<-ifelse(nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],gsub("sin","cos",sinInd)],nc[[i]][meanind[[i]],sinInd])
      }
    }
  }
  
  parCount<- lapply(fullDM,ncol)
  for(i in distnames[!unlist(lapply(inputs$circularAngleMean,isFALSE))]){
    parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(fullDM[[i]])))))
  }
  parindex <- c(0,cumsum(unlist(parCount))[-length(fullDM)])
  names(parindex) <- distnames
  
  if(is.null(workBounds)) {
    workBounds <- vector('list',length(distnames))
    names(workBounds) <- distnames
  }
  
  checkNames(distnames,workBounds=workBounds)
  
  wpar <- Par
  for(i in distnames){
    if(!is.null(DM[[i]])){
      if(DMind[[i]] & nrow(unique(fullDM[[i]]))==ncol(unique(fullDM[[i]])) & !inputs$circularAngleMean[[i]]){
        if(length(wpar[[i]])!=nrow(fullDM[[i]])) stop('Par$',i,' should be of length ',nrow(fullDM[[i]]))
        bounds<-inputs$p$bounds[[i]]
        if(any(wpar[[i]]<=bounds[,1] | wpar[[i]]>=bounds[,2])) stop('Par$',i,' must be within parameter bounds')
        bndInd <- which(!duplicated(getboundInd(fullDM[[i]])))
        a<-bounds[bndInd,1]
        b<-bounds[bndInd,2]
        par <- wpar[[i]][bndInd]
        if(any(wpar[[i]]!=par[getboundInd(fullDM[[i]])])) stop('Par$',i,' values are not consistent with DM$',i)
        piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
        ind1<-which(piInd)
        ind2<-which(!piInd)
        
        p<-numeric(length(bndInd))
        if(length(ind1)){
          if(isFALSE(inputs$circularAngleMean[[i]])) p[ind1] <- solve(unique(fullDM[[i]])[ind1,ind1],tan(par[ind1]/2))
          else stop("sorry, circular angle mean parameters are not supported by getParDM")
        }
        
        ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
        ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
        ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
        ind24<-ind2[which(is.infinite(a[ind2]) & is.infinite(b[ind2]))]
        
        if(length(ind21)) p[ind21]<-solve(unique(fullDM[[i]])[ind21,ind21],log(par[ind21]-a[ind21]))
        if(length(ind22)) p[ind22]<-solve(unique(fullDM[[i]])[ind22,ind22],stats::qlogis((par[ind22]-a[ind22])/(b[ind22]-a[ind22])))
        if(length(ind23)) p[ind23]<-solve(unique(fullDM[[i]])[ind23,ind23],-log(-par[ind23]+b[ind23]))
        if(length(ind24)) p[ind24]<-solve(unique(fullDM[[i]])[ind24,ind24],par[ind24])
 
      } else {
        
        if(length(wpar[[i]])!=nrow(fullDM[[i]])) stop('Par$',i,' should be of length ',nrow(fullDM[[i]]))
        bounds<-inputs$p$bounds[[i]]
        if(any(wpar[[i]]<=bounds[,1] | wpar[[i]]>=bounds[,2])) stop('Par$',i,' must be within parameter bounds')
        if(is.list(inputs$DM[[i]])){
          if(isFALSE(inputs$circularAngleMean[[i]]))
            gbInd <- getboundInd(fullDM[[i]])
          else
            gbInd <- c(1:nbStates,getboundInd(fullDM[[i]][1:nbStates+nbStates,])+nbStates)
        } else {
          if(isFALSE(inputs$circularAngleMean[[i]]))
            gbInd <- getboundInd(inputs$DM[[i]])
          else 
            gbInd <- c(1:nbStates,getboundInd(inputs$DM[[i]][1:nbStates+nbStates,])+nbStates)
        }
        bndInd <- which(!duplicated(gbInd))
        a<-bounds[bndInd,1]
        b<-bounds[bndInd,2]
        par <- wpar[[i]][bndInd]
        if(any(wpar[[i]]!=par[gbInd])) stop('Par$',i,' values are not consistent with DM$',i)
        piInd<-(abs(a- -pi)<1.e-6 & abs(b - pi)<1.e-6)
        ind1<-which(piInd)
        ind2<-which(!piInd)
        
        ind21<-ind2[which(is.finite(a[ind2]) & is.infinite(b[ind2]))]
        ind22<-ind2[which(is.finite(a[ind2]) & is.finite(b[ind2]))]
        ind23<-ind2[which(is.infinite(a[ind2]) & is.finite(b[ind2]))]
        ind24<-ind2[which(is.infinite(a[ind2]) & is.infinite(b[ind2]))]
         
        if(length(ind1)){
          if(inputs$estAngleMean[[i]]){
            
            p<-numeric(parCount[[i]])
            
            if(isFALSE(inputs$circularAngleMean[[i]])) {
              
              meanind<-which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],2,function(x) !all(unlist(x)==0))))
              
              asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind1,meanind,drop=FALSE])
              adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
              p[meanind] <- asvd$v %*% adiag %*% t(asvd$u) %*% tan(par[ind1]/2)
              #p[meanind] <- solve(unique(fullDM[[i]])[ind1,meanind],tan(par[ind1]/2))

              if(length(ind21)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind21,-meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[-meanind] <- asvd$v %*% adiag %*% t(asvd$u) %*% log(par[ind21]-a[ind21])
              }
              
              if(length(ind22)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind22,-meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[-meanind] <- asvd$v %*% adiag %*% t(asvd$u) %*% stats::qlogis((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))
              }
              
              if(length(ind23)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind23,-meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[-meanind] <- asvd$v %*% adiag %*% t(asvd$u) %*% (-log(-par[ind23]+b[ind23]))
              }
              
              if(length(ind24)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind24,-meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[-meanind] <- asvd$v %*% adiag %*% t(asvd$u) %*% par[ind24]
              }
            } else {
              
              meanind1<-which((apply(fullDM[[i]][1:nbStates,,drop=FALSE],1,function(x) !all(unlist(x)==0))))
              if(!is.list(inputs$DM[[i]])){
                meanind1<-meanind1[!duplicated(inputs$DM[[i]][meanind1,])]
              }
              meanind2<-which(colSums(nc[[i]][1:nbStates,,drop=FALSE])>0)#which(match(gsub("cos","",gsub("sin","",colnames(fullDM[[i]]))),gsub("cos","",names(which((apply(fullDM[[i]][meanind1,,drop=FALSE],2,function(x) !all(unlist(x)==0)))))),nomatch=0)>0)
              ##if(length(meanind2)) meanind2 <- sort(c(meanind2,meanind2-1))
              #xmat <- fullDM[[i]][gbInd,,drop=FALSE][meanind1,meanind2,drop=FALSE]
              ##nc<-apply(xmat,1:2,function(x) !all(unlist(x)==0))
              
              solveatan2<-function(x,theta,covs,nbStates,nc,meanind,oparms,circularAngleMean){
                XB<-getXB(covs,1,c(x,rep(1,oparms)),TRUE,circularAngleMean,FALSE,nbStates,nc,meanind)[meanind,]
                c(abs(theta - XB),rep(0,max(0,length(x)-length(theta))))
              }

              if(length(meanind1)) {
                if (!requireNamespace("nleqslv", quietly = TRUE)) {
                  stop("Package \"nleqslv\" needed for this function to work. Please install it.",
                       call. = FALSE)
                }
                p[1:(length(meanind2)/2)] <- nleqslv::nleqslv(x=rep(1,length(meanind2)/2),fn=solveatan2,theta=par[meanind1],covs=fullDM[[i]],nbStates=nbStates,nc=nc[[i]],meanind=meanind1,oparms=parCount[[i]]-length(meanind2)/2,circularAngleMean=inputs$circularAngleMean[[i]],control=list(allowSingular=TRUE))$x[1:(length(meanind2)/2)]
              }
              meanind<-which((apply(fullDM[[i]][nbStates+1:nbStates,,drop=FALSE],2,function(x) !all(unlist(x)==0))))
              
              if(length(ind21)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind21,meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[meanind-length(meanind2)/2] <- asvd$v %*% adiag %*% t(asvd$u) %*% log(par[ind21]-a[ind21])
              }
              
              if(length(ind22)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind22,meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[meanind-length(meanind2)/2] <- asvd$v %*% adiag %*% t(asvd$u) %*% stats::qlogis((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))
              }
              
              if(length(ind23)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind23,meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[meanind-length(meanind2)/2] <- asvd$v %*% adiag %*% t(asvd$u) %*% (-log(-par[ind23]+b[ind23]))
              }
              
              if(length(ind24)){
                asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind24,meanind,drop=FALSE])
                adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
                p[meanind-length(meanind2)/2] <- asvd$v %*% adiag %*% t(asvd$u) %*% par[ind24]
              }
            
            }
          } else if(!length(ind2)){
            p <- solve(fullDM[[i]][gbInd,],tan(par/2))
          } else stop("sorry, the parameters for ",i," cannot have different bounds")
        } else if(((length(ind21)>0) + (length(ind22)>0) + (length(ind23)>0) + (length(ind24)>0))>1){ 
          stop("sorry, getParDM requires the parameters for ",i," to have identical bounds when covariates are included in the design matrix")
        } else {
          if(length(ind21)){
            asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind21,,drop=FALSE])
            adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
            p <- asvd$v %*% adiag %*% t(asvd$u) %*% log(par[ind21]-a[ind21])
          }
          
          if(length(ind22)){
            asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind22,,drop=FALSE])
            adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
            p <- asvd$v %*% adiag %*% t(asvd$u) %*% stats::qlogis((par[ind22]-a[ind22])/(b[ind22]-a[ind22]))
          }
          
          if(length(ind23)){
            asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind23,,drop=FALSE])
            adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
            p <- asvd$v %*% adiag %*% t(asvd$u) %*% (-log(-par[ind23]+b[ind23]))
          }

          if(length(ind24)){
            asvd<-svd(fullDM[[i]][gbInd,,drop=FALSE][ind24,,drop=FALSE])
            adiag <- diag(x=1/asvd$d,nrow=ifelse(length(asvd$d)>1,length(asvd$d),1),ncol=ifelse(length(asvd$d)>1,length(asvd$d),1))
            p <- asvd$v %*% adiag %*% t(asvd$u) %*% par[ind24]
          }
        }
      }
      if(any(!is.finite(p))) stop(i," working scale parameters are not finite. Check natural parameter values, bounds, and constraints.")
      
      workBounds[i] <- getWorkBounds(workBounds[i],i,p,parindex[i]-parindex[[i]],parCount,inputs$DM)[i]
      if(any(p<workBounds[[i]][,1]) | any(p>workBounds[[i]][,2])) stop("could not find valid working scale parameters for ",i," that satisfy workBounds")
      
      p <- nw2w(p,workBounds[[i]])
      
      if(any(!is.finite(p))) stop("could not find valid working scale parameters for ",i," that satisfy workBounds")
      
      wpar[[i]]<-c(p)      
    }
  }
  wpar
}

#' @rdname getParDM
#' @method getParDM hierarchical
#' @param hierStates A hierarchical model structure \code{\link[data.tree]{Node}} for the states.  See \code{\link{fitHMM}}.
#' @param hierDist A hierarchical data structure \code{\link[data.tree]{Node}} for the data streams. See \code{\link{fitHMM}}.
#'
#' @export
getParDM.hierarchical<-function(data=data.frame(), hierStates, hierDist,
                       Par,
                       zeroInflation=NULL,
                       oneInflation=NULL,
                       estAngleMean=NULL,
                       circularAngleMean=NULL,
                       DM=NULL, userBounds=NULL, workBounds=NULL, ...){
  
  ## check that the data is a momentuHierHMMData object or valid data frame
  if(!is.momentuHierHMMData(data)){
    if(!is.data.frame(data)) stop('data must be a data.frame')
  } else {
    data <- momentuHMMData(data)
  }
  
  chkDots(...)
  
  installDataTree()
  
  nbStates <- nbHierStates(hierStates)$nbStates
  dist <- getHierDist(hierDist,data=NULL,checkData=FALSE)
  
  return(getParDM.default(data,nbStates,dist,Par,zeroInflation,oneInflation,estAngleMean,circularAngleMean,DM,userBounds,workBounds))
}
