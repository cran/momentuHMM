
#' Print \code{miSum}
#' @method print miSum
#'
#' @param x A \code{miSum} object.
#' @param ... Currently unused. For compatibility with generic method.
#'
#' @examples
#' \dontshow{
#' set.seed(3,kind="Mersenne-Twister",normal.kind="Inversion")
#' }
#' \dontrun{
#' # Extract data from miExample
#' obsData <- miExample$obsData
#' 
#' # error ellipse model
#' err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)
#' 
#' # Fit crawl to obsData
#' crwOut <- crawlWrap(obsData,theta=c(4,0),fixPar=c(1,1,NA,NA),
#'                     err.model=err.model)
#'                     
#' # Fit four imputations
#' bPar <- miExample$bPar
#' HMMfits <- MIfitHMM(crwOut,nSims=4,poolEstimates=FALSE,
#'                    nbStates=2,dist=list(step="gamma",angle="vm"),
#'                    Par0=bPar$Par,beta0=bPar$beta,
#'                    formula=~cov1+cos(cov2),
#'                    estAngleMean=list(angle=TRUE),
#'                    covNames=c("cov1","cov2"))
#'                    
#' # Pool estimates
#' miSum <- MIpool(HMMfits)
#' print(miSum)
#' }
#' @export

print.miSum <- function(x,...)
{
  m <- x # the name "x" is for compatibility with the generic method
  m <- delta_bc(m)
  
  m <- formatmiSum(m)
  m$CIbeta <- m$Par$beta
  m$CIreal <- m$Par$real
  
  if(inherits(m,"hierarchical")) print.momentuHierHMM(m)
  else print.momentuHMM(m)
}
