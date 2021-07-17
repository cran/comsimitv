#' Intraspecific Trait Variation
#'
#' This function adds a random noise to mother's trait values of each seed
#'
#' The function uses parameters of \code{\link{Gener.species.pool}}. First it
#' transforms back mother's trait values to multivariate normal distribution.
#' Then random noise was added to this values. Random noise has multivariate
#' normal distribution, with zero means and the same  \strong{correlation} structure
#' as specified in parameter \emph{sigma}. \strong{Note} that \emph{sigma} specifies
#' covariance matrix, not correlation structure \emph{per se}. Variances in the random noise
#' are diagonals (i.e. variance componens) of parameter \emph{sigma} multiplied by
#' \emph{ITV.ratio}. The non-diagonal elements of covariance matrix were specified
#'  to conserve the correlation structure among traits.
#'
#'
#'@param seeds Matrix of produced seeds (with mother'trait values) as produced by
#'             \code{\link{SeedProduction}} function
#'@param n.traits Number of traits
#'@param distribs Types of the distributions of traits (see \code{\link{Gener.species.pool}})
#'@param distr.parms Parameters of distribution (see \code{\link{Gener.species.pool}})
#'@param sigma Matrix of variance-covariance matrix of traits (see \code{\link{Gener.species.pool}})
#'@param ITV.ratio Ratio of within/between species variances of traits
#'@param ... Any additional parameters
#'@return Matrix of produced seeds as produced by
#'             \code{\link{SeedProduction}} function

randomITV<-function(seeds=matrix(),
                    n.traits=3,
                    distribs=rep("unif",n.traits),
                    distr.parms=list(),
                    sigma=diag(1,n.traits,n.traits),ITV.ratio=0.01,...)
{
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$ITV.params))
    parameters$ITV.params<-list(ITV.ratio=ITV.ratio)
  assign("parameters",parameters,envir = comsimitvEnv)
  
  if (!length(distribs)==n.traits) distribs<-rep(distribs[1],n.traits)
  if (length(distr.parms)==1) 
    {
    w<-vector("list",n.traits)
    for (i in 1:n.traits)  w[[i]]<-dist.parms
    dist.parms<-w
    }

  
  
  f<-paste("p",distribs,sep="")
  for (i in 1:n.traits) 
    {
    if (length(distr.parms)>0) seeds[,i+2]<-do.call(f[i],c(list(q=seeds[,i+2]),distr.parms[[i]]))
    else seeds[,i+2]<-do.call(f[i],list(q=seeds[,i+2]))
    seeds[,i+2]<-stats::qnorm(seeds[,i+2])
    }
  if (mean(sigma==diag(1,n.traits,n.traits))==1)
    for (i in 1:n.traits)  seeds[,i+2]<-seeds[,i+2]+stats::rnorm(nrow(seeds),0,ITV.ratio)
  else
    {
    x<-MASS::mvrnorm(n=nrow(seeds),rep(0,n.traits),sigma*ITV.ratio)
    for (i in 1:n.traits) seeds[,i+2]<-seeds[,i+2]+x[,i]
    }
  
  f<-paste("q",distribs,sep="")
  
  for (i in 1:n.traits) 
    {
    seeds[,i+2]<-stats::pnorm(seeds[,i+2])
    if (length(distr.parms)>0) seeds[,i+2]<-do.call(f[i],c(list(p=seeds[,i+2]),distr.parms[[i]]))
    else seeds[,i+2]<-do.call(f[i],list(p=seeds[,i+2]))
    }

  return(seeds)
}
