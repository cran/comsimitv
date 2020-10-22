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
#'@param distribs Types of the distributions of traits (see \code{\link{Gener.species.pool}})
#'@param distr.parms Parameters of distribution (see \code{\link{Gener.species.pool}})
#'@param sigma Matrix of variance-covariance matrix of traits (see \code{\link{Gener.species.pool}})
#'@param ITV.ratio Ratio of within/between species variances of traits
#'@param ... Any additional parameters
#'@return Matrix of produced seeds as produced by
#'             \code{\link{SeedProduction}} function

randomITV<-function(seeds=matrix(),
              distribs=rep("unif",3),
              distr.parms=list(a=list(min=0,max=1),b=list(min=0,max=1),c=list(min=0,max=1)),
              sigma=diag(1,3,3),ITV.ratio=0.01,...)
{
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$ITV.params))
    parameters$ITV.params<-list(ITV.ratio=ITV.ratio)
  assign("parameters",parameters,envir = comsimitvEnv)
  if (!length(distribs)==3) distribs<-rep(distribs[1],3)
  if (length(distr.parms)==1) distr.parms<-list(a=distr.parms[[1]],b=distr.parms[[1]],c=distr.parms[[1]])

  f<-paste("p",distribs,sep="")
  seeds[,"trait.a"]<-do.call(f[1],c(list(q=seeds[,"trait.a"]),distr.parms$a))
  seeds[,"trait.b"]<-do.call(f[2],c(list(q=seeds[,"trait.b"]),distr.parms$b))
  seeds[,"trait.c"]<-do.call(f[3],c(list(q=seeds[,"trait.c"]),distr.parms$c))

  seeds[,"trait.a"]<-stats::qnorm(seeds[,"trait.a"])
  seeds[,"trait.b"]<-stats::qnorm(seeds[,"trait.b"])
  seeds[,"trait.c"]<-stats::qnorm(seeds[,"trait.c"])

  if (sum(sigma==diag(1,3,3))==9)
  {
    seeds[,"trait.a"]<-seeds[,"trait.a"]+stats::rnorm(nrow(seeds),0,ITV.ratio)
    seeds[,"trait.b"]<-seeds[,"trait.b"]+stats::rnorm(nrow(seeds),0,ITV.ratio)
    seeds[,"trait.c"]<-seeds[,"trait.c"]+stats::rnorm(nrow(seeds),0,ITV.ratio)
  }
  else
  {

    x<-MASS::mvrnorm(n=nrow(seeds),rep(0,3),sigma*ITV.ratio)
    seeds[,"trait.a"]<-seeds[,"trait.a"]+x[,1]
    seeds[,"trait.b"]<-seeds[,"trait.b"]+x[,2]
    seeds[,"trait.c"]<-seeds[,"trait.c"]+x[,3]
  }

  seeds[,"trait.a"]<-stats::pnorm(seeds[,"trait.a"])
  seeds[,"trait.b"]<-stats::pnorm(seeds[,"trait.b"])
  seeds[,"trait.c"]<-stats::pnorm(seeds[,"trait.c"])

  f<-paste("q",distribs,sep="")
  seeds[,"trait.a"]<-do.call(f[1],c(list(p=seeds[,"trait.a"]),distr.parms$a))
  seeds[,"trait.b"]<-do.call(f[2],c(list(p=seeds[,"trait.b"]),distr.parms$b))
  seeds[,"trait.c"]<-do.call(f[3],c(list(p=seeds[,"trait.c"]),distr.parms$c))


  return(seeds)
}
