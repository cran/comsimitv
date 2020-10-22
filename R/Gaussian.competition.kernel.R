#' Gaussian competition kernel
#'
#' It calculates pairwise competition coefficients as overlap of Gaussian resource utilization curve
#'
#' It assumes that each species has Gaussian resource utilization curve:
#' \deqn{\exp(\frac{(x-trait.value)^2}{sigma.b})}{exp([(x-trait.value)^2]/(sigma.b/2))}
#' where:      x = quality of resource (e.g. seed size or rooting depth\cr
#'
#' Optima of curves depend on trait value related to resource use,
#' while standard deviation is the same for all species (note that for technical reason
#' parameter \code{sigma.b} is twice of the common sqared s.d.).
#' Pairwise competition coefficients are calculated as overlap of
#' resource utilization functions (MacArthur & Levins 1967).See details in
#' \code{vignette("competition")}
#'
#'@param trait.values Values of trait related to resource use
#'@param sigma.b Width of Gaussian kernel
#'@param ... Any additional parameters
#'
#'@references	MacArthur R, Levins R (1967) The Limiting Similarity,
#'Convergence, and Divergence of Coexisting Species.
#'\emph{The American Naturalist} \bold{101}: 377-385.
#'\url{http://dx.doi.org/10.1086/282505}
#'@seealso \code{\link{competition.kernel}}
Gaussian.competition.kernel<-function(trait.values,sigma.b=0.03,...)
{
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$competition.kernel.params))
    parameters$competition.kernel.params<-list(sigma.b=sigma.b)
  assign("parameters",parameters,envir = comsimitvEnv)
  S=length(trait.values)
  compet <- matrix(0,S,S)
  if (sigma.b==0) diag(compet) <- 1
  if (sigma.b==Inf) compet <- matrix(1,S,S)
  if ((sigma.b>0) & (sigma.b<Inf))
  {
    trait.dist <- as.matrix(stats::dist(trait.values))
    compet <- exp(-trait.dist^2/sigma.b)
  }
  return(compet)
}
