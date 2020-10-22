#' Bell-shaped tolerance function
#'
#' It calculates probability of seedling's survival from their trait related to habitat filtering and the local environment.
#'
#' It assumes that probability of seedling's survival is maximal if the local
#' environment has the same value as its trait. Survival probability decrease as
#' environmental value departs from the optimum according to a Gaussian (bell-shaped)
#' curve. The speed of decrease depends on the tolerance width parameter (\code{sigma.a}).
#'
#'@param trait.values Values of trait related to habitat preference
#'@param env Vector of environmental conditions in the local communities
#'@param sigma.a Tolerance width (same for all species)
#'@param ... Any additional parameters
#'@return A matrix of survival probabilities, communities in rows, species/individuals in columns
#'@seealso \code{\link{tolerance}}
Gaussian.tolerance<-function(trait.values,env,sigma.a=0.001,...)
{
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$Survive.params))
    parameters$Survive.params<-list(sigma.a=sigma.a)
  assign("parameters",parameters,envir = comsimitvEnv)
  S<-length(trait.values)
  n<-length(env)
  X<-matrix(rep(env,S),ncol=S)
  A<-t(matrix(rep(trait.values,n),ncol=n))
  survive <- if (sigma.a<Inf) exp(-((X-A)^2)/sigma.a) else matrix(1,nrow(X),ncol(X))
  for (i in 1:n)
    if (max(survive[i,])<=0.01) survive[i,which.max(survive[i,])]<-1
  survive<-pmax(survive-0.01,0)
  return(survive)
}

