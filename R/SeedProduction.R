#' Calculating number of produced seeds
#'
#' Number of seeds calculated following the formula used by Botta-Dukat & Czucz (2016).
#' This built-in function can be replaced by a user-defined one.
#'
#' Each individual produces one seed or does not produce seed at all.
#' Probability of seed production is a decreasing sigmoid function of strength of
#' competition  (sum of abundances weighted by competition coefficients).
#' If strength of competition is higher than parameter K, probability is set to
#' zero. See \code{vignette("competition")} for formulas
#'
#'@param compet Matrix of pairwise competition coefficients
#'@param b0 Probability of producing seed, if no competition
#'@param K Critical level of competition (See Details)
#'@param ... any additional parameters
#'@return Matrix of produced seeds
SeedProduction<-function(compet,b0=1,K=200,...)
{
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$SeedProduction.params))
    parameters$SeedProduction.params<-list(b0=b0,
                                         K=K)
  assign("parameters",parameters,envir = comsimitvEnv)
  seed<-rep(0,nrow(compet))
  NE<-colSums(compet)
  birth.limit<-b0*(K-NE)/K
  birth.limit[birth.limit<0]<-0
  seed<-stats::rbinom(nrow(compet),size=1,birth.limit)
  if (sum(seed)==0) seed[which.max(birth.limit)]<-1
  return(seed)
}

