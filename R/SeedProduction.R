#' Calculating number of produced seeds
#'
#' Number of seeds calculated following the formula used by Botta-Dukat & Czucz (2016).
#' This built-in function can be replaced by a user-defined one.
#'
#' Expected value of produced seeds is a decreasing sigmoid function of strength of
#' competition  (sum of abundances weighted by competition coefficients).
#' If strength of competition is higher than parameter K, probability is set to
#' zero. See \code{vignette("competition")} for formulas
#'
#' In simulation of Botta-Dukat & Czucz (2016) each individual produces one seed or does 
#' not produce seed at all. In this case number of seeds follows binomial distribution 
#' (i.e. distrib="binom"). A more realistic alternative is using Poisson distribution 
#' (distrib="pois").
#'
#'@param compet Matrix of pairwise competition coefficients
#'@param b0 Probability of producing seed, if no competition
#'@param K Critical level of competition (See Details)
#'@param seed.distrib Distribution of seed numbers (See Details)
#'@param ... any additional parameters
#'@return Matrix of produced seeds
SeedProduction<-function(compet,b0=1,K=200,seed.distrib=c("pois","binom"),...)
{
  seed.distrib <- match.arg(seed.distrib)
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$SeedProduction.params))
    parameters$SeedProduction.params<-list(b0=b0,
                                           K=K,
                                           seed.distrib=seed.distrib)
  assign("parameters",parameters,envir = comsimitvEnv)
  seed<-rep(0,nrow(compet))
  NE<-colSums(compet)
  birth.limit<-b0*(K-NE)/K
  birth.limit[birth.limit<0]<-0
  if (seed.distrib=="binom")  seed<-stats::rbinom(n=nrow(compet),size=1,prob=birth.limit)
  else seed<-stats::rpois(n=nrow(compet),lambda=birth.limit)
  if (sum(seed)==0) seed[which.max(birth.limit)]<-1
  return(seed)
}

