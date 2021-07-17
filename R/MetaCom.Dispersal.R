#' Seed dispersion in a metacommunity
#'
#' Seeds can disperse to any other local community with the same probability;
#' i.e. probability to disperse other subcommunity/(number local communities - 1). Each seed is dispersed independently.
#'
#'Both input and output is a matrix where seeds are in the rows, and their
#'attributes (i.e. location, species identity and trait values) are in the
#'columns.
#'
#'
#'@param before A matrix of seed's attributes; seeds in rows, their location,
#'species identity and traits are in columns. Column that contains information on
#' locality has to be called 'site'
#'@param m probability that a seed are dispersed into other (sub)community
#'@param n number of local (sub)communities
#'@param ... Additional parameters. It necessary for thechnical reasons: the framework don't know the current list of parameters when call this function
#'@return Same type of matrix as \code{before}
#'@seealso \code{\link{fDispersal}}
  MetaCom.Dispersal<-function(n,before,m=0.1,...)
{
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$Dispersal.params))
  parameters$Dispersal.params<-list(m=m)
  assign("parameters",parameters,envir = comsimitvEnv)
  after<-before
  for (u in 1:nrow(before))
  {
    if (sum(after[,"site"]==after[u,"site"])>1)
      {
      prob<-rep(m/(n-1),n)
      prob[before[u,"site"]]<-(1-m)
      after[u,"site"]<-sample(1:n,size=1,prob=prob)
      }
  }
  return(after)
}

