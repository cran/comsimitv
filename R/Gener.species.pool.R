#' Generating trait values for the species pool
#'
#' It generates random trait values for species. Each species (individual)
#' are characterized by three traits.
#'
#' Each species are characterized by three traits called trait A,
#'  B and C. Trait A describes the habitat preference,
#'  trait B influences the competitive interactions, while trait C
#'  is a completely neutral trait.
#'
#' Any standard distribution of stats package can be used for generating
#' the random numbers. For list of these distribution see \code{\link[stats]{Distributions}}
#' In stats package the functions for the density/mass function are named
#' in the form dxxx."xxx" (without d!) as string (i.e. between quatation marks)
#'  should be supplied for parameter \code{distribs}.
#'
#'  In this step single value of each trait is generated for each species,
#'  i.e. there is no intraspecific trait variation.
#'
#'  If traits are independent (it is the default option), random number
#'  generating functions are called with parameters specified by the user.
#'
#'  Otherwise, a variance-covariance matrix has to be given. First, triplets of
#'  random numbers are drawn from multivariate normal distribution with zero
#'  means and the supplied variance-covariance matrix as parameters.
#'  Then these random numbers are converted to probability by standard normal
#'  probability function, and then these probabilities converted to trait values
#'   using quantile function of selected distribution with parameters given
#'   by the user.
#'
#'@param S Species pool size
#'@param distribs Types of the distributions of traits
#'@param distr.parms Parameters of distribution (see Details)
#'@param sigma Matrix of variance-covariance matrix of traits
#'@param ... Any additional parameters
#'@return A list with three elemets:
#'@return               a = vector of fist trait,
#'@return               b = vector of second trait,
#'@return               c = vector of third trait
Gener.species.pool<-function(S,distribs=rep("unif",3),
                             distr.parms=list(a=list(min=0,max=1),b=list(min=0,max=1),c=list(min=0,max=1)),
                             sigma=diag(1,3,3),...)
{
  if (!length(distribs)==3) distribs<-rep(distribs[1],3)
  if (length(distr.parms)==1) distr.parms<-list(a=distr.parms[[1]],b=distr.parms[[1]],c=distr.parms[[1]])
  if (sum(sigma==diag(1,3,3))==9)
  {
    f<-paste("r",distribs,sep="")
    trait.a<-do.call(f[1],c(list(n=S),distr.parms$a))
    trait.b<-do.call(f[2],c(list(n=S),distr.parms$b))
    trait.c<-do.call(f[3],c(list(n=S),distr.parms$c))
  }
  else
  {
    x<-MASS::mvrnorm(n=S,rep(0,3),sigma)
    f<-paste("q",distribs,sep="")
    trait.a<-do.call(f[1],c(list(p=stats::pnorm(x[,1])),distr.parms$a))
    trait.b<-do.call(f[2],c(list(p=stats::pnorm(x[,2])),distr.parms$b))
    trait.c<-do.call(f[3],c(list(p=stats::pnorm(x[,3])),distr.parms$c))
  }
  traits<-list(a=trait.a,
               b=trait.b,
               c=trait.c)
  return(traits)
}

