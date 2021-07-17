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
#'@param n.traits Number of traits
#'@param distribs Types of the distributions of traits
#'@param distr.parms Parameters of distribution (see Details)
#'@param sigma Matrix of variance-covariance matrix of traits
#'@param ... Any additional parameters
#'@return A data frame with traits as columns


Gener.species.pool<-function(S,n.traits=3,distribs=rep("unif",n.traits),
                             distr.parms=list(),
                             sigma=diag(1,n.traits,n.traits),...)
{
  if (length(distribs)==1) distribs<-rep(distribs,n.traits)
  if (!(length(distribs)==n.traits)) stop("Number of distributions should be equal to number of traits!")
  if (length(distribs)>0)
    {
    if (length(distr.parms)==1) 
      {
      w<-vector("list",n.traits)
      for (i in 1:n.traits)  w[[i]]<-dist.parms
      dist.parms<-w
      }
    if (length(distr.parms)<n.traits)
      {
      w<-vector("list",n.traits)
      w[[1:length(distr.parms)]]<-distr.parms
      distr.parms<-w
      warning("Less parameter set than traits! Default parameters will be used.")
      }
    }
  
  
  traits<-matrix(NA,nrow=S,ncol=n.traits)
  colnames(traits)<-paste("trait",letters[1:n.traits],sep=".")
  traits<-as.data.frame(traits)
  if (mean(sigma==diag(1,n.traits,n.traits))==1)
  {
    f<-paste("r",distribs,sep="")
    for (i in 1:n.traits) traits[,i]<-do.call(f[i],c(list(n=S),distr.parms[[i]]))
  }
  else
  {
    x<-MASS::mvrnorm(n=S,rep(0,n.traits),sigma)
    f<-paste("q",distribs,sep="")
    traits[,i]<-do.call(f[i],c(list(p=stats::pnorm(x[,i])),distr.parms[[i]]))
  }
  
return(traits)
}

