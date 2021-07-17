#' Asymmetric competition kernels
#'
#' It calculates asymmetric competition coefficients
#'
#'
#'Depending on value of ac.type the convex-concave function from Kisdi (1999)
#'or  smooth function suggested by Nattrass et al (2012) are used.
#'
#'For formulas and meaning of parameters see the \code{vignette("competition")}
#'
#'@param trait.values Dataframe of all traits 
#'@param trait.compet Name of trait related to resource use
#'@param ac.type Type of the function (see \code{vignette("competition")})
#'@param sigma.b steepness of competition kernel
#'@param ac.C parameter influencing shape of the function (has to be positive)
#'@param ac.v parameter influencing shape of the function (has to be positive)
#'@param ... Any additional parameters
#'@references	Kisdi, E. (1999) Evolutionary Branching under Asymmetric
#'Competition
#'\emph{Journal of Theoretical Biology} \bold{197}(2): 149-162.
#'\doi{10.1006/jtbi.1998.0864}
#'@references Nattrass, S., Baigent, S., & Murrell, D. J. (2012) Quantifying
#'the Likelihood of Co-existence for Communities with Asymmetric Competition.
#'\emph{Bulletin of Mathematical Biology}, \bold{74}(10): 2315–2338.
#'\doi{10.1007/s11538-012-9755-8}
#'@seealso \code{\link{competition.kernel}}
#'
asymmetric.competition.kernel<-function(trait.values, trait.compet="trait.b",ac.type=c("Kisdi","Nattrass"),
                                      sigma.b=0.03,ac.C=1,ac.v=1,...)
{
  ac.type<-match.arg(ac.type)
  if (ac.C<=0) stop("Parameter ac.C has to be positive")
  if (ac.v<=0) stop("Parameter ac.C has to be positive")
  parameters<-get("parameters",envir = comsimitvEnv)
  if (is.null(parameters$competition.kernel.params))
    parameters$competition.kernel.params<-list(trait.compet=trait.compet,sigma.b=sigma.b, ac.type=ac.type,
                                               ac.C=ac.C,ac.v=ac.v)
  assign("parameters",parameters,envir = comsimitvEnv)
  trait.values<-trait.values[,trait.compet]
  S=length(trait.values)
  compet <- matrix(0,S,S)
  for (i in 1:(S-1))
    for (j in i:S)
    {
    d<-trait.values[i]-trait.values[j]
    if (ac.type=="Kisdi")
      {
      if (sigma.b==0) stop("Parameter sigma.b shouldn't be zero")
      if (sigma.b>0)
        {
        compet[i,j]<-ac.C*(1-(1/(1+ac.v*exp(-d/sigma.b))))
        compet[j,i]<-ac.C*(1-(1/(1+ac.v*exp(d/sigma.b))))
        }
      else
       {
        if (d<0)
          {
          compet[i,j]<-ac.C
          compet[j,i]<-0
          }
        if (d>0)
         {
           compet[i,j]<-0
           compet[j,i]<-ac.C
         }
        if (d==0)
         {
           compet[i,j]<-ac.C*ac.v/(ac.v+1)
           compet[j,i]<-ac.C*ac.v/(ac.v+1)
         }
        }
      }
    if (ac.type=="Nattrass")
      {
      if (sigma.b>0)
        {
        compet[i,j]<-1+ac.C-(2*ac.C/(1+exp(-2*d/sigma.b)))
        compet[j,i]<-1+ac.C-(2*ac.C/(1+exp(-2*d/sigma.b)))
        }
      else
       {
        if (d<0)
          {
          compet[i,j]<-1+ac.C
          compet[j,i]<-1-ac.C
          }
        if (d>0)
         {
           compet[i,j]<-1-ac.C
           compet[j,i]<-1+ac.C
         }
        if (d==0)
         {
           compet[i,j]<-1
           compet[j,i]<-1
         }
        }
      }

    }

  return(compet)
}


