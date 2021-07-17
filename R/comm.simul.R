#' Framework for community assembly simulation
#'
#' Flexible framework of individual-based simulation of community assembly
#' following framework proposed by Botta-Dukat & Czucz (2016), but allowing
#' intraspecific trait variation (ITV)
#'
#'This function is a framework for simulation of assembly in a meta-community.
#'The simulation consists of a community initialization followed by an
#'iterative simulation of a "disturbance–regeneration" cycle.
#'During initialization a species pool is created defining each species by
#'its trait values. Each locality is characterized by an environmental variable.
#'Initial composition of local communities is a random selection from the species
#'pool: species identity is selected independently for each individal with
#'probability of seedling survival (that depends on local environment and trait value).
#'
#'
#'The "disturbance-regeneration" cycle consists of the following steps:
#'\enumerate{
#'       \item disturbance event: some randomly selected individuals
#'             die in each community
#'       \item survivors produce seeds. Seed production depends on fertility
#'             of the locality and competition among coexisting individuals
#'       \item seeds are dispersed among localities
#'       \item all seeds germinate and seedlings struggle for survival. The number of
#'             adults in local communities is fixed, thus number of seedlings that can
#'             survive and grow up equals to the number of individulas died in the
#'             disturbance event (in the recent version one individual dies, but planed development
#'             is introducing a disturbance severity/number of deaths parameter)
#'}
#'
#'
#' It is a flexible framework that calls funcions for:
#'
#' \itemize{
#'    \item generating species pool
#'          (\code{\link{Gener.species.pool}})
#'    \item calculating pairwise competition coefficients
#'          (\code{\link{competition.kernel}})
#'    \item calculating seedling's survival probabilities
#'          (\code{\link{tolerance}})
#'    \item calculating number of produced seeds
#'          (\code{\link{SeedProduction}})
#'    \item calculating trait values of offsprings
#'          (\code{\link{fITV}})
#'    \item seed dispersal among localities
#'          (\code{\link{fDispersal}})
#'         }
#'Functions available in the package can be easily replaced by user-defined
#'functions.
#'
#'
#'@param x Vector of environmental values in communities. If not given, 40
#'         communities are created, with environmental variable equally
#'         spacing from 0.11 to 0.89
#'@param S Species pool size
#'@param n.traits Number of traits
#'@param J Number of individuals in each community
#'@param rand.seed
#'       Random seed number. Setting the same value allows repeating
#'       the same simulation
#'@param sim.length
#'       Length of simulation. \code{sim.length*S} cycle (disturbance-seed
#'       production-dispersal-establishment) will happen.
#'@param fSpecPool
#'       Name of (the user defined) function that generates the
#'       species pool. See \code{\link{Gener.species.pool}}
#'@param fSurvive
#'       Name of the (user defined) function for calculating survival
#'       probability of seeds. See more details
#'       in available functions and specification of your own
#'       function in \code{\link{tolerance}}
#'@param competition.kernel
#'       Name of the (user defined) function for calculating
#'       pairwise competition coefficients. See more details
#'       in available functions and specification of your own
#'       function in \code{\link{competition.kernel}}
#'@param fSeedProduction
#'       Name of the user defined function for calculating
#'       number of produced seeds See \code{\link{SeedProduction}}
#'@param fDispersal
#'       Name of the user defined function for dispersal of
#'       produced seeds among local communities.
#'       See more details in available functions and specification
#'       of your own function in \code{\link{fDispersal}}
#'@param fITV
#'       Name of the function that define seeds trait values, possibly
#'       considering mother's trait and mothers environment.
#'       If "noITV", there is no intraspecific trait variation.
#'       See more details in available functions and specification
#'       of your own function in \code{\link{fITV}}
#'@param verbose
#'       Runing may take long time. If \code{verbose} set to \code{TRUE},
#'       it writes messages into the screen indicating the progress.
#'@param ... Additional parameters of functions called by the framework.
#'@return A list with two elements:
#'@return \code{$final.community} a dataframe containing data on individuals in the final meta-community.
#'                         Each individual represented by a row; columns are: sub-community,
#'                         species identity, trait values.
#'@return \code{$parameters} list of simulation parameters (including parameters of functions called by the framework function)
#'
#'@references	Botta-Dukat Z, Czucz B (2016) Testing the ability of
#'functional diversity indices to detect trait convergence and divergence
#'using individual-based simulation.
#'\emph{Methods in Ecology and Evolution} \bold{7}(1): 114-126.
#'\doi{10.1111/2041-210X.12450}
#'
#'@export
#'
#'@examples
#' w<-comm.simul(S=20, J=30)
#' str(w)
#'
#' set.seed(1)
#' w<-comm.simul(S=20, J=30, fITV=NULL)$final.community
#' w[w[,2]==1,] # Each individuals belonging to Species1 has the same trait values
comm.simul<-function(x=vector(),S=200, n.traits=3, J=300,rand.seed=NULL, sim.length=1,
                         fSpecPool="Gener.species.pool",
                         competition.kernel="Gaussian.competition.kernel",
                         fSurvive="Gaussian.tolerance",
                         fSeedProduction="SeedProduction",
                         fDispersal="MetaCom.Dispersal",
                         fITV="randomITV",
                         verbose=FALSE,
                         ...)
{
  if (length(x)==0) x <- seq(0.1,0.9,0.8/49)
  n<-length(x)

  n.trait<-3
  
  parameters<-list(x=x,
                   s=S,
                   J=J,
                   sim.length=sim.length,
                   rand.seed=rand.seed,
                   fSpecPool=fSpecPool,
                   fSpecPool.params=NULL,
                   competition.kernel=competition.kernel,
                   competition.kernel.params=NULL,
                   fSurvive=fSurvive,
                   Survive.params=NULL,
                   fSeedProduction=fSeedProduction,
                   SeedProduction.params=NULL,
                   fDispersal=fDispersal,
                   Dispersal.params=NULL,
                   fITV=fITV,
                   ITV.params=NULL)
  assign("parameters",parameters,envir = comsimitvEnv)
  set.seed(rand.seed)

  if (verbose) cat("Generating species pool... \n")
  traits.modus<-do.call(fSpecPool,c(list(S=S,n.traits=n.traits),list(...)))

  if (verbose) cat("Generating starting community composition...\n")
  survive <-do.call(fSurvive,c(list(trait.values=traits.modus,env=x),list(...)))

  Y<-array(NA,dim=c(n,J,n.traits+1)) # species abundances
  dimnames(Y)[[3]]<-c("species",paste("trait",letters[1:n.traits],sep="."))
  for (i in 1:n)
    Y[i,,"species"]<-sample(1:S,J,replace=TRUE,prob=survive[i,])
  
  for (i in 1:n) Y[i,,2:(n.traits+1)]<-as.matrix(traits.modus[Y[i,,"species"],])

  if (verbose)
    {
    cat("Community assembly... \n")
    pb <- utils::txtProgressBar (min = 0, max = sim.length, char = ".", width = 45, style = 3)
  }

  for (epoch in 1:sim.length)
  {
    for (j in 1:J)
    {
      for (i in 1:n) Y[i,,]<-Y[i,sample(1:J),]

      seed<-vector()
      for (i in 1:n)
      {
        compet <-do.call(competition.kernel,c(list(trait.values=Y[i,-1,-1]),list(...)))
        w<-do.call(fSeedProduction,c(list(compet=compet,abund=rep(1,J-1)),list(...)))
        # w is a vector with number of seeds produced by each induividuals
        w<-c(0,w) # the first individual died, so it does not produce seeds
        for (kk in 1:length(w))
          if (w[kk]>0)
            seed<-rbind(seed,c(i,Y[i,kk,]))
      }

      colnames(seed)<-c("site","species",paste("trait",letters[1:n.traits],sep="."))

      if (!is.null(fITV)) seed<-do.call(fITV,c(list(seeds=seed,n.traits=n.traits),list(...)))

      seed<-do.call(fDispersal,c(list(n=n,before=seed),list(...)))
      for (i in 1:n)
      {
       ns<-sum(seed[,"site"]==i)
       if (ns>1)
        {
          seed.local<-seed[seed[,"site"]==i,-1]
          survive <-do.call(fSurvive,c(list(trait.values=seed.local,env=x[i]),list(...)))
          Y[i,1,]<-seed.local[sample(1:nrow(seed.local),size=1,prob=survive),]
       }
       if (ns==1) Y[i,1,]<-seed[seed[,"site"]==i,-1]
      }

    }
    if (verbose) utils::setTxtProgressBar (pb, epoch)
  }
  if (verbose) close(pb)
  final.community<-vector()
  for (i in 1:n) final.community<-rbind(final.community,cbind(rep(i,J),Y[i,,]))
  colnames(final.community)<-c("site","species",paste("trait",letters[1:n.traits],sep="."))
  final.community<-as.data.frame(final.community)

  final.community[,"species"]<-paste("sp_",
                                     sprintf("%03d", final.community[,"species"]),
                                     sep="")
  final.community[,"site"]<-paste("site_",
                                     sprintf("%03d", final.community[,"site"]),
                                     sep="")
  parameters<-get("parameters",envir = comsimitvEnv)
  ret<-list(final.community=final.community,
            parameters=parameters,
            traits.modus=traits.modus)
  return(ret)
}
