% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gener.species.pool.R
\name{Gener.species.pool}
\alias{Gener.species.pool}
\title{Generating trait values for the species pool}
\usage{
Gener.species.pool(
  S,
  n.traits = 3,
  distribs = rep("unif", n.traits),
  distr.parms = list(),
  sigma = diag(1, n.traits, n.traits),
  ...
)
}
\arguments{
\item{S}{Species pool size}

\item{n.traits}{Number of traits}

\item{distribs}{Types of the distributions of traits}

\item{distr.parms}{Parameters of distribution (see Details)}

\item{sigma}{Matrix of variance-covariance matrix of traits}

\item{...}{Any additional parameters}
}
\value{
A data frame with traits as columns
}
\description{
It generates random trait values for species. Each species (individual)
are characterized by three traits.
}
\details{
Each species are characterized by three traits called trait A,
 B and C. Trait A describes the habitat preference,
 trait B influences the competitive interactions, while trait C
 is a completely neutral trait.

Any standard distribution of stats package can be used for generating
the random numbers. For list of these distribution see \code{\link[stats]{Distributions}}
In stats package the functions for the density/mass function are named
in the form dxxx."xxx" (without d!) as string (i.e. between quatation marks)
 should be supplied for parameter \code{distribs}.

 In this step single value of each trait is generated for each species,
 i.e. there is no intraspecific trait variation.

 If traits are independent (it is the default option), random number
 generating functions are called with parameters specified by the user.

 Otherwise, a variance-covariance matrix has to be given. First, triplets of
 random numbers are drawn from multivariate normal distribution with zero
 means and the supplied variance-covariance matrix as parameters.
 Then these random numbers are converted to probability by standard normal
 probability function, and then these probabilities converted to trait values
  using quantile function of selected distribution with parameters given
  by the user.
}
