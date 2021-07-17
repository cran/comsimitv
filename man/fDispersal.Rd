\name{fDispersal}
\alias{fDispersal}
\title{User defined functions for dispersal}
\arguments{
\item{before}{matrix where each seed is represented by one row, and seeds's attributes (location, species, trait values) are in the columns}

\item{...}{Additional parameters of functions called by the framework.}
}
\value{
Same type as the first argument.
}
\description{
These functions define how seeds can spread among local communities.
}
\details{
User can define any rule for seed dispersal. The only requirement is that both first argument and value of the function should be a matrix where each seed is represented by one row, and seeds's attributes (location, species, trait values) are in the columns. The locality information has to be stored in column named 'site'.


Available function in the package:

\code{\link{MetaCom.Dispersal}}

}
