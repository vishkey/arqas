\name{MarkovianModel}
\alias{MarkovianModel}
\title{Define a queueing model}
\usage{
  MarkovianModel(arrivalDistribution = Exp(1),
    serviceDistribution = Exp(1))
}
\arguments{
  \item{arrdistr}{Arrival distribution (object of S4-class
  \code{distr} defined in \pkg{distr} package)}

  \item{serviceDistribution}{Service distribution (object
  of S4-class \code{distr} defined in \pkg{distr} package)}
}
\value{
  An object of class \code{MarkovianModel}, a list with the
  following components: \item{arrivalDistribution}{Arrival
  distribution} \item{serviceDistribution}{Service
  distribution}
}
\description{
  Constructor for \code{MarkovianModel} class.
}

