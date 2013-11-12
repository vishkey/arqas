\name{MarkovianModel}
\alias{MarkovianModel}
\title{Define a queueing model}
\usage{
  MarkovianModel(arr.distr = Exp(1), serv.distr = Exp(1))
}
\arguments{
  \item{arr.distr}{Arrival distribution (object of S4-class
  \code{distr} defined in \pkg{distr} package)}

  \item{serv.distr}{Service distribution (object of
  S4-class \code{distr} defined in \pkg{distr} package)}
}
\value{
  An object of class \code{MarkovianModel}, a list with the
  following components: \item{arr.distr}{Arrival
  distribution} \item{serv.distr}{Service distribution}
}
\description{
  Constructor for \code{MarkovianModel} class.
}

