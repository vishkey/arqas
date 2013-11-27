\name{G_G_S_K}
\alias{G_G_S_K}
\title{Obtains the main characteristics of a G/G/s/K model by simulation}
\usage{
  G_G_S_K(arrivalDistribution = Exp(1),
    serviceDistribution = Exp(1), s = 2, K = 3,
    staClients = 100, nClients = 1000, historic = FALSE)
}
\arguments{
  \item{arrivalDistribution}{Arrival distribution}

  \item{serviceDistribution}{Service distribution}

  \item{s}{Number of servers}

  \item{K}{Maximun size of the queue}

  \item{staClients}{Number of customers used in the
  stabilization stage}

  \item{nClients}{Number of customers used in the
  simulation stage}

  \item{historic}{Parameter to activate/deactivate the
  historic information}
}
\value{
  Returns the next information of a G/G/S/K model:
  \item{pn}{Vector of steady-state probabilities of having
  n customers in the system \eqn{P(n)}} \item{l}{Expected
  number of customers in the system \eqn{L}}
  \item{lq}{Expected number of customers in the queue
  \eqn{L_{q}}} \item{w}{Expected waiting time in the system
  \eqn{W}} \item{wq}{Expected waiting time in the queue
  \eqn{W_{q}}} \item{eff}{Efficiency of the system \eqn{Eff
  = W/(W-Wq)}} \item{rho}{Traffic intensity \eqn{\rho}}
  \item{historic}{Optional parameter that stores the
  evolution of \eqn{L}, \eqn{Lq}, \eqn{W}and \eqn{Wq}
  during the simulation}
}
\description{
  Obtains the main characteristics of a G/G/s/K model by
  simulation
}
\seealso{
  Other SimulatedModels: \code{\link{ClosedNetwork}},
  \code{\link{G_G_1}}, \code{\link{G_G_1_INF_H}},
  \code{\link{G_G_1_K}}, \code{\link{G_G_INF}},
  \code{\link{G_G_S}}, \code{\link{G_G_S_INF_H}},
  \code{\link{G_G_S_INF_H_Y}}, \code{\link{OpenNetwork}}
}

