\name{G_G_INF}
\alias{G_G_INF}
\title{Obtains the main characteristics of a G/G/\eqn{\infty} model by simulation}
\usage{
  G_G_INF(arrivalDistribution = Exp(3),
    serviceDistribution = Exp(6), staClients = 100,
    nClients = 1000, historic = FALSE, nsim = 1, nproc = 1)
}
\arguments{
  \item{arrivalDistribution}{Arrival distribution}

  \item{serviceDistribution}{Service distribution}

  \item{staClients}{Number of customers used in
  stabilization stage}

  \item{nClients}{Number of customers used in the
  simulation stage}

  \item{historic}{Parameter to activate/deactivate the
  historic information}
}
\value{
  Returns the next information of a G/G/\eqn{\infty} model:
  \item{pn}{Vector of steady-state probabilities of having
  n customers in the system \eqn{P(n)}} \item{l}{Expected
  number of customers in the system \eqn{L}}
  \item{lq}{Expected number of customers in the queue
  \eqn{L_{q}}} \item{w}{Expected waiting time in the system
  \eqn{W}} \item{wq}{Expected waiting time in the queue
  \eqn{W_{q}}} \item{eff}{Efficiency of the system \eqn{Eff
  = W/(W-Wq)}} \item{rho}{Traffic intensity \eqn{\rho}}
  \item{historic}{Optional parameter that stores the
  evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq}
  during the simulation}
}
\description{
  Obtains the main characteristics of a G/G/\eqn{\infty}
  model by simulation
}
\seealso{
  Other SimulatedModels: \code{\link{ClosedNetwork}},
  \code{\link{G_G_1}}, \code{\link{G_G_1_INF_H}},
  \code{\link{G_G_1_K}}, \code{\link{G_G_S}},
  \code{\link{G_G_S_INF_H}}, \code{\link{G_G_S_INF_H_Y}},
  \code{\link{G_G_S_K}}, \code{\link{OpenNetwork}}
}

