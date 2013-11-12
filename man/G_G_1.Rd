\name{G_G_1}
\alias{G_G_1}
\title{Obtains the main characteristics of a G/G/1 model by simulation}
\usage{
  G_G_1(arr.distr = Exp(1), serv.distr = Exp(1),
    staClients = 100, nClients = 1000, historic = FALSE)
}
\arguments{
  \item{arr.distr}{Arrival distribution}

  \item{serv.distr}{Service distribution}

  \item{staClients}{Number of customers used in the
  stabilization stage}

  \item{nClients}{Number of customers used in the
  simulation stage}

  \item{historic}{Parameter used to activate/deactivate the
  historic information}
}
\value{
  Returns the next information of a G/G/1 model:
  \item{pn}{Stores all the positives steady-state
  probabilities of having n customers, with n from 0 to
  staClients+nClients} \item{l}{Expected number of
  customers in the system \eqn{L}} \item{lq}{Expected
  number of customers in the queue \eqn{L_{q}}}
  \item{w}{Expected waiting time in the system \eqn{W}}
  \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
  \item{eff}{Efficiency of the system \eqn{Eff = W/(W-Wq)}}
  \item{rho}{Traffic intensity \eqn{\rho}}
  \item{historic}{Optional parameter that stores the
  evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and \eqn{Wq}
  during the simulation}
}
\description{
  Obtains the main characteristics of a G/G/1 model by
  simulation
}
\seealso{
  Other SimulatedModels: \code{\link{ClosedNetwork}},
  \code{\link{G_G_1_INF_H}}, \code{\link{G_G_1_K}},
  \code{\link{G_G_INF}}, \code{\link{G_G_S}},
  \code{\link{G_G_S_INF_H}}, \code{\link{G_G_S_INF_H_Y}},
  \code{\link{G_G_S_K}}, \code{\link{OpenNetwork}}
}

