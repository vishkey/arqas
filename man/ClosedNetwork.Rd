\name{ClosedNetwork}
\alias{ClosedNetwork}
\title{Obtains the main characteristics of a Closed Network model by simulation}
\usage{
  ClosedNetwork(serv.distr, s, p, staClients, nClients,
    transitions, historic = FALSE)
}
\arguments{
  \item{serv.distr}{Service distributions for the nodes of
  the network}

  \item{s}{Vector of servers at each node}

  \item{p}{Routing matrix, where \eqn{p_{ij}} is the
  routing probability from node i to node j}

  \item{staClients}{Number of customers used in the
  stabilization stage}

  \item{nClients}{Number of customers in the system}

  \item{transitions}{Number of transitions between nodes
  used in the simulation stage}

  \item{historic}{Parameter to activate/deactivate the
  historic information}
}
\value{
  Returns the next information of a Closed Network model:
  \item{pn}{Vector of steady-state probabilities of having
  n customers in the system \eqn{P(n)}} \item{l}{Vector of
  expected number of customers in the nodes \eqn{L}}
  \item{lq}{Vector of expected number of customers in the
  queues of the nodes \eqn{L_{q}}} \item{lqt}{Expected
  number of customers in the all queues} \item{w}{Vector of
  expected waiting times in the nodes \eqn{W}}
  \item{wq}{Vector of expected waiting times in the queues
  of the nodes \eqn{W_{q}}} \item{eff}{Efficiency of the
  system \eqn{Eff = W/(W-Wq)}} \item{rho}{Traffic intensity
  \eqn{\rho}} \item{historic}{Optional parameter that
  stores the evolution of \eqn{L}, \eqn{Lq}, \eqn{W} and
  \eqn{Wq} during the simulation.}
}
\description{
  Obtains the main characteristics of a Closed Network
  model by simulation
}
\seealso{
  Other SimulatedModels: \code{\link{G_G_1}},
  \code{\link{G_G_1_INF_H}}, \code{\link{G_G_1_K}},
  \code{\link{G_G_INF}}, \code{\link{G_G_S}},
  \code{\link{G_G_S_INF_H}}, \code{\link{G_G_S_INF_H_Y}},
  \code{\link{G_G_S_K}}, \code{\link{OpenNetwork}}
}

