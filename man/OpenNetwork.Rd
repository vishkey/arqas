\name{OpenNetwork}
\alias{OpenNetwork}
\title{Obtains the main characteristics of an Open Network model by simulation}
\usage{
  OpenNetwork(arr.distr, serv.distr, s, p, staClients,
    transitions, historic = FALSE)
}
\arguments{
  \item{arr.distr}{PairList indicating the arrival
  distribution and the node that uses it.}

  \item{serv.distr}{Vector of service distribution in each
  node}

  \item{s}{Vector of servers in each node}

  \item{p}{Routing matrix, where \eqn{p_{ij}} is the
  routing probability from node i to node j}

  \item{staClients}{Number of customers used in the
  stabilization stage}

  \item{transitions}{Number of transitions between nodes
  used in the simulation stage}

  \item{historic}{Parameter to activate/deactivate the
  historic information}
}
\value{
  Returns the next information of an Open network model:
  \item{pn}{Vector of steady-state probabilities of having
  n customers in the system \eqn{P(n)}} \item{l}{Vector of
  expected number of customers in the nodes \eqn{L}}
  \item{lq}{Vector of expected number of customers in the
  queues of the nodes \eqn{L_{q}}} \item{lqt}{Expected
  number of customers in all queues} \item{w}{Vector of
  expected waiting times in the nodes \eqn{W}}
  \item{wq}{Vector of expected waiting time in the queues
  of the nodes \eqn{W_{q}}} \item{eff}{Efficiency of the
  system \eqn{Eff = W/(W-Wq)}} \item{rho}{Traffic intensity
  \eqn{\rho}} \item{historic}{Optional parameter that
  stores the evolution of L, Lq, W and Wq during the
  simulation.}
}
\description{
  Obtains the main characteristics of an Open Network model
  by simulation
}
\seealso{
  Other SimulatedModels: \code{\link{ClosedNetwork}},
  \code{\link{G_G_1}}, \code{\link{G_G_1_INF_H}},
  \code{\link{G_G_1_K}}, \code{\link{G_G_INF}},
  \code{\link{G_G_S}}, \code{\link{G_G_S_INF_H}},
  \code{\link{G_G_S_INF_H_Y}}, \code{\link{G_G_S_K}}
}

