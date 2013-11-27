\name{ClosedJacksonNetwork}
\alias{ClosedJacksonNetwork}
\title{Obtains the main characteristics of a Closed Jackson Network model}
\usage{
  ClosedJacksonNetwork(mu = c(5, 5, 10, 15),
    s = c(2, 2, 1, 1),
    p = array(c(0.25, 0.15, 0.5, 0.4, 0.15, 0.35, 0.25, 0.3, 0.2, 0.2, 0.15, 0.25, 0.4, 0.3, 0.1, 0.05), dim = c(4, 4)),
    n = 10)
}
\arguments{
  \item{mu}{Vector of mean service rates}

  \item{s}{Vector of servers at each node}

  \item{p}{Routing matrix, where \eqn{p_{ij}} is the
  routing probability from node i to node j}

  \item{n}{Number of customers in the network}
}
\value{
  Returns the next information of a Closed Jackson Network
  model: \item{rho}{Traffic intensity \eqn{\rho}}
  \item{l}{Expected number of customers in the system
  \eqn{L}} \item{lq}{Expected number of customers in the
  queue \eqn{L_{q}}} \item{w}{Expected waiting time in the
  system \eqn{W}} \item{wq}{Expected waiting time in the
  queue \eqn{W_{q}}} \item{eff}{Efficiency of the system
  \eqn{Eff = W/(W-W_q)}}
}
\description{
  Obtains the main characteristics of a Closed Jackson
  Network model
}
\seealso{
  Other AnaliticalModels: \code{\link{M_M_1}},
  \code{\link{M_M_1_INF_H}}, \code{\link{M_M_1_K}},
  \code{\link{M_M_INF}}, \code{\link{M_M_S}},
  \code{\link{M_M_S_INF_H}}, \code{\link{M_M_S_INF_H_Y}},
  \code{\link{M_M_S_K}}, \code{\link{OpenJacksonNetwork}}
}

