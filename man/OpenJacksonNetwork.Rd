\name{OpenJacksonNetwork}
\alias{OpenJacksonNetwork}
\title{Obtains the main characteristics of an Open Jackson network model}
\usage{
  OpenJacksonNetwork(lambda = c(20, 30), mu = c(100, 25),
    s = c(1, 2),
    p = matrix(c(0.2, 0.25, 0.1, 0), nrow = 2, ncol = 2))
}
\arguments{
  \item{lambda}{Vector of arrival rates at each node}

  \item{mu}{Vector of mean service rates}

  \item{s}{Vector with the number of servers at each node}

  \item{p}{Routing matrix, where \eqn{p_{ij}} is the
  routing probability from node i to node j}
}
\value{
  Returns the next information of an Open Jackson network
  model: \item{rho}{Traffic intensity \eqn{\rho}}
  \item{l}{Vector with the expected number of customers in
  the nodes \eqn{L}} \item{lq}{Vector of expected number of
  customers in the queue of each node \eqn{L_{q}}}
  \item{w}{Vector of expected waiting time in each node
  \eqn{W}} \item{wq}{Vector of expected waiting time in the
  queue of each node \eqn{W_{q}}} \item{lt}{Expected number
  of customers in the network} \item{lqt}{Expected number
  of customers in all queues} \item{wt}{Total expected
  waiting time in the network} \item{wqt}{Total expected
  waiting time in all queues} \item{eff}{Efficiency
  \eqn{Eff = W/(W-W_q)}}
}
\description{
  Obtains the main characteristics of an Open Jackson
  network model
}
\seealso{
  Other AnaliticalModels:
  \code{\link{ClosedJacksonNetwork}}, \code{\link{M_M_1}},
  \code{\link{M_M_1_INF_H}}, \code{\link{M_M_1_K}},
  \code{\link{M_M_INF}}, \code{\link{M_M_S}},
  \code{\link{M_M_S_INF_H}}, \code{\link{M_M_S_INF_H_Y}},
  \code{\link{M_M_S_K}}
}

