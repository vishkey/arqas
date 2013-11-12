\name{M_M_S_K}
\alias{M_M_S_K}
\title{Obtains the main characteristics of a M/M/S/k queueing model}
\usage{
  M_M_S_K(lambda = 3, mu = 6, s = 2, k = 3)
}
\arguments{
  \item{lambda}{Mean arrival rate}

  \item{mu}{Mean service rate}

  \item{s}{Number of servers}

  \item{k}{Maximun size of the queue}
}
\value{
  Returns the next information of a M/M/S/K model:
  \item{rho}{Constant coefficient \eqn{\lambda/\rho}}
  \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
  \item{barlambda}{Effective arrival rate
  \eqn{\bar{\lambda}}} \item{cn}{Constant coefficients used
  in the computation of \eqn{P(n)}} \item{pks}{Probability
  of \eqn{K+s} customers in the system \eqn{P_{K+s}}}
  \item{p0}{Probability of empty system \eqn{P_{0}}}
  \item{l}{Expected number of customers in the system
  \eqn{L}} \item{lq}{Expected number of customers in the
  queue\eqn{L_{q}}} \item{w}{Expected waiting time in the
  system \eqn{W}} \item{wq}{Expected waiting time in the
  queue \eqn{W_{q}}} \item{eff}{Efficiency of the system
  \eqn{Eff = W/(W-W_q)}}
}
\description{
  Obtains the main characteristics of a M/M/S/k queueing
  model
}
\seealso{
  Other AnaliticalModels:
  \code{\link{ClosedJacksonNetwork}}, \code{\link{M_M_1}},
  \code{\link{M_M_1_INF_H}}, \code{\link{M_M_1_K}},
  \code{\link{M_M_INF}}, \code{\link{M_M_S}},
  \code{\link{M_M_S_INF_H}}, \code{\link{M_M_S_INF_H_Y}},
  \code{\link{OpenJacksonNetwork}}
}

