\name{M_M_1_K}
\alias{M_M_1_K}
\title{Obtains the main characteristics of a M/M/1/K queueing model}
\usage{
  M_M_1_K(lambda = 3, mu = 6, k = 2)
}
\arguments{
  \item{lambda}{Mean arrival rate}

  \item{mu}{Mean service rate}

  \item{k}{Maximun size of the queue}
}
\value{
  Returns the next information of a M/M/1/K model:
  \item{rho}{Constant coefficient \eqn{\lambda/\rho}}
  \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
  \item{barlambda}{Effective arrival rate
  \eqn{\bar{\lambda}}} \item{l}{Expected mean number of
  customers in the system \eqn{L}} \item{lq}{Expected mean
  number of customers in the queue \eqn{L_{q}}}
  \item{w}{Expected waiting time in the system \eqn{W}}
  \item{wq}{Expected waiting time in the queue \eqn{W_{q}}}
  \item{eff}{Efficiency of the system \eqn{Eff =
  W/(W-W_q)}}
}
\description{
  Obtains the main characteristics of a M/M/1/K queueing
  model
}
\seealso{
  Other AnaliticalModels:
  \code{\link{ClosedJacksonNetwork}}, \code{\link{M_M_1}},
  \code{\link{M_M_1_INF_H}}, \code{\link{M_M_INF}},
  \code{\link{M_M_S}}, \code{\link{M_M_S_INF_H}},
  \code{\link{M_M_S_INF_H_Y}}, \code{\link{M_M_S_K}},
  \code{\link{OpenJacksonNetwork}}
}

