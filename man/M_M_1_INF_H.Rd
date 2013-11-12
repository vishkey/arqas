\name{M_M_1_INF_H}
\alias{M_M_1_INF_H}
\title{Obtains the main characteristics of a M/M/1/\eqn{\infty}/H queueing model}
\usage{
  M_M_1_INF_H(lambda = 3, mu = 6, h = 5)
}
\arguments{
  \item{lambda}{Mean arrival rate}

  \item{mu}{Mean service rate}

  \item{h}{Population size}
}
\value{
  Returns the next information of a M/M/1/\eqn{\infty}/H
  model: \item{rho}{Constant \eqn{\lambda/\rho}}
  \item{barrho}{Traffic intensity \eqn{\bar{\rho}}}
  \item{barlambda}{Mean service rate \eqn{\bar{\lambda}}}
  \item{cn}{Constant coefficients used in the computation
  of \eqn{P(n)}} \item{p0}{Probability of empty system
  \eqn{P_{0}}} \item{l}{Expected number of customers in the
  system \eqn{L}} \item{lq}{Expected number of customers in
  the queue \eqn{L_{q}}} \item{w}{Expected waiting time in
  the system \eqn{W}} \item{wq}{Expected waiting time in
  the queue \eqn{W_{q}}} \item{eff}{Efficiency of the
  system \eqn{Eff = W/(W-W_q)}}
}
\description{
  Obtains the main characteristics of a
  M/M/1/\eqn{\infty}/H queueing model
}
\seealso{
  Other AnaliticalModels:
  \code{\link{ClosedJacksonNetwork}}, \code{\link{M_M_1}},
  \code{\link{M_M_1_K}}, \code{\link{M_M_INF}},
  \code{\link{M_M_S}}, \code{\link{M_M_S_INF_H}},
  \code{\link{M_M_S_INF_H_Y}}, \code{\link{M_M_S_K}},
  \code{\link{OpenJacksonNetwork}}
}

