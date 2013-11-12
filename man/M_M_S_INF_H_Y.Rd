\name{M_M_S_INF_H_Y}
\alias{M_M_S_INF_H_Y}
\title{Obtains the main characteristics of a M/M/s/\eqn{\infty}/H with Y replacements queueing model}
\usage{
  M_M_S_INF_H_Y(lambda = 3, mu = 6, s = 3, h = 5, y = 3)
}
\arguments{
  \item{lambda}{Mean arrival rate}

  \item{mu}{Mean service rate}

  \item{s}{Number of servers}

  \item{h}{Population size}

  \item{y}{Number of replacements}
}
\value{
  Returns the next information of a M/M/s/\eqn{\infty}/H/Y
  model: \item{rho}{Constant coefficient
  \eqn{\lambda/\rho}} \item{barrho}{Traffic intensity
  \eqn{\bar{\rho}}} \item{barlambda}{Effective arrival rate
  \eqn{\bar{\lambda}}} \item{cn}{Constant coefficients used
  in the computation of \eqn{P(n)} \eqn{C_{n}}}
  \item{p0}{Probability of 0 customers in the system
  \eqn{P_{0}}} \item{l}{Expected number of customers in the
  system \eqn{L}} \item{lq}{Expected number of customers in
  the queue \eqn{L_{q}}} \item{w}{Expected waiting time in
  the system \eqn{W}} \item{wq}{Expected waiting time in
  the queue \eqn{W_{q}}} \item{eff}{Efficiency of the
  system \eqn{Eff = W/(W-W_q)}}
}
\description{
  Obtains the main characteristics of a
  M/M/s/\eqn{\infty}/H with Y replacements queueing model
}
\seealso{
  Other AnaliticalModels:
  \code{\link{ClosedJacksonNetwork}}, \code{\link{M_M_1}},
  \code{\link{M_M_1_INF_H}}, \code{\link{M_M_1_K}},
  \code{\link{M_M_INF}}, \code{\link{M_M_S}},
  \code{\link{M_M_S_INF_H}}, \code{\link{M_M_S_K}},
  \code{\link{OpenJacksonNetwork}}
}

