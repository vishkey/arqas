\name{maxCustomers}
\alias{maxCustomers}
\alias{maxCustomers.M_M_1_INF_H}
\alias{maxCustomers.M_M_1_K}
\alias{maxCustomers.M_M_S_INF_H}
\alias{maxCustomers.M_M_S_INF_H_Y}
\alias{maxCustomers.M_M_S_K}
\alias{maxCustomers.MarkovianModel}
\title{Returns the maximun value for n that satisfies P(n) > 0}
\usage{
  maxCustomers(qm)

  \method{maxCustomers}{MarkovianModel} (qm)

  \method{maxCustomers}{M_M_1_K} (qm)

  \method{maxCustomers}{M_M_S_K} (qm)

  \method{maxCustomers}{M_M_1_INF_H} (qm)

  \method{maxCustomers}{M_M_S_INF_H} (qm)

  \method{maxCustomers}{M_M_S_INF_H_Y} (qm)
}
\arguments{
  \item{qm}{object MarkovianModel}
}
\description{
  Returns the maximun value for n that satisfies P(n) > 0
}
\details{
  \code{maxCustomers.MarkovianModel} implements the default
  method. Returns infinite.

  \code{maxCustomers.M_M_1_K} implements the method for a
  M/M/1/K queueing model

  \code{maxCustomers.M_M_S_K} implements the method for a
  M/M/S/K queueing model

  \code{maxCustomers.M_M_1_INF_H} implements the method for
  a M/M/1/\eqn{\infty}/H queueing model

  \code{maxCustomers.M_M_S_INF_H} implements the method for
  a M/M/s/\eqn{\infty}/H queueing model

  \code{maxCustomers.M_M_S_INF_H_Y} implements the method
  for a M/M/s/\eqn{\infty}/H/Y queueing model
}

