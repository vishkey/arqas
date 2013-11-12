\name{Qn}
\alias{Qn}
\alias{Qn.M_M_1_INF_H}
\alias{Qn.M_M_1_K}
\alias{Qn.M_M_S_INF_H}
\alias{Qn.M_M_S_INF_H_Y}
\alias{Qn.M_M_S_K}
\alias{Qn.MarkovianModel}
\title{Steady-state probability of finding n customers in the system when a new customer arrives}
\usage{
  Qn(qm, n)

  \method{Qn}{MarkovianModel} (qm, n)

  \method{Qn}{M_M_1_K} (qm, n)

  \method{Qn}{M_M_S_K} (qm, n)

  \method{Qn}{M_M_1_INF_H} (qm, n)

  \method{Qn}{M_M_S_INF_H} (qm, n)

  \method{Qn}{M_M_S_INF_H_Y} (qm, n)
}
\arguments{
  \item{qm}{Queueing model}

  \item{n}{Customers}
}
\value{
  \deqn{Q(n)}
}
\description{
  Returns the probability of n customers in the system in
  the moment of the arrival of a customer.
}
\details{
  \code{Qn.MarkovianModel} implements the default method
  (generates a message).

  \code{Qn.M_M_1_K} implements the method for a M/M/1/K
  queueing model

  \code{Qn.M_M_S_K} implements the method for a M/M/S/K
  queueing model

  \code{Qn.M_M_1_INF_H} implements the method for a
  M/M/1/\eqn{\infty}/H queueing model

  \code{Qn.M_M_S_INF_H} implements the method for a
  M/M/s/\eqn{\infty}/H queueing model

  \code{Qn.M_M_S_INF_H_Y} implements the method for a
  M/M/s/\eqn{\infty}/H with Y replacements queueing model
}

