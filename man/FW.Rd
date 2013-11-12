\name{FW}
\alias{FW}
\alias{FW.M_M_1}
\alias{FW.M_M_1_INF_H}
\alias{FW.M_M_1_K}
\alias{FW.M_M_INF}
\alias{FW.M_M_S}
\alias{FW.M_M_S_INF_H}
\alias{FW.M_M_S_INF_H_Y}
\alias{FW.M_M_S_K}
\alias{FW.MarkovianModel}
\title{Distribution function of the waiting time in the system}
\usage{
  FW(qm, x)

  \method{FW}{MarkovianModel} (qm, x)

  \method{FW}{M_M_1} (qm, x)

  \method{FW}{M_M_S} (qm, x)

  \method{FW}{M_M_1_K} (qm, x)

  \method{FW}{M_M_S_K} (qm, x)

  \method{FW}{M_M_1_INF_H} (qm, x)

  \method{FW}{M_M_S_INF_H} (qm, x)

  \method{FW}{M_M_S_INF_H_Y} (qm, x)

  \method{FW}{M_M_INF} (qm, x)
}
\arguments{
  \item{qm}{Queueing model}

  \item{x}{Time}
}
\value{
  \deqn{W(x)}
}
\description{
  Returns the value of the cumulative distribution function
  of the waiting time in the system for a queueing model
  \deqn{W(x) = P(W \le x)}
}
\details{
  \code{FW.MarkovianModel} implements the default method
  (generates a message)

  \code{FW.M_M_1} implements the method for a M/M/1
  queueing model

  \code{FW.M_M_S} implements the method for a M/M/S
  queueing model

  \code{FW.M_M_1_K} implements the method for a M/M/1/K
  queueing model

  \code{FW.M_M_S_K} implements the method for a M/M/S/K
  queueing model

  \code{FW.M_M_1_INF_H} implements the method for a
  M/M/1/\eqn{\infty}/H queueing model

  \code{FW.M_M_S_INF_H} implements the method for a
  M/M/s/\eqn{\infty}/H queueing model

  \code{FW.M_M_S_INF_H_Y} implements the method for a
  M/M/s/\eqn{\infty}/H/ with Y replacements queueing model

  \code{FW.M_M_INF} implements the method for a
  M/M/\eqn{\infty} queueing model
}

