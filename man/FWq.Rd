\name{FWq}
\alias{FWq}
\alias{FWq.M_M_1}
\alias{FWq.M_M_1_INF_H}
\alias{FWq.M_M_1_K}
\alias{FWq.M_M_INF}
\alias{FWq.M_M_S}
\alias{FWq.M_M_S_INF_H}
\alias{FWq.M_M_S_INF_H_Y}
\alias{FWq.M_M_S_K}
\alias{FWq.MarkovianModel}
\title{Distribution function of the waiting time in the queue}
\usage{
  FWq(qm, x)

  \method{FWq}{MarkovianModel} (qm, x)

  \method{FWq}{M_M_1} (qm, x)

  \method{FWq}{M_M_S} (qm, x)

  \method{FWq}{M_M_1_K} (qm, x)

  \method{FWq}{M_M_S_K} (qm, x)

  \method{FWq}{M_M_1_INF_H} (qm, x)

  \method{FWq}{M_M_S_INF_H} (qm, x)

  \method{FWq}{M_M_S_INF_H_Y} (qm, x)

  \method{FWq}{M_M_INF} (qm, x)
}
\arguments{
  \item{qm}{Queueing model}

  \item{x}{Time}
}
\value{
  \deqn{W_{q}(x)}
}
\description{
  Returns the value of the cumulative distribution function
  of waiting time in the queue \deqn{W_{q} = P(W_{q} \le
  x)}
}
\details{
  \code{FWq.MarkovianModel} implements the default method
  (generates a message)

  \code{FWq.M_M_1} implements the method for a M/M/1
  queueing model

  \code{FWq.M_M_S} implements the method for a M/M/S
  queueing model

  \code{FWq.M_M_1_K} implements the method for a M/M/1/K
  queueing model

  \code{FWq.M_M_S_K} implements the method for a M/M/S/K
  queueing model

  \code{FWq.M_M_1_INF_H} implements the method for a
  M/M/1/\eqn{\infty}/H queueing model

  \code{FWq.M_M_S_INF_H} implements the method for a
  M/M/s/\eqn{\infty}/H queueing model

  \code{FWq.M_M_S_INF_H_Y} implements the method for a
  M/M/s/\eqn{\infty}/H with Y replacements queueing model

  \code{FWq.M_M_INF} implements the method for a
  M/M/\eqn{\infty} queueing model
}

