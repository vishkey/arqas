\name{summary.MarkovianModel}
\alias{summary.MarkovianModel}
\title{Shows the main graphics of the parameters of a Markovian Model}
\usage{
  \method{summary}{MarkovianModel} (object,
    t = list(range = seq(object$out$w, object$out$w * 3, length.out = 100)),
    n = c(0:5), ...)
}
\arguments{
  \item{object}{Markovian Model}

  \item{t}{Range of t}

  \item{n}{Range of n}

  \item{...}{Further arguments passed to or from other
  methods.}
}
\description{
  Shows the main graphics of the parameters of a Markovian
  Model
}

