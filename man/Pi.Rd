\name{Pi}
\alias{Pi}
\alias{Pi.ClosedJackson}
\title{Steady-state probability of n customers in the node i of a network.}
\usage{
  Pi(net, n, node)

  \method{Pi}{ClosedJackson} (net, n, node)
}
\arguments{
  \item{net}{Network}

  \item{n}{Customers}

  \item{node}{Node}
}
\value{
  P(n) in the selected node
}
\description{
  Returns the value \eqn{P_{i}(n)} in the node i of a
  network
}
\details{
  \code{Pi.ClosedJackson} implements the method for a
  Closed Jackson Network model
}

