\name{summarySimple}
\alias{summarySimple}
\title{Shows a plot of the evolution of a variable during the simulation}
\usage{
  summarySimple(object, minrange, maxrange, var,
    graphics = "ggplot2")
}
\arguments{
  \item{object}{Simulated Model}

  \item{minrange}{Start client of simulation}

  \item{maxrange}{Last client of simulation}

  \item{var}{Variable to show graphic (L, Lq, W, Wq)}

  \item{graphics}{Type of graphics: "graphics" use the
  basic R plot and "ggplot2" the library ggplot2}
}
\description{
  Shows a plot of the evolution of a variable during the
  simulation
}

