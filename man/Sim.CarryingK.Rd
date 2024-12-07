% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Sim.CarryingK.R
\name{Sim.CarryingK}
\alias{Sim.CarryingK}
\title{Generate Random Carrying Capacities per Species}
\usage{
Sim.CarryingK(n_spec = 20, k_range = c(200, 200), seed = 42)
}
\arguments{
\item{n_spec}{Numeric. Number of species / size of node set.}

\item{k_range}{Numeric vector of length 2. Minimum and maximum values of carrying capacity range.}

\item{seed}{Numeric. Seed for random processes.}
}
\value{
A named vector of random carrying capacity values for each species.
}
\description{
This function generates a vector of random numbers representing carrying capacities for the demographic simulation.
}
\examples{
Sim.CarryingK(n_spec = 20, k_range = c(200,200), seed = 42)
Sim.CarryingK(n_spec = 20, k_range = c(100,300), seed = 42)

}
