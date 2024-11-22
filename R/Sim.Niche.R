#' Generate Random Carrying Capacities per Species
#' 
#' This function generates a vector of random numbers representing carrying capacities for the demographic simulation.
#' 
#' @param n_spec Numeric. Number of species / size of node set.
#' @param Env_range Numeric vector of length 2. Minimum and maximum values of environmental gradient.
#' @param seed Numeric. Seed for random processes.
#' 
#' @return A named vector of random carrying capacity values for each species.
#' 
#' @examples
#' Sim.Niche(n_spec = 20, Env_range = c(0,10), seed = 42)
#' Sim.Niche(n_spec = 15, Env_range = c(11,33), seed = 42)
#' 
#' @export
Sim.Niche <- function(n_spec = 20, Env_range = c(0,10), seed = 42){
  set.seed(seed)
  Trait_means <- runif(n = n_spec, min = Env_range[1], max = Env_range[2])
  names(Trait_means) <- paste("Sp", 1:n_spec, sep = "_")
  return(Trait_means)
}

