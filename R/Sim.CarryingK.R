#' Generate Random Carrying Capacities per Species
#'
#' This function generates a vector of random numbers representing carrying capacities for the demographic simulation.
#'
#' @param n_spec Numeric. Number of species / size of node set.
#' @param k_range Numeric vector of length 2. Minimum and maximum values of carrying capacity range.
#' @param seed Numeric. Seed for random processes.
#'
#' @return A named vector of random carrying capacity values for each species.
#'
#' @importFrom stringr str_pad
#'
#' @author Erik Kusch, Natural History Museum, University of Oslo, Norway.
#'
#' @examples
#' Sim.CarryingK(n_spec = 20, k_range = c(200, 200), seed = 42)
#' Sim.CarryingK(n_spec = 20, k_range = c(100, 300), seed = 42)
#'
#' @export
Sim.CarryingK <- function(n_spec = 20, k_range = c(200, 200), seed = 42) {
  set.seed(seed)
  k_vec <- as.integer(runif(n = n_spec, min = k_range[1], max = k_range[2]))
  names(k_vec) <- paste0("Sp_", stringr::str_pad(1:n_spec, width = 2, pad = "0"))
  return(k_vec)
}
