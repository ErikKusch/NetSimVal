#' Carrying Capacities
#'
#' A vector of carrying capacities.
#'
#' @format ## `CarryingK_vec`
#' A vector with 10 values (all values = 200).
"CarryingK_vec"

#' Initialising Information
#'
#' A list object containing data frame of initialising individuals and their bioclimatic niche preferences.
#'
#' @format ## `Initialise_ls`
#' List of 2. (1) data frame containing 400 initialising individuals, (2) containing 10 values of biolcimatic niche preferences.
"Initialise_ls"

#' Simulated Network
#'
#' Random network used to inform simulation framework.
#'
#' @format ## `Network_igraph`
#' igraph object containing 10 nodes and 23 links.
"Network_igraph"

#' Simulation Output
#'
#' Output of simulation framework.
#'
#' @format ## `SimulationOutput`
#' List of 99, each containing a dataframe of individuals alive at their respective simulation timestep.
"SimulationOutput"

#' Inferred Network
#'
#' Random network treated as though it was an inferred network for package documentation purposes.
#'
#' @format ## `Inferred_igraph`
#' igraph object containing 10 nodes and 23 links.
"Inferred_igraph"