#' Carrying Capacities
#'
#' A vector of carrying capacities.
#'
#' @format ## `CarryingK_vec`
#' A vector with 10 values (all values = 200).
"CarryingK_vec"

#' Bioclimatic Niche Preferences
#'
#' A vector of bioclimatic niche preferences.
#'
#' @format ## `Niches_vec`
#' A vector with 10 values.
"Niches_vec"

#' Initialising Information
#'
#' A data frame object containing initialising individuals.
#'
#' @format ## `Initialise_df`
#' Data frame containing 400 initialising individuals.
"Initialise_df"

#' Simulated Network
#'
#' Random network used to inform simulation framework.
#'
#' @format ## `Network_igraph`
#' igraph object containing 10 nodes and 23 links.
"Network_igraph"

#' Environmental Matrix
#'
#' An environmental matrix used for demographic simulations.
#'
#' @format ## `Env_mat`
#' Environmental matrix of 1e3 rows and columns with environmental range 0-10 in either dimension
"Env_mat"

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