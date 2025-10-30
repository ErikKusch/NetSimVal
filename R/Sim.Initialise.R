#' Generate Initialising Individuals for Demographic Simulation
#'
#' This function generates a set of individuals placed at random across the simulation landscape in x and y dimensions. It also applies species specific niche preferences to each individual with respect to their species membership.
#'
#' @param n_spec Numeric. Number of species / size of node set.
#' @param n_individuals Numeric. How many initialising individuals to randomly create.
#' @param n_mode Character. How to use n_individuals. If "each" the the number of initialising individuals will be equal to n_individual multiplied by n_spec. If "total", n_individuals are created with random species memberships within n_spec.
#' @param Env_range Numeric vector of length 2. Minimum and maximum values of environmental gradient.
#' @param Trait_means Numeric. Mean for bioclimatic niche preference of each species. Must be same length as number specified in n_spec.
#' @param Trait_sd Numeric. Standard deviation of trait distributions from which trait values of each individual are drawn.
#' @param seed Numeric. Seed for random processes.
#'
#' @return A data frame - the data frame of initialising individuals with columns ID, Trait, X, Y, and Species.
#' 
#' @author Erik Kusch, Natural History Museum, University of Oslo, Norway.
#' 
#' @examples
#' Niches_vec <- Sim.Niche(n_spec = 10, Env_range = c(0, 10), seed = 42)
#' Sim.Initialise(
#'   n_spec = 10, n_individuals = 4e2, n_mode = "each",
#'   Env_range = c(0, 10), Trait_means = Niches_vec, Trait_sd = 1, seed = 42
#' )
#' Sim.Initialise(
#'   n_spec = 20, n_individuals = 4e2, n_mode = "total",
#'   Env_range = c(10, 100), Trait_means = Niches_vec, Trait_sd = 2.5, seed = 42
#' )
#'
#' @export
Sim.Initialise <- function(
    n_spec = 10,
    n_individuals = 4e2,
    n_mode = "each", # or "total"
    Env_range = c(0, 10),
    Trait_means,
    Trait_sd = 1,
    seed = 42) {
  set.seed(seed)

  if (!(n_mode %in% c("each", "total"))) {
    stop('n_mode needs to be either "each" or "total"')
  }

  ## generating IDs and species memeberships
  if (n_mode == "each") {
    ID <- 1:(n_individuals * n_spec)
    Species <- rep(paste("Sp", 1:n_spec, sep = "_"), each = n_individuals)
  } else {
    ID <- 1:n_individuals
    Species <- sample(paste("Sp", 1:n_spec, sep = "_"), size = n_individuals, replace = TRUE)
  }

  ## creating data frame to hold individuals
  ID_df <- data.frame(ID = ID, Trait = NA, X = NA, Y = NA, Species = Species)

  ## sampling trait values and assigning to data frame
  for (x in names(Trait_means)) {
    ID_df$Trait[ID_df$Species == x] <- rnorm(n = length(x %in% ID_df$Species), mean = Trait_means[x], sd = Trait_sd)
  }

  ## individual locations
  ID_df$X <- runif(n = nrow(ID_df), min = Env_range[1], max = Env_range[2])
  ID_df$Y <- runif(n = nrow(ID_df), min = Env_range[1], max = Env_range[2])

  ## return object
  return(ID_df)
}
