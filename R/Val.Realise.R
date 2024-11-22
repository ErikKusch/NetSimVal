#' Generate Random Association/Interaction Network Matrices
#' 
#' This function generates a random network matrix representing either an undirected (association) or directed (interaction) network matrix according to user specifications.
#' 
#' @param n_spec Numeric. Number of species / size of node set.
#' @param NetworkType Character. Whether to create an interaction or association matrix. In an interaction matrix, columns act on rows. In an association matrix, entries are mirrored around the diagonal. Possible values: "Association" or "Interaction".
#' @param Sparcity Numeric. Proportion of unique associations or interactions which are exactly 0.
#' @param MaxStrength Numeric. Maximum value of associations or interactions.
#' @param seed Numeric. Seed for random processes.
#' 
#' @return An igraph object with association/interaction strength stored as "weight" attribute of edges.
#' 
#' @importFrom randcorr randcorr
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph E
#' 
#' @examples
#' Sim.Network(n_spec = 20,
#'             NetworkType = "Association",
#'             Sparcity = 0.5,
#'             MaxStrength = 1,
#'             seed = 42
#' )
#' 
#' Sim.Network(n_spec = 10,
#'             NetworkType = "Interaction",
#'             Sparcity = 0,
#'             MaxStrength = 20,
#'             seed = 21
#' )
#' @export
Sim.Network <- function(n_spec = 20,
                        NetworkType = "Association", # or "Association"
                        Sparcity = 0.5,
                        MaxStrength = 1,
                        seed = 42){
  set.seed(seed)
  
  if(!(NetworkType %in% c("Interaction", "Association"))){stop('NetworkType needs to be either "Interaction" or "Association"')}
  
  Rand_corr <- randcorr(n_spec) # establish random correlation matrix
  if(NetworkType == "Interaction"){
    Rand_corr[lower.tri(Rand_corr)] <- randcorr(n_spec)[lower.tri(randcorr(n_spec))]
    Rand_corr[sample(which(!is.na(Rand_corr)), as.integer(sum(!is.na(Rand_corr))*Sparcity))] <- 0
    Rand_corr <- graph_from_adjacency_matrix(adjmatrix = Rand_corr,
                                             mode = "directed",
                                             weighted = TRUE,
                                             diag = FALSE)
  }
  if(NetworkType == "Association"){
    Rand_corr[lower.tri(Rand_corr)] <- NA # make into undirected adjacency matrix representation
    diag(Rand_corr) <- NA
    Rand_corr[sample(which(!is.na(Rand_corr)), as.integer(sum(!is.na(Rand_corr))*Sparcity))] <- 0
    Rand_corr <- graph_from_adjacency_matrix(adjmatrix = Rand_corr,
                                             mode = "undirected",
                                             weighted = TRUE,
                                             diag = FALSE)
  }
  igraph::E(Rand_corr)$weight <-  igraph::E(Rand_corr)$weight*MaxStrength
  return(Rand_corr)
}
