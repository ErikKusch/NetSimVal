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
#' @return An matrix object with association/interaction strength stored as in cells.
#'
#' @importFrom randcorr randcorr
#' @importFrom stringr str_pad
#'
#' @author Erik Kusch, Natural History Museum, University of Oslo, Norway.
#'
#' @examples
#' Sim.Network(
#'   n_spec = 20,
#'   NetworkType = "Association",
#'   Sparcity = 0.5,
#'   MaxStrength = 1,
#'   seed = 42
#' )
#'
#' Sim.Network(
#'   n_spec = 10,
#'   NetworkType = "Interaction",
#'   Sparcity = 0,
#'   MaxStrength = 20,
#'   seed = 21
#' )
#' @export
Sim.Network <- function(n_spec = 5,
                        NetworkType = "Association", # or "Association"
                        Sparcity = 0.5,
                        MaxStrength = 1,
                        seed = 42) {
  set.seed(seed)

  if (!(NetworkType %in% c("Interaction", "Association"))) {
    stop('NetworkType needs to be either "Interaction" or "Association"')
  }

  Rand_corr <- randcorr(n_spec) * MaxStrength # establish random correlation matrix
  if (NetworkType == "Interaction") {
    diag(Rand_corr) <- 0
    Rand_corr[lower.tri(Rand_corr)] <- randcorr(n_spec)[lower.tri(randcorr(n_spec))]
    Rand_corr[sample(which(Rand_corr != 0), as.integer(sum(!is.na(Rand_corr)) * Sparcity))] <- 0
  }
  if (NetworkType == "Association") {
    # Rand_corr[lower.tri(Rand_corr)] <- 0 # make into undirected adjacency matrix representation
    diag(Rand_corr) <- 0
    idx <- which(upper.tri(Rand_corr, diag = FALSE), arr.ind = TRUE)
    # how many to zero
    k <- floor(nrow(idx) * Sparcity)
    # randomly choose positions
    zero_idx <- idx[sample(seq_len(nrow(idx)), k), , drop = FALSE]
    # set both (i,j) and (j,i) to zero
    for (r in seq_len(nrow(zero_idx))) {
      i <- zero_idx[r, 1]
      j <- zero_idx[r, 2]
      Rand_corr[i, j] <- 0
      Rand_corr[j, i] <- 0
    }
  }
  colnames(Rand_corr) <- rownames(Rand_corr) <- paste0("Sp_", stringr::str_pad(1:n_spec, width = 2, pad = "0"))
  # igraph::E(Rand_corr)$weight <- igraph::E(Rand_corr)$weight * MaxStrength
  return(Rand_corr)
}
