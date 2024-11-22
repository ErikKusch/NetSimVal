#' Network Inference Error Rates Validation
#' 
#' Calculate accuracy of network inference by sign of underlying true link.
#' 
#' @param Network1 igraph object. True network in comparison.
#' @param Network2 igraph object. Inferred network in comparison.
#' 
#' @return A dataframe containing percentages and error rate types:
#'    - TP ... true positive associations; how many of the inferred positive associations are correctly identified as such? - NaN means no positives in the inferred matrix
#'    - TN ... true negative associations; how many of the inferred negative associations are correctly identified as such? - NaN means no negatives in the inferred matrix
#'    - TA ... true absent; how many of the inferred absent interactions are correctly identified as such? - NaN means no absents in the inferred matrix
#'    - FP ... falsely inferred positive; how many of the inferred positive associations are incorrectly identified as such?
#'    - FN ... falsely inferred negative; how many of the inferred negative associations are incorrectly identified as such?
#'    - FA ... false absent; how many of the inferred absent interactions are correctly identified as such?
#'    - MP ... true positive links which were not inferred; how many of the true positive associations were not inferred?
#'    - MN ... true negative links which were not inferred; how many of the true negative associations were not inferred?
#'    - MA ... true absent links which were not inferred; how many of the true absent associations were not inferred?
#'  
#' @importFrom igraph V
#' @importFrom igraph as_adjacency_matrix
#' @importFrom igraph is_directed
#' 
#' @examples
#' data("Network_igraph")
#' data("Inferred_igraph")
#' Val.ErrorRates(Network1 = Network_igraph, Network2 = Inferred_igraph)
#' 
#' @export
Val.ErrorRates <- function(Network1, Network2){
  ## matrices
  matrices_ls <- lapply(list(Network1, Network2), FUN = function(NetworkI){
    if(is.null(igraph::V(NetworkI)$names)){
      igraph::V(NetworkI)$names <- paste0("Sp_", igraph::V(NetworkI))
    }
    mat <- as.matrix(as_adjacency_matrix(NetworkI, attr = "weight"))
    colnames(mat) <- rownames(mat) <- igraph::V(NetworkI)$names
    diag(mat) <- NA
    if(!is_directed(NetworkI)){
      mat[lower.tri(mat) ] <- NA
    }
    sign(mat)
  })
  True_mat <- matrices_ls[[1]]
  Inference_mat <- matrices_ls[[2]]
  
  ## error rates  
  errorrates_df <- data.frame(
    Values = c(
      # true positive associations; how many of the inferred positive associations are correctly identified as such? - NaN means no positives in the inferred matrix
      TP = sum((True_mat + Inference_mat) == 2, na.rm = TRUE)/
        sum(Inference_mat == 1, na.rm = TRUE), 
      # true negative associations; how many of the inferred negative associations are correctly identified as such? - NaN means no negatives in the inferred matrix
      TN = sum((True_mat + Inference_mat) == -2, na.rm = TRUE)/
        sum(Inference_mat == -1, na.rm = TRUE),
      # true absent; how many of the inferred absent interactions are correctly identified as such? - NaN means no absents in the inferred matrix
      TA = sum((True_mat == 0) + (Inference_mat == 0) == 2, na.rm = TRUE)/
        sum(Inference_mat == 0, na.rm = TRUE),
      
      # falsely inferred positive; how many of the inferred positive associations are incorrectly identified as such?
      FP = sum((True_mat != 1) + (Inference_mat == 1) == 2, na.rm = TRUE)/
        sum(Inference_mat == 1, na.rm = TRUE), 
      # falsely inferred negative; how many of the inferred negative associations are incorrectly identified as such?
      FN = sum((True_mat != -1) + (Inference_mat == -1) == 2, na.rm = TRUE)/
        sum(Inference_mat == -1, na.rm = TRUE),
      # false absent; how many of the inferred absent interactions are correctly identified as such?
      FA = sum((True_mat != 0) + (Inference_mat == 0) == 2, na.rm = TRUE)/
        sum(Inference_mat == 0, na.rm = TRUE),
      
      # true positive links which were not inferred; how many of the true positive associations were not inferred?
      MP = (sum((True_mat == 1) + (Inference_mat != 1) == 2, na.rm = TRUE)/
              sum(True_mat == 1, na.rm = TRUE)), 
      # true negative links which were not inferred; how many of the true negative associations were not inferred?
      MN = (sum((True_mat == -1) + (Inference_mat != -1) == 2, na.rm = TRUE)/
              sum(True_mat == -1, na.rm = TRUE)),
      # true absent links which were not inferred; how many of the true absent associations were not inferred?
      MA = (sum((True_mat == 0) + (Inference_mat != 0) == 2, na.rm = TRUE)/
              sum(True_mat == 0, na.rm = TRUE))
      
    ),
    Metric = c("TP", "TN", "TA", "FP", "FN", "FA", "MP", "MN", "MA")
  )
  return(errorrates_df)
}
