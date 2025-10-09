#' Generate Spatial Environment as Matrix
#'
#' This function generates a matrix that serves as a space-definition for demographic simulations.
#'
#' @param x_range Numeric vector of length 2 specifying the minimum and maximum x-coordinates
#' @param y_range Numeric vector of length 2 specifying the minimum and maximum y-coordinates
#' @param ncol Integer specifying the number of columns (resolution in x direction). Default is 100.
#' @param nrow Integer specifying the number of rows (resolution in y direction). Default is 100.
#' @param x_gradient Function that defines how values change along the x-axis. Default is the identity function resulting in a linear gradient.
#' @param y_gradient Function that defines how values change along the y-axis. Default is the identity function resulting in a linear gradient.
#'
#' @return A matrix where cells report environmental values while column and row names report X and Y coordinates, respectively.
#'
#' @examples
#' # 1. Simple linear gradient matrix
#' mat <- Sim.Space(
#'     x_range = c(0, 10), y_range = c(0, 10),
#'     ncol = 1e4, nrow = 1e4,
#'     x_gradient = function(x) x,
#'     y_gradient = function(y) y
#' )
#'
#' # 2. Homogenuous matrix
#' mat <- Sim.Space(
#'     x_range = c(0, 10), y_range = c(0, 10),
#'     ncol = 1e4, nrow = 1e4,
#'     x_gradient = function(x) 1,
#'     y_gradient = function(y) 1
#' )
#' @export
Sim.Space <- function(
    x_range, y_range, ncol = 100, nrow = 100,
    x_gradient = function(x) x,
    y_gradient = function(y) y) {

    # Create sequences for X and Y axes
    x_seq <- seq(x_range[1], x_range[2], length.out = ncol)
    y_seq <- seq(y_range[1], y_range[2], length.out = nrow)

    # Apply gradient functions
    x_vals <- x_gradient(x_seq)
    if (length(x_vals) == 1) {
        x_vals <- rep(x_vals, nrow)
    }
    y_vals <- y_gradient(y_seq)
    if (length(y_vals) == 1) {
        y_vals <- rep(y_vals, ncol)
    }

    # Outer product to create 2D grid: each cell is x + y (or x * y, etc.)
    grid <- outer(y_vals, x_vals, "+") # or "*", depending on desired effect
    colnames(grid) <- x_seq
    rownames(grid) <- y_seq

    return(grid)
}
