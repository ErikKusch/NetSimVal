create_gradient_grid <- function(
    x_range, y_range, ncol = 100, nrow = 100,
    x_gradient = function(x) x,
    y_gradient = function(y) y) {
    # Validate output type
    output <- match.arg(output)

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


# 1. Simple linear gradient matrix
mat1 <- create_gradient_grid(
    x_range = c(0, 10), y_range = c(0, 10),
    ncol = 1e4, nrow = 1e4,
    x_gradient = function(x) x,
    y_gradient = function(y) y
)
image(z = mat1)


get_matrix_index_from_xy <- function(x, y, env_mat) {
    col <- which.min(abs(x_target - as.numeric(colnames(env_mat))))[1]
    row <- which.min(abs(y_target - as.numeric(rownames(env_mat))))[1]
    return(list(row = row, col = col))
}

x_target <- 9
y_target <- 9

# Find matrix indices
idx <- get_matrix_index_from_xy(
    x_target,
    y_target,
    env_mat = mat1
)
print(idx)

# Access value from matrix
value <- mat1[idx$row, idx$col]
print(value)
