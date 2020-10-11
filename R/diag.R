
# Diagonal matrices operations --------------------------------------------

diag_mul_vec <- function(diag_vec, vec) {
    diag_vec <- matrix(diag_vec, ncol = 1)
    vec <- matrix(vec, ncol = 1)
    return(diag_vec * vec)
}

diag_mul_diag <- function(diag1_vec, diag2_vec, mode = "vec") {
    if (mode == "vec") return(diag1_vec * diag2_vec)
    return(diag(diag1_vec * diag2_vec))
}

diag_solve <- function(diag_vec, mode = "vec") {
    stopifnot(prod(diag_vec) != 0)
    if (mode == "vec") return(1 / diag_vec)
    return(diag(1/diag_vec))
}

large_vec_norm <- function(vec, mode = "Gelder", p = 2) {
    if (mode == "Gelder") return(sum(vec^p)^(1/p))
    return(max(vec))
}

circulante_matrix <- function(first_row, direction = "r") {
    len <- length(first_row)
    matrix1 <- matrix(nrow = len, ncol = len)
    if (direction == "l") { 
        for (i in 1:len) {
            matrix1[i,] <- first_row[(((1:len) + (len - 2 + i)) %% len) + 1]
        }
    } else {
        for (i in 1:len) {
            matrix1[i,] <- first_row[(((1:len) - i) %% len) + 1]
        }
    }
    
    return(matrix1)
}
