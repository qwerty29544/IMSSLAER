
# Dependencies ------------------------------------------------------------

require(torch)

# Diagonal matrices operations --------------------------------------------

#' Title
#'
#' @param diag_vec 
#' @param vec 
#'
#' @return
#' @export
#'
#' @examples
diag_mul_vec <- function(diag_vec, vec) {
    diag_vec <- matrix(diag_vec, ncol = 1)
    vec <- matrix(vec, ncol = 1)
    return(diag_vec * vec)
}

#' Title
#'
#' @param diag1_vec 
#' @param diag2_vec 
#' @param mode 
#'
#' @return
#' @export
#'
#' @examples
diag_mul_diag <- function(diag1_vec, diag2_vec, mode = "vec") {
    if (mode == "vec") return(diag1_vec * diag2_vec)
    return(diag(diag1_vec * diag2_vec))
}

#' Title
#'
#' @param diag_vec 
#' @param mode 
#'
#' @return
#' @export
#'
#' @examples
diag_solve <- function(diag_vec, mode = "vec") {
    stopifnot(prod(diag_vec) != 0)
    if (mode == "vec") return(1 / diag_vec)
    return(diag(1/diag_vec))
}

#' Title
#'
#' @param vec 
#' @param mode 
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
large_vec_norm <- function(vec, mode = "Gelder", p = 2) {
    if (mode == "Gelder") return(sum(vec^p)^(1/p))
    return(max(vec))
}

#' Title
#'
#' @param first_row 
#' @param direction 
#'
#' @return
#' @export
#'
#' @examples
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



# Torch algos -------------------------------------------------------------


#' Title
#'
#' @param diag_vec 
#' @param vec 
#'
#' @return
#' @export
#'
#' @examples
diag_mul_vec_torch <- function(diag_vec, vec) {
    require(torch)
    if (torch::cuda_is_available()) {
        device <- "cuda:0"
    } else {
        device <- "cpu:0"
    }
    diag_vec <- torch::torch_tensor(diag_vec)$reshape(c(1, length(diag_vec)))$to(device)
    vec <- torch::torch_tensor(vec)$to(device)
    return(as_array((diag_vec * vec)$to(device = "cpu:0")))
}

#' Title
#'
#' @param diag1_vec 
#' @param diag2_vec 
#' @param mode 
#'
#' @return
#' @export
#'
#' @examples
diag_mul_diag_torch <- function(diag1_vec, diag2_vec, mode = "vec") {
    require(torch)
    if (torch::cuda_is_available()) {
        device <- "cuda:0"
    } else {
        device <- "cpu:0"
    }
    diag1_vec <- torch::torch_tensor(diag1_vec)$reshape(c(1, length(diag1_vec)))$to(device)
    diag2_vec <- torch::torch_tensor(diag2_vec)$reshape(c(1, length(diag2_vec)))$to(device)
    if (mode == "vec") return(as_array((diag1_vec * diag2_vec)$to(device = "cpu:0")))
    return(diag(as_array((diag1_vec * diag2_vec)$to(device = "cpu:0"))))
}

#' Title
#'
#' @param diag_vec 
#' @param mode 
#'
#' @return
#' @export
#'
#' @examples
diag_solve_torch <- function(diag_vec, mode = "vec") {
    stopifnot(prod(diag_vec) != 0)
    require(torch)
    if (torch::cuda_is_available()) {
        device <- "cuda:0"
    } else {
        device <- "cpu:0"
    }
    diag_vec <- torch::torch_tensor(diag_vec)$reshape(c(1, length(diag_vec)))$to(device)
    if (mode == "vec") return(as_array((diag_vec^(-1))$to(device = "cpu:0")))
    return(diag(as_array((diag_vec^(-1))$to(device = "cpu:0"))))
}

#' Title
#'
#' @param vec 
#' @param mode 
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
large_vec_norm_torch <- function(vec, mode = "Gelder", p = 2) {
    require(torch)
    if (torch::cuda_is_available()) {
        device <- "cuda:0"
    } else {
        device <- "cpu:0"
    }
    vec <- torch::torch_tensor(vec)$reshape(c(1, length(vec)))$to(device)
    
    if (mode == "Gelder") return(as_array((sum(abs(vec)^p)^(1/p))$to(device = "cpu:0")))
    return(as_array(max(vec)$to(device = "cpu:0")))
}
