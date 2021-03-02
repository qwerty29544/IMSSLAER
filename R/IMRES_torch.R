require(torch)

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
#' @export
#' @import torch
#'
#' @examples
IMRES_torch <- function(A, f, u, eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    if (dim(A)[1] < 2) stop("Linear operator has dim less than 2x2")
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be a quadratic")
    if (cuda_is_available()) device = "cuda:0"
    else device = "cpu:0"
    A_torch <- torch::torch_tensor(A)$to(device = device)
    u_torch <- torch::torch_tensor(u)$reshape(c(ncol(A), 1))$to(device = device)
    f_torch <- torch::torch_tensor(f)$reshape(c(ncol(A), 1))$to(device = device)
    h <- A_torch$matmul(u_torch) - f_torch
    repeat {
        ut <- u_torch
        h <- A_torch$matmul(u_torch) - f_torch
        tau <- (h$t()$matmul(A_torch$matmul(h))) / (A_torch$matmul(h)$t()$matmul(A_torch$matmul(h)))
        u_torch <- u_torch - tau[1,1] * h
        if(as_array((u_torch - ut)$abs()$max()$to(device = "cpu:0"))[1] < eps) break
    }
    return(as.matrix(as_array(u_torch$to(device = "cpu:0"))))
}
