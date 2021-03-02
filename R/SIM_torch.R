
#' Simple Iterations Method {torch}
#'
#' @param A 
#' @param u 
#' @param f 
#' @param eps 
#' @param iterations 
#' @import torch
#'
#' @return
#' @export
#'
#' @examples
SIM_torch <- function(A, u, f, eps = 10e-4, iterations = 1000) {
#     stopifnot(is.matrix(A) || is.vector(A), 
#               is.numeric(A) || is.complex(A), 
#               is.numeric(f) || is.complex(f), 
#               is.numeric(u) || is.complex(u), 
#               is.numeric(eps), length(eps) == 1, 
#               is.atomic(eps), nrow(A) == ncol(A), 
#               ncol(A) == length(f), length(f) == length(u), 
#               ncol(A) >= 2, 
#               is.numeric(iterations), 
#               length(iterations) == 1, is.atomic(iterations))
#     require(torch)
#     if (torch::cuda_is_available()) {
#         device = "cuda:0"
#     } else {
#         device = "cpu:0"
#     }
#     
#     u_torch <- torch::torch_tensor(u)$reshape(c(length(u), 1))$to(device = device)
#     f_torch <- torch::torch_tensor(f)$reshape(c(length(f), 1))$to(device = device)
#     
#     if (is.vector(A)) {
#         A_torch <- torch::torch_tensor(A)$reshape(c(length(A), 1))$to(device = device)    
#         B_torch <- 1 - A_torch
#         i <- 0
#         repeat {
#             u_torch <- B_torch * u_torch + f_torch
#             i <- i + 1
#             if (abs())
#             if (i > iterations) {
#                 message("Iterations of the method may
#                             not come close to the final result
#                             / allowed number of iterations is
#                             exceeded / check the spectrum of
#                             operator A: sigma(A) must be less than 1")
#                 break
#         }
#                 
#         # dimA <- dim(A)[1]
#         # B <- diag(1, nrow = dimA, ncol = dimA) - A
#         # i <- 0
#         # repeat {
#         #     u <- B %*% u + f
#         #     i <- i + 1
#         #     if (abs((sqrt(t(A %*% u - f) 
#         #                   %*% Conj(A %*% u - f))) / 
#         #             (sqrt(t(f) %*% Conj(f)))) < eps) break
#         #     if (i > iterations) {
#         #         message("Iterations of the method may 
#         #             not come close to the final result 
#         #             / allowed number of iterations is 
#         #             exceeded / check the spectrum of 
#         #             operator A: sigma(A) must be less than 1")
#         #         break
#         #     }
#         # }
#         # return(u)
#         
#     } else if (is.matrix(A)) {
#         A <- torch::torch_tensor(A)$reshape(c(nrow(A), ncol(A)))$to(device = device)
#         
#     }
#     
}
