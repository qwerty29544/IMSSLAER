# 
# # Biconjugate gradient stabilized method, BiCGMStab -----------------------
# 
# # https://ru.wikipedia.org/wiki/%D0%A1%D1%82%D0%B0%D0%B1%D0%B8%D0%BB%D0%B8%D0%B7%D0%B8%D1%80%D0%BE%D0%B2%D0%B0%D0%BD%D0%BD%D1%8B%D0%B9_%D0%BC%D0%B5%D1%82%D0%BE%D0%B4_%D0%B1%D0%B8%D1%81%D0%BE%D0%BF%D1%80%D1%8F%D0%B6%D1%91%D0%BD%D0%BD%D1%8B%D1%85_%D0%B3%D1%80%D0%B0%D0%B4%D0%B8%D0%B5%D0%BD%D1%82%D0%BE%D0%B2
# 
# 
# # Устойчивый алгоритм, итерации проверять не нужно
# BiCGMStab <- function(A, f, u, eps = 0.00001) {
#     stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
#     
#     dimA <- dim(A)[1]
#     
#     # Проверка на n >= 2
#     if (dimA[1] < 2) stop("Operator must have dim >= 2")
#     
#     # Проверка на размерность матрицы опер`атора
#     if (dim(A)[1] != dim(A)[2]) stop("Operator must be quadratic")
#     
#     r <- f - A %*% u
#     rtild <- r
#     rho <- alpha <- omega <- 1
#     v <- p <- 0
#     
#     repeat {
#         rho1 <- (t(rtild) %*% r)[1, 1]
#         beta <- (rho1 / rho) * (alpha / omega)
#         p <- r + beta * (p - omega * v)
#         v <- A %*% p
#         alpha <- (rho1 / (t(rtild) %*% v)[1, 1])
#         s <- rho - alpha * v
#         t1 <- (A %*% s)
#         omega <- ((Conj(t(t1)) %*% s) / (Conj(t(t1)) %*% t1))[1, 1]
#         u <- u + omega * s + alpha * p
#         r <- s - omega * t1
#         if (sqrt(t(r) %*% r)/sqrt(t(f) %*% f) < eps) break
#         rho <- rho1
#     }
#     
#     return(u)
# }
# 
# 
# # Testing -----------------------------------------------------------------
# 
# 
# A <- diag(rnorm(25), nrow = 5, ncol = 5)
# f <- rnorm(5)
# u <- rnorm(5)
# 
# BiCGMStab(A, f, u)
