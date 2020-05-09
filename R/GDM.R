# GDM ---------------------------------------------------------------------

#' Gradient Descent method [GDM]
#' (метод градиентного спуска)
#' @description The non-stationary iterative gradient descent method defined for both active and complex matrices. The cornerstone is the search for a minimum of a parabolic functional. A minimum of functionality will correspond to an unknown vector.
#' (Нестационарный итерационный метод градиентного спуска, определённый как для действительных, так и для комплексных матриц. В основе метода лежит поиск минимума параболического функционала. Минимум функционала будет соответствовать неизвестному вектору.)
#' @param A - the original matrix of the operator equation - numeric or complex matrix (исходная матрица операторного уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric or complex vector (начальное приближение неизвестного вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when the method diverges (ограничение сверху на число итераций при расхождении метода)
#'
#' @return u - unknown vector in some approximation (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples A <- matrix(rnorm(100) + 1i * rnorm(100), ncol = 10, nrow = 10)
#' f <- rnorm(10) + 1i * rnorm(10)
#' u <- rnorm(10)
#' print(IMSSLAER::GDM(A, f, u))
GDM <- function(A, f, u, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A) || is.complex(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, is.atomic(eps), 
              nrow(A) == ncol(A), ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, 
              is.numeric(iterations), length(iterations) == 1, is.atomic(iterations))
    if (is.numeric(A)) {
        i <- 0
        repeat {
            r <- A %*% u - f
            u <- u - ((t(t(A) %*% r) %*% (t(A) %*% r))/(t(A %*% t(A) %*% r) %*% (A %*% t(A) %*% r)))[1,1] * (t(A) %*% r)
            i <- i + 1
            if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / (sqrt(t(f) %*% Conj(f))))[1,1] < eps) break
            if (i > iterations) {
                message("Iterations of the method may not come close to the final result / allowed number of iterations is exceeded")
                break
            }
        }
        return(u)
    } else if (is.complex(A)) {
        HA <- Re(A)
        HB <- Im(A)
        Hv <- A
        Hb <- t(HA) - 1i * t(HB)
        rm(HA)
        rm(HB)
        i <- 0
        repeat {
            r <- Hv %*% u - f
            u <- u - ((t(Hb %*% r) %*% Conj(Hb %*% r)) / (t(Hv %*% Hb %*% r) %*% Conj(Hv %*% Hb %*% r)))[1, 1] * 
                (Hb %*% r)
            i <- i + 1
            if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
            if (i > iterations) {
                message("Iterations of the method may not come close to the final result / allowed number of iterations is exceeded")
                break
            }
        }
        return(u)
    }
}


# GDM история -------------------------------------------------------------

#' Gradient Descent method [GDM]
#' (метод градиентного спуска)
#' @description The non-stationary iterative gradient descent method defined for both active and complex matrices. The cornerstone is the search for a minimum of a parabolic functional. A minimum of functionality will correspond to an unknown vector.
#' (Нестационарный итерационный метод градиентного спуска, определённый как для действительных, так и для комплексных матриц. В основе метода лежит поиск минимума параболического функционала. Минимум функционала будет соответствовать неизвестному вектору.)
#' @details This method is necessary to preserve the history of sequential calculation of an unknown vector in order to visualize the convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - numeric or complex matrix (исходная матрица операторного уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric or complex vector (начальное приближение неизвестного вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when the method diverges (ограничение сверху на число итераций при расхождении метода) 
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples A <- matrix(rnorm(100) + 1i * rnorm(100), ncol = 10, nrow = 10)
#' f <- rnorm(10) + 1i * rnorm(10)
#' u <- rnorm(10)
#' print(IMSSLAER::GDM.history(A, f, u))
GDM.history <- function(A, f, u, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A) || is.complex(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, is.atomic(eps), 
              nrow(A) == ncol(A), ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, 
              is.numeric(iterations), length(iterations) == 1, is.atomic(iterations))
    if (is.numeric(A)) {
        i <- 0
        iterate <- 0
        iterate2 <- 0
        u.hist <- matrix(u, nrow = nrow(A))
        t1 <- Sys.time()
        repeat {
            r <- A %*% u - f
            u <- u - ((t(t(A) %*% r) %*% (t(A) %*% r))/(t(A %*% t(A) %*% r) %*% (A %*% t(A) %*% r)))[1,1] * (t(A) %*% r)
            u.hist <- cbind(u.hist, u)
            i <- i + 1
            iterate <- iterate + 8
            iterate2 <- c(iterate2, iterate)
            if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / (sqrt(t(f) %*% Conj(f))))[1,1] < eps) break
            if (i > iterations) {
                message("Iterations of the method may not come close to the final result / allowed number of iterations is exceeded")
                break
            }
        }
        t2 <- Sys.time()
        return(list(num.iter = iterate2, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
    } else if (is.complex(A)) {
        HA <- Re(A)
        HB <- Im(A)
        Hv <- A
        Hb <- t(HA) - 1i * t(HB)
        rm(HA)
        rm(HB)
        i <- 0
        iterate <- 0
        iterate2 <- 0
        u.hist <- matrix(u, nrow = nrow(A))
        t1 <- Sys.time()
        repeat {
            r <- Hv %*% u - f
            u <- u - ((t(Hb %*% r) %*% Conj(Hb %*% r)) / 
                          (t(Hv %*% Hb %*% r) %*% 
                               Conj(Hv %*% Hb %*% r)))[1, 1] * (Hb %*% r)
            u.hist <- cbind(u.hist, u)
            i <- i + 1
            iterate <- iterate + 8
            iterate2 <- c(iterate2, iterate)
            if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
            if (i > iterations) {
                message("Iterations of the method may not come close to the final result / allowed number of iterations is exceeded")
                break
            }
        }
        t2 <- Sys.time()
        return(list(num.iter = iterate2, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
    }
}
