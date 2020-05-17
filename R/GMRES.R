# GMRES -------------------------------------------------------------------

#' Generalized minimal residual algorithm [GMRES]
#' (Обобщенный метод минимальных невязок)
#' @description Non-stationary iterative numerical method for solving 
#' systems of linear algebraic equations. The projection of the vector 
#' onto the Krylov subspace of arbitrary order is used as a residual.
#' (Нестационарный итерационный численный метод решения систем 
#' линейных алгебраических уравнений. В качестве невязки используется 
#' проекция вектора на подпространство Крылова произвольного порядка.)
#' @param A - the original matrix of the operator equation - numeric 
#' or complex matrix (исходная матрица операторного уравнения - 
#' вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных членов 
#' вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric 
#' or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param layers - Krylov subspace order (порядок подпространства 
#' Крылова)
#' @param eps - accuracy of calculation of the desired vector - numeric 
#' (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return u - unknown vector in some approximation 
#' (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples A <- diag(rnorm(25, 5) + 1i * rnorm(25, 1), ncol = 25, nrow = 25)
#' f <- rnorm(25, 2)
#' u <- rnorm(25) 
#' print(GMRES(A, f, u, layers = 5, eps = 10e-07))
#' print(IMRES(A, f, u, eps = 10e-7))
GMRES <- function(A, f, u, layers = 2, eps = 10e-04, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps),
              nrow(A) == ncol(A), ncol(A) == length(f), 
              length(f) == length(u),
              ncol(A) >= 2,
              is.numeric(iterations), length(iterations) == 1, 
              is.atomic(iterations),
              is.atomic(layers), layers >= 1, layers <= 1000)
    BetaLK <- function(h, nu) {
        return(
            ((t(A %*% h) %*% Conj(A %*% nu))
             / (t(A %*% nu) %*% Conj(A %*% nu)))[1, 1]
        )
    }
    sumBetaLK <- function(l, h, nu) {
        cumsumR <- vector(length = length(h))
        for (i in 1:(l - 1)) {
            cumsumR <- cumsumR + BetaLK(h, nu[, i]) * nu[, i]
        }
        return(as.vector(cumsumR))
    }
    i <- 0      # iteration num
    repeat {
        h <- matrix(0, ncol = layers + 1, nrow = length(u))
        h[, 1] <- A %*% u - f
        nu <- matrix(0, ncol = layers, nrow = length(u))
        for (l in 1:layers) {
            if (l == 1) {
                nu[, l] <- h[, l]
            } else {
                nu[, l] <- as.vector(h[, l]) - sumBetaLK(l, as.vector(h[, l]), nu)
            }
            tau <- ((t(as.vector(h[, l])) %*% 
                         Conj(A %*% as.vector(h[, l]))) / 
                        (t(A %*% as.vector(nu[, l])) %*% 
                             Conj(A %*% as.vector(nu[, l]))))[1,1]
            u <- u - tau * as.vector(nu[, l])
            h[, (l + 1)] <- h[, l] - tau * (A %*% as.vector(nu[, l]))
        }
        i <- i + 1
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    return(u)

}

# GMRES.history -----------------------------------------------------------

#' Generalized minimal residual algorithm history [GMRES.history]
#' (Обобщенный метод минимальных невязок)
#' @description Non-stationary iterative numerical method for solving 
#' systems of linear algebraic equations. The projection of the vector 
#' onto the Krylov subspace of arbitrary order is used as a residual.
#' (Нестационарный итерационный численный метод решения систем 
#' линейных алгебраических уравнений. В качестве невязки используется 
#' проекция вектора на подпространство Крылова произвольного порядка.)
#' @details This method is necessary to preserve the history of sequential 
#' calculation of an unknown vector in order to visualize the 
#' convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного 
#' вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - numeric 
#' or complex matrix (исходная матрица операторного уравнения 
#' - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных 
#' членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param layers - Krylov subspace order (порядок подпространства 
#' Крылова)
#' @param eps - accuracy of calculation of the desired vector 
#' - numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations 
#' when the method diverges (ограничение сверху на число итераций 
#' при расхождении метода)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история вычисления 
#' неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples A <- diag(rnorm(25, 5) + 1i * rnorm(25, 1), ncol = 25, nrow = 25)
#' f <- rnorm(25, 2)
#' u <- rnorm(25) 
#' print(GMRES.history(A, f, u, layers = 5, eps = 10e-07))
#' print(IMRES.history(A, f, u, eps = 10e-7))
GMRES.history <- function(A, f, u, layers = 2, eps = 10e-04, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps),
              nrow(A) == ncol(A), ncol(A) == length(f), 
              length(f) == length(u),
              ncol(A) >= 2,
              is.numeric(iterations), length(iterations) == 1, 
              is.atomic(iterations),
              is.atomic(layers), layers >= 1, layers <= 1000)
    BetaLK <- function(h, nu) {
        return(
            ((t(A %*% h) %*% Conj(A %*% nu))
             / (t(A %*% nu) %*% Conj(A %*% nu)))[1, 1]
        )
    }
    sumBetaLK <- function(l, h, nu) {
        cumsumR <- vector(length = length(h))
        for (i in 1:(l - 1)) {
            cumsumR <- cumsumR + BetaLK(h, nu[, i]) * nu[, i]
        }
        return(as.vector(cumsumR))
    }
    LayerIterations <- function(x) {
        x1 <- 3
        if (x >= 3) {
            for (i in 3:x) {
                x2 <- x1 + 3 * (i - 1)
                x1 <- x2
            }
        }
        return(ifelse(x == 1, 0, ifelse(x == 2, 3, x2)))
    }
    t1 <- Sys.time()
    i <- 0      # iteration num
    iterate <- 0
    iterate2 <- 0
    u.hist <- matrix(u, nrow = nrow(A))
    repeat {
        h <- matrix(0, ncol = layers + 1, nrow = length(u))
        h[, 1] <- A %*% u - f
        nu <- matrix(0, ncol = layers, nrow = length(u))
        for (l in 1:layers) {
            if (l == 1) {
                nu[, l] <- h[, l]
            } else {
                nu[, l] <- as.vector(h[, l]) - 
                    sumBetaLK(l, as.vector(h[, l]), nu)
            }
            tau <- ((t(as.vector(h[, l])) %*% 
                         Conj(A %*% as.vector(h[, l]))) / 
                        (t(A %*% as.vector(nu[, l])) %*% 
                             Conj(A %*% as.vector(nu[, l]))))[1,1]
            u <- u - tau * as.vector(nu[, l])
            h[, (l + 1)] <- h[, l] - tau * (A %*% as.vector(nu[, l]))
        }
        iterate <- iterate + 4 + LayerIterations(layers)
        iterate2 <- c(iterate2, iterate)
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    t2 <- Sys.time()
    return(list(num.iter = iterate2, var = u, var.hist = u.hist, 
                systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}
