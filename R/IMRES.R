# an iterative minimal residual algorithm -------------------------------------------------------------------

#' An iterative minimal residual algorithm [IMRES]
#' (Итерационный метод минимальных невязок)
#' @description Non-stationary iterative method for solving systems of linear algebraic equations. The essence of the method is to search at each iteration of the iterative parameter by the projection method. This method is a special case of the generalized minimal residual algorithm (GMRES) with the Krylov subspace of the first degree.
#' (Нестационарный итерационный метод решения систем линейных алгебраических уравнений. Суть метода состоит в поиске на каждой итерации итерационного параметра проекционным методом. Данный метод является частным случаем обобщенного метода минимальных невязок с подпространством Крылова первой степени.)
#' @param A - the original matrix of the operator equation - numeric or complex matrix (исходная матрица операторного уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric or complex vector (начальное приближение неизвестного вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when the method diverges (ограничение сверху на число итераций при расхождении метода)
#'
#' @return u - unknown vector in some approximation (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples Reductor <- 10
#' AN <- diag(seq(0.1, 0.9, 0.8/Reductor))   # Creating diag numeric matrix
#' fN <- rnorm(Reductor + 1)
#' uN <- rnorm(Reductor + 1)
#' result <- IMRES(AN, fN, uN, eps = 10e-5)
#' print(result)
#' all.equal(solve(AN) %*% fN, result) 
#' 
#' AC <- diag(seq(0.1, 0.9, 0.8/Reductor) + 1i * seq(0.1, 0.9, 0.8/Reductor))
#' fC <- rnorm(Reductor + 1) + 1i * rnorm(Reductor + 1)
#' uC <- rnorm(Reductor + 1) + 1i * rnorm(Reductor + 1)
#' result <- IMRES(AC, fC, uC, eps = 10e-5)
#' print(result)
#' all.equal(solve(AC) %*% fC, result) 
#' 
#' M <- matrix(rnorm((Reductor + 1) ^ 2), nrow = Reductor + 1, ncol = Reductor + 1)
#' print(IMRES(M, fC, uN, eps = 10e-6, iterations = 10000))
IMRES <- function(A, f, u, eps = 10e-4, iterations = 10000) {
    
    stopifnot(is.matrix(A), 
              is.numeric(A) || is.complex(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, is.atomic(eps), 
              nrow(A) == ncol(A), ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, 
              is.numeric(iterations), length(iterations) == 1, is.atomic(iterations))
    
    h <- A %*% u - f
    if (abs((1 - (abs(t(h) %*%  Conj(A %*% h))^2 / ((t(h) %*% Conj(h)) * (t(A %*% h) %*% Conj(A %*% h)))))[1, 1]) < 0)
        stop("q >= 1, method is growing")
    i <- 0
    repeat {
        h <- A %*% u - f
        tau <- (t(h) %*% Conj(A %*% h)) / (t(A %*% h) %*% Conj(A %*% h))
        u <- u - tau[1,1] * h
        i <- i + 1
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / (sqrt(t(f) %*% Conj(f)))) < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to the final result / allowed number of iterations is exceeded")
            break
        }
    }
    return(u)
}


# IMRES history -----------------------------------------------------------

#' An iteration minimal residual algorithm history [IMRES.history]
#' (Итерационный метод минимальных невязок)
#' @description Non-stationary method for solving systems of linear algebraic equations. The essence of the method is to search at each iteration of the iterative parameter by the projection method. This method is a special case of the generalized minimal residual algorithm (GMRES) with the Krylov subspace of the first degree.
#' (Нестационарный метод решения систем линейных алгебраических уравнений. Суть метода состоит в поиске на каждой итерации итерационного параметра проекционным методом. Данный метод является частным случаем обобщенного метода минимальных невязок с подпространством Крылова первой степени.)
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
#' @examples Reductor <- 10
#' AN <- diag(seq(0.1, 0.9, 0.8/Reductor))   # Creating diag numeric matrix
#' fN <- rnorm(Reductor + 1)
#' uN <- rnorm(Reductor + 1)
#' result <- IMRES.history(AN, fN, uN, eps = 10e-5)
#' print(result)
#' all.equal(solve(AN) %*% fN, result$var) 
#' 
#' AC <- diag(seq(0.1, 0.9, 0.8/Reductor) + 1i * seq(0.1, 0.9, 0.8/Reductor))
#' fC <- rnorm(Reductor + 1) + 1i * rnorm(Reductor + 1)
#' uC <- rnorm(Reductor + 1) + 1i * rnorm(Reductor + 1)
#' result <- IMRES.history(AC, fC, uC, eps = 10e-5)
#' print(result)
#' all.equal(solve(AC) %*% fC, result$var) 
#' M <- matrix(rnorm((Reductor + 1) ^ 2), nrow = Reductor + 1, ncol = Reductor + 1)
#' print(IMRES.history(M, fC, uN, eps = 10e-6, iterations = 10000))
IMRES.history <- function(A, f, u, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A) || is.complex(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, is.atomic(eps), 
              nrow(A) == ncol(A), ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, 
              is.numeric(iterations), length(iterations) == 1, is.atomic(iterations))
    iterate <- 0
    iterate2 <- 0
    h <- A %*% u - f
    if (abs((1 - (abs(t(h) %*%  Conj(A %*% h))^2 / ((t(h) %*% Conj(h)) * (t(A %*% h) %*% Conj(A %*% h)))))[1, 1]) < 0)
        stop("q >= 1, method is growing")
    i <- 0
    u.hist <- matrix(u, nrow = dim(A)[1])
    t1 <- Sys.time()
    repeat {
        h <- A %*% u - f
        tau <- (t(h) %*% Conj(A %*% h)) / (t(A %*% h) %*% Conj(A %*% h))
        u <- u - tau[1,1] * h
        i <- i + 1
        iterate <- iterate + 4
        iterate2 <- c(iterate2, iterate)
        u.hist <- cbind(u.hist, u)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / (sqrt(t(f) %*% Conj(f)))) < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to the final result / allowed number of iterations is exceeded")
            break
        }
    }
    t2 <- Sys.time()
    return(list(num.iter = iterate2, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}
