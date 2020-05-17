# Conjugate gradient method -----------------------------------------------

#' Conjugate gradient method [CGM]
#' (Метод сопряженных градиентов)
#' @description A numerical method for solving systems 
#' of linear algebraic equations is an unsteady iterative 
#' method of the Krylov type. 
#' (Численный метод решения систем линейных алгебраических 
#' уравнений, является нестационарным итерационным методом 
#' Крыловского типа.)
#' @param A - the original matrix of the operator equation 
#' - numeric or complex matrix (исходная матрица операторного 
#' уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных 
#' членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного 
#' вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector 
#' - numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of 
#' iterations when the method diverges (ограничение сверху 
#' на число итераций при расхождении метода)
#'
#' @return u - unknown vector in some approximation 
#' (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples 
#' A <- diag(rnorm(5, 20) + 1i * rnorm(5, 9), nrow = 5, ncol = 5)
#' u <- rnorm(5, 12)
#' f <- rnorm(5, 17)
#' solve(A) %*% f - CGM(A, f, u, iterations = 10000)
#' all.equal(solve(A) %*% f, CGM(A, f, u, iterations = 10000))
CGM <- function(A, f, u, eps = 10e-04, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A) || is.complex(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps), nrow(A) == ncol(A), 
              ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, is.numeric(iterations), 
              length(iterations) == 1, is.atomic(iterations))
    i <- 0
    r <- f - A %*% u
    z <- r
    repeat {
        alpha = (t(r) %*% Conj(r)) / (t(A %*% z) %*% Conj(z))
        u = u + alpha[1,1] * z
        r1 = r - alpha[1,1] * (A %*% z)
        beta <- (t(r1) %*% Conj(r1)) / (t(r) %*% Conj(r))
        z = r1 + beta[1,1] * z
        i <- i + 1
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / 
                (sqrt(t(f) %*% Conj(f)))) < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come 
                    close to the final result / allowed 
                    number of iterations is exceeded")
            break
        }
        r = r1
        rm(r1)
    }
    return(u)
}


# CGM.history -------------------------------------------------------------

#' Conjugate gradient method history [CGM.history]
#' (Метод сопряженных градиентов)
#' @description A numerical method for solving systems 
#' of linear algebraic equations is an unsteady iterative 
#' method of the Krylov type. 
#' (Численный метод решения систем линейных алгебраических 
#' уравнений, является нестационарным итерационным методом 
#' Крыловского типа.)
#' @details This method is necessary to preserve the 
#' history of sequential calculation of an unknown 
#' vector in order to visualize the convergence of 
#' the method 
#' (Данный метод необходим для сохранения истории 
#' последовательного вычисления неизвестного вектора 
#' с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation 
#' - numeric or complex matrix (исходная матрица операторного 
#' уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение 
#' неизвестного вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector 
#' - numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of 
#' iterations when the method diverges (ограничение сверху 
#' на число итераций при расхождении метода)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного 
#' вектора); 
#' var.hist - history of computing an unknown vector (история 
#' вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время 
#' вычисления); 
#' @export
#'
#' @examples 
#' A <- diag(rnorm(5, 20) + 1i * rnorm(5, 9), nrow = 5, ncol = 5)
#' u <- rnorm(5, 12)
#' f <- rnorm(5, 17)
#' CGM.history(A, f, u, eps = 0.00001, iterations = 10000)
CGM.history <- function(A, f, u, eps = 10e-04, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A) || is.complex(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps), nrow(A) == ncol(A), 
              ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, is.numeric(iterations), 
              length(iterations) == 1, is.atomic(iterations))
    i <- 0
    iterate <- 0
    iterate2 <- 0
    u.hist <- matrix(u, nrow = nrow(A))
    t1 <- Sys.time()

    r <- f - A %*% u
    z <- r
    iterate <- iterate + 1
    repeat {
        alpha = (t(r) %*% Conj(r)) / (t(A %*% z) %*% Conj(z))
        u = u + alpha[1,1] * z
        u.hist <- cbind(u.hist, u)
        r1 = r - alpha[1,1] * (A %*% z)
        beta <- (t(r1) %*% Conj(r1)) / (t(r) %*% Conj(r))
        z = r1 + beta[1,1] * z
        i <- i + 1
        iterate <- iterate + 2
        iterate2 <- c(iterate2, iterate)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / 
                (sqrt(t(f) %*% Conj(f)))) < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come 
                    close to the final result / allowed number 
                    of iterations is exceeded")
            break
        }
        r = r1
        rm(r1)
    }
    t2 <- Sys.time()
    return(list(num.iter = iterate2, var = u, var.hist = u.hist, 
                systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}
