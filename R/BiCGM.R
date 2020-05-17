# Biconjugate gradient method, BiCGM --------------------------------------

# Данный алгоритм предназначен только для действительных матриц

#' Biconjugate gradient method [BiCGM]
#' (Метод бисопряженных градиентов)
#' @description Non-stationary iterative numerical method 
#' for solving SLAEs of the Krylov type. It is a generalization 
#' of the conjugate gradient method.
#' (Нестационарный итерационный численный метод решения 
#' СЛАУ крыловского типа. Является обобщением метода 
#' сопряжённых градиентов.)
#' @param A - the original matrix of the operator equation 
#' - numeric matrix only (исходная матрица операторного уравнения 
#' - вещественная только)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение 
#' неизвестного вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - 
#' numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations 
#' when the method diverges (ограничение сверху на число итераций 
#' при расхождении метода)
#'
#' @return u - unknown vector in some approximation 
#' (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples 
#' A <- diag(rnorm(5, 2), nrow = 5, ncol = 5)
#' u <- rnorm(5, 12)
#' f <- rnorm(5, 17)
#' solve(A) %*% f - BiCGM(A, f, u, iterations = 10000)
#' all.equal(solve(A) %*% f, BiCGM(A, f, u, iterations = 10000))
BiCGM <- function(A, f, u, eps = 10e-04, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps), nrow(A) == ncol(A), 
              ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, is.numeric(iterations), 
              length(iterations) == 1, is.atomic(iterations))
    r <- f - A %*% u
    p <- r
    z <- r
    s <- r
    i <- 0
    repeat {
        alpha <- ((t(p) %*% Conj(r)) / 
                      (t(s) %*% Conj(A %*% z)))[1, 1]
        u <- u + alpha * z
        i <- i + 1
        r1 <- r - alpha * (A %*% z)
        p1 <- p - alpha * (t(A) %*% s)
        beta <- ((t(p1) %*% Conj(r1)) / 
                     (t(p) %*% Conj(r)))[1, 1]
        z <- r1 + beta * z
        s <- p1 + beta * s
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / 
                (sqrt(t(f) %*% Conj(f)))) < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come 
                    close to the final result / allowed number 
                    of iterations is exceeded")
            break
        }
        r <- r1
        p <- p1
        rm(r1)
        rm(p1)
    }
    return(u)
}

# BiCGM.history -----------------------------------------------------------

#' Biconjugate gradient method [BiCGM]
#' (Метод бисопряженных градиентов)
#' @description Non-stationary iterative numerical method 
#' for solving SLAEs of the Krylov type. It is a 
#' generalization of the conjugate gradient method.
#' (Нестационарный итерационный численный метод решения 
#' СЛАУ крыловского типа. Является обобщением метода 
#' сопряжённых градиентов.)
#' @details This method is necessary to preserve the history 
#' of sequential calculation of an unknown vector in order 
#' to visualize the convergence of the method 
#' (Данный метод необходим для сохранения истории 
#' последовательного вычисления неизвестного вектора с целью 
#' визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - 
#' numeric matrix only (исходная матрица операторного уравнения 
#' - вещественная только)
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
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного 
#' вектора); 
#' var.hist - history of computing an unknown vector 
#' (история вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время 
#' вычисления); 
#' @export
#'
#' @examples 
#' A <- diag(rnorm(5, 2), nrow = 5, ncol = 5)
#' u <- rnorm(5, 12)
#' f <- rnorm(5, 17)
#' print(BiCGM.history(A, f, u, iterations = 10000))
BiCGM.history <- function(A, f, u, eps = 10e-04, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps), nrow(A) == ncol(A), 
              ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, is.numeric(iterations), 
              length(iterations) == 1, is.atomic(iterations))
    iterate <- 0
    iterate2 <- 0
    i <- 0
    u.hist <- matrix(u, nrow = nrow(A))
    t1 <- Sys.time()
    r <- f - A %*% u
    iterate <- iterate + 1
    p <- r
    z <- r
    s <- r
    repeat {
        alpha <- ((t(p) %*% Conj(r)) / 
                      (t(s) %*% Conj(A %*% z)))[1, 1]
        u <- u + alpha * z
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        r1 <- r - alpha * (A %*% z)
        p1 <- p - alpha * (t(A) %*% s)
        beta <- ((t(p1) %*% Conj(r1)) / 
                     (t(p) %*% Conj(r)))[1, 1]
        z <- r1 + beta * z
        s <- p1 + beta * s
        iterate <- iterate + 3
        iterate2 <- c(iterate2, iterate)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / 
                (sqrt(t(f) %*% Conj(f)))) < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come 
                    close to the final result / allowed number 
                    of iterations is exceeded")
            break
        }
        r <- r1
        p <- p1
        rm(r1)
        rm(p1)
    }
    t2 <- Sys.time()
    return(list(num.iter = iterate2, var = u, var.hist = u.hist, 
                systime.iter = difftime(t2, t1, units = "secs")[[1]]))

}
