# Biconjugate gradient stabilized method, BiCGMStab -----------------------

#' Biconjugate gradient method stabilized [BiCGMStab]
#' (Метод бисопряженных градиентов стабилизированный)
#' @description Iterative method for solving SLAE of the 
#' Krylov type. Developed by van der Vorst to solve systems 
#' with asymmetric matrices. It converges faster than the usual 
#' method of bis conjugate gradients, which is unstable, and 
#' therefore is used more often.
#' (Нестационарный итерационный метод решения СЛАУ крыловского типа. 
#' Разработан Ван дэр Ворстом для решения систем с несимметричными 
#' матрицами. Сходится быстрее, чем обычный метод бисопряженных 
#' градиентов, который является неустойчивым, и поэтому применяется чаще.)
#' @param A - the original matrix of the operator equation - numeric 
#' or complex matrix (исходная матрица операторного уравнения - 
#' вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных 
#' членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric 
#' or complex vector (начальное приближение неизвестного вектора - 
#' вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - 
#' numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return u - unknown vector in some approximation (неизвестный 
#' вектор в некотором приближении)
#' @export
#'
#' @examples 
#' A <- diag(rnorm(5, 2), nrow = 5, ncol = 5)
#' u <- rnorm(5, 12)
#' f <- rnorm(5, 17)
#' solve(A) %*% f - BiCGMStab(A, f, u, iterations = 10000)
#' all.equal(solve(A) %*% f, BiCGMStab(A, f, u, iterations = 10000))
BiCGMStab <- function(A, f, u, eps = 10e-04, iterations = 10000) {
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
    rtild <- r
    rho <- alpha <- omega <- 1
    v <- p <- 0
    repeat {
        rho1 <- (Conj(t(rtild)) %*% r)[1, 1]
        beta <- (rho1 / rho) * (alpha / omega)
        p <- r + beta * (p - omega * v)
        v <- A %*% p
        alpha <- (rho1 / (Conj(t(rtild)) %*% v)[1, 1])
        s <- r - alpha * v
        t1 <- (A %*% s)
        omega <- ((Conj(t(t1)) %*% s) / 
                      (Conj(t(t1)) %*% t1))[1, 1]
        u <- u + omega * s + alpha * p
        r <- s - omega * t1
        i <- i + 1
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / 
                (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close 
                    to the final result / allowed number of 
                    iterations is exceeded")
            break
        }
        rho <- rho1
    }

    return(u)
}


# BiCGMStab history -------------------------------------------------------

#' Biconjugate gradient method stabilized history [BiCGMStab]
#' (Метод бисопряженных градиентов стабилизированный)
#' @description Iterative method for solving SLAE of the Krylov 
#' type. Developed by van der Vorst to solve systems with 
#' asymmetric matrices. It converges faster than the usual 
#' method of bis conjugate gradients, which is unstable, and 
#' therefore is used more often.
#' (Нестационарный итерационный метод решения СЛАУ крыловского 
#' типа. Разработан Ван дэр Ворстом для решения систем с 
#' несимметричными матрицами. Сходится быстрее, чем обычный метод 
#' бисопряженных градиентов, который является неустойчивым, и 
#' поэтому применяется чаще.)
#' @details This method is necessary to preserve the history of 
#' sequential calculation of an unknown vector in order to 
#' visualize the convergence of the method 
#' (Данный метод необходим для сохранения истории 
#' последовательного вычисления неизвестного вектора с целью 
#' визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - 
#' numeric or complex matrix (исходная матрица операторного 
#' уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных 
#' членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного 
#' вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - 
#' numeric (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations 
#' when the method diverges (ограничение сверху на число итераций 
#' при расхождении метода)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного 
#' вектора); 
#' var.hist - history of computing an unknown vector (история 
#' вычисления неизвестного вектора); 
#' @export
#'
#' @examples 
#' A <- diag(rnorm(5, 2), nrow = 5, ncol = 5)
#' u <- rnorm(5, 12)
#' f <- rnorm(5, 17)
#' print(BiCGMStab.history(A, f, u, iterations = 10000))
BiCGMStab.history <- function(A, f, u, eps = 10e-04, iterations = 10000) {
    stopifnot(is.matrix(A), 
              is.numeric(A) || is.complex(A), 
              is.numeric(f) || is.complex(f), 
              is.numeric(u) || is.complex(u), 
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps), nrow(A) == ncol(A), 
              ncol(A) == length(f), length(f) == length(u), 
              ncol(A) >= 2, is.numeric(iterations), 
              length(iterations) == 1, is.atomic(iterations))
    u.hist <- matrix(u, nrow = nrow(A))
    i <- 0
    iterate <- 0
    iterate2 <- 0
    r <- f - A %*% u
    rtild <- r
    rho <- alpha <- omega <- 1
    v <- p <- 0
    repeat {
        rho1 <- (Conj(t(rtild)) %*% r)[1, 1]
        beta <- (rho1 / rho) * (alpha / omega)
        p <- r + beta * (p - omega * v)
        v <- A %*% p
        alpha <- (rho1 / (Conj(t(rtild)) %*% v)[1, 1])
        s <- r - alpha * v
        t1 <- (A %*% s)
        omega <- ((Conj(t(t1)) %*% s) / 
                      (Conj(t(t1)) %*% t1))[1, 1]
        u <- u + omega * s + alpha * p
        u.hist <- cbind(u.hist, u)
        r <- s - omega * t1
        i <- i + 1
        iterate <- iterate + 2
        iterate2 <- c(iterate2, iterate)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) / 
                (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come 
                    close to the final result / allowed number 
                    of iterations is exceeded")
            break
        }
        rho <- rho1
    }
    return(list(num.iter = iterate2, var = u, var.hist = u.hist))
    
}
