# IMRES -------------------------------------------------------------------

#' Title
#'
#' @param A - the original matrix of the operator equation (исходная матрица операторного уравнения)
#' @param f - bias (вектор свободных членов)
#' @param u - initial approximation of an unknown vector (начальное приближение неизвестного вектора)
#' @param eps - accuracy of calculation of the desired vector (точность вычисления искомого вектора)
#'
#' @return u - unknown vector in some approximation (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples
IMRES <- function(A, f, u, eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    if (dim(A)[1] < 2) stop("Linear operator has dim less than 2x2")
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be a quadratic")
    h <- A %*% u - f
    if ((1 - ((t(h) %*%  (A %*% h))^2 / ((t(h) %*% h) * (t(A %*% h) %*% (A %*% h))))) < 0)
        stop("q >= 1, method is growing")
    repeat {
        ut <- u
        h <- A %*% u - f
        tau <- (t(h) %*% (A %*% h)) / (t(A %*% h) %*% (A %*% h))
        u <- u - tau[1,1] * h
        if (max(abs(u - ut)) < eps) break
    }
    return(u)
}


# IMRES история -----------------------------------------------------------

#' Title
#' 
#' @details This method is necessary to preserve the history of sequential calculation of an unknown vector in order to visualize the convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation (исходная матрица операторного уравнения)
#' @param f - bias (вектор свободных членов)
#' @param u - initial approximation of an unknown vector (начальное приближение неизвестного вектора)
#' @param eps - accuracy of calculation of the desired vector (точность вычисления искомого вектора)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples
IMRES.history <- function(A, f, u, eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    if (dim(A)[1] < 2) stop("Linear operator has dim less than 2x2")
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be a quadratic")
    h <- A %*% u - f
    if ((1 - ((t(h) %*%  (A %*% h))^2 / ((t(h) %*% h) * (t(A %*% h) %*% (A %*% h))))) < 0)
        stop("q >= 1, method is growing")
    i <- 0
    u.hist <- matrix(u, nrow = dim(A)[1])
    t1 <- Sys.time()
    repeat {
        ut <- u
        h <- A %*% u - f
        tau <- (t(h) %*% (A %*% h)) / (t(A %*% h) %*% (A %*% h))
        u <- u - tau[1,1] * h
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        if (max(abs(u - ut)) < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}