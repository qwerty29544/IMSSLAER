# GDM ---------------------------------------------------------------------

#' Title
#'
#' @param A - the original matrix of the operator equation (исходная матрица операторного уравнения)
#' @param f - bias (вектор свободных членов)
#' @param u - initial approximation of an unknown vector (начальное приближение неизвестного вектора)
#' @param eps - accuracy of calculation of the desired vector (точность вычисления искомого вектора)
#'
#' @return u - unknown vector in some approximation (неизвестный вектор в некотором приближении)
#'
#' @return
#' @export
#'
#' @examples
GDM <- function(A, f, u, eps = 10e-4) {
    stopifnot(is.numeric(A), is.matrix(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    
    repeat {
        r <- A %*% u - f
        u <- u - ((t(t(A) %*% r) %*% (t(A) %*% r))/(t(A %*% t(A) %*% r) %*% (A %*% t(A) %*% r)))[1,1] * (t(A) %*% r)
        if ((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f)) < eps) break
    }
    return(u)
}


# GDM история -------------------------------------------------------------

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
GDM.history <- function(A, f, u, eps = 10e-4) {
    stopifnot(is.numeric(A), is.matrix(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    i <- 0
    u.hist <- matrix(u, nrow = dimA)
    t1 <- Sys.time()
    repeat {
        r <- A %*% u - f
        u <- u - ((t(t(A) %*% r) %*% (t(A) %*% r))/(t(A %*% t(A) %*% r) %*% (A %*% t(A) %*% r)))[1,1] * (t(A) %*% r)
        u.hist <- cbind(u.hist, u)
        i <- i + 10
        if ((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f)) < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}