# Метод является неустойчивым, необходимо сверху ограничивать число итераций
# Также со всеми методами проделать, обязательно

# Conjugate gradient method -----------------------------------------------
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
CGM <- function(A, f, u, eps = 0.00001) {
    # https://en.wikipedia.org/wiki/Conjugate_gradient_method
    
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) stop("Operator must have dim >= 2")
    
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be quadratic")
    
    r <- f - A %*% u
    z <- r
    repeat {
        alpha = (t(r) %*% r) / (t(A %*% z) %*% z)
        u = u + alpha[1,1] * z
        r1 = r - alpha[1,1] * (A %*% z)
        beta <- (t(r1) %*% r1) / (t(r) %*% r)
        z = r1 + beta[1,1] * z
        if ((sqrt(t(r1) %*% r1) / sqrt(t(f) %*% f)) < eps) break
        r = r1
        rm(r1)
    }
    return(u)
}


# CGM.history -------------------------------------------------------------

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
CGM.history <- function(A, f, u, eps = 0.00001) {
    # https://en.wikipedia.org/wiki/Conjugate_gradient_method
    
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) stop("Operator must have dim >= 2")
    
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be quadratic")
    
    i <- 0
    u.hist <- matrix(u, nrow = dimA)
    t1 <- Sys.time()

    r <- f - A %*% u
    z <- r
    repeat{
        alpha = (t(r) %*% r) / (t(A %*% z) %*% z)
        u = u + alpha[1,1] * z
        
        i <- i + 2
        u.hist <- cbind(u.hist, u)
        
        r1 = r - alpha[1,1] * (A %*% z)
        beta <- (t(r1) %*% r1) / (t(r) %*% r)
        z = r1 + beta[1,1] * z
        if ((sqrt(t(r1) %*% r1) / sqrt(t(f) %*% f)) < eps) break
        r = r1
        rm(r1)
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}


# testing -----------------------------------------------------------------

A <- diag(rnorm(25), nrow = 5, ncol = 5)
u <- rnorm(5)
f <- rnorm(5)

CGM.history(A, f, u, eps = 0.000001)
solve(A) %*% f


(solve(A) %*% f) - CGM(A, u, f)
