# Метод является неустойчивым, необходимо сверху ограничивать число итераций
# Также со всеми методами проделать, обязательно

# Conjugate gradient method -----------------------------------------------
#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
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
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
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
