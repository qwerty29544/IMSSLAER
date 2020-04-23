# GDM ---------------------------------------------------------------------

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
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
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