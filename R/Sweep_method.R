#' Title
#'
#' @param A - input triadiagonal matrix
#' @param f - input vector of free elements
#' @param U_g - variable of boundaries
#'
#' @return
#' @export
#'
#' @examples 
#' A <- cbind(c(4, -3, 0, 0), c(5, 4, 5, 0), c(0, 1, 6, 4), c(0, 0, 7, 3))
#' f <- c(3, 2, 4, 8)
#' Sweep_me(A, f)
Sweep_me <- function(A, f, U_g = 0) {
    stopifnot(is.numeric(A) || is.integer(A))
    stopifnot(is.matrix(A))
    stopifnot(is.vector(f))
    stopifnot(is.numeric(f) || is.integer(f))
    
    a <- numeric(nrow(A) - 1)
    b <- numeric(nrow(A) - 1)
    
    # Первые элементы a и b
    a[1] <- -A[1, 2] / A[1, 1]
    b[1] <- f[1] / A[1, 1]
    
    # Прямой проход метода
    for (i in 2:(nrow(A) - 1)) {
        a[i] <- -A[i, i + 1] / (A[i, i - 1] * a[i - 1] + A[i, i])
        b[i] <- (f[i] - A[i, i-1] * b[i - 1]) / (A[i, i - 1] * a[i - 1] + A[i, i])
    }
    
    u <- numeric(nrow(A))
    lenu <- length(u)
    
    # Последний элемент вектора u
    u[lenu] <- (f[lenu] - A[lenu, lenu - 1] * b[lenu - 1]) / 
        (A[lenu, lenu - 1] * a[lenu - 1] + A[lenu, lenu])
    # Обратный проход
    for (i in (lenu - 1):1) {
        u[i] <- a[i] * u[i + 1] + b[i]
    }
    rm("lenu")
    
    
    # Присоединение границ, если таковые заданы
    if ((length(U_g) == 2) & 
        is.vector(U_g) & 
        (is.numeric(U_g) || is.integer(U_g))) {
        u <- c(U_g[1], u, U_g[2])
    }
    
    return(u)
}



