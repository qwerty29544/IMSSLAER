
# Biconjugate gradient method, BiCGM --------------------------------------

# https://ru.wikipedia.org/wiki/%D0%9C%D0%B5%D1%82%D0%BE%D0%B4_%D0%B1%D0%B8%D1%81%D0%BE%D0%BF%D1%80%D1%8F%D0%B6%D1%91%D0%BD%D0%BD%D1%8B%D1%85_%D0%B3%D1%80%D0%B0%D0%B4%D0%B8%D0%B5%D0%BD%D1%82%D0%BE%D0%B2

# Алгоритм для действительных матриц
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
BiCGM <- function(A, f, u, eps = 0.000001) {
    
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    
    dimA <- dim(A)[1]
    
    # Проверка на n >= 2
    if (dimA[1] < 2) stop("Operator must have dim >= 2")
    
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be quadratic")
    
    r <- f - A %*% u
    p <- r
    z <- r
    s <- r
    
    repeat {
        alpha <- ((t(p) %*% r) / (t(s) %*% (A %*% z)))[1, 1]
        u <- u + alpha * z
        r1 <- r - alpha * (A %*% z)
        p1 <- p - alpha * (t(A) %*% s)
        beta <- ((t(p1) %*% r1) / (t(p) %*% r))[1, 1]
        z <- r1 + beta * z
        s <- p1 + beta * s
        if ((sqrt(t(r1) %*% r1) / sqrt(t(f) %*% f) < eps) & (sqrt(t(p1) %*% p1) / sqrt(t(f) %*% f) < eps)) break
        r <- r1
        p <- p1
        rm(r1)
        rm(p1)
    }
    
    return(u)

}



# BiCGM.history -----------------------------------------------------------

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
BiCGM.history <- function(A, f, u, eps = 0.0000001) {
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
    p <- r
    z <- r
    s <- r
    
    repeat {
        alpha <- ((t(p) %*% r) / (t(s) %*% (A %*% z)))[1, 1]
        u <- u + alpha * z
        r1 <- r - alpha * (A %*% z)
        p1 <- p - alpha * (t(A) %*% s)
        beta <- ((t(p1) %*% r1) / (t(p) %*% r))[1, 1]
        
        i <- i + 3
        u.hist <- cbind(u.hist, u)
        
        z <- r1 + beta * z
        s <- p1 + beta * s
        if ((sqrt(t(r1) %*% r1) / sqrt(t(f) %*% f) < eps) & (sqrt(t(p1) %*% p1) / sqrt(t(f) %*% f) < eps)) break
        r <- r1
        p <- p1
        rm(r1)
        rm(p1)
    }
    
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))

}


# testing -----------------------------------------------------------------


A <- diag(rnorm(25), nrow = 5, ncol = 5)
f <- rnorm(5)
u <- rnorm(5)

BiCGM.history(A, f, u)
solve(A) %*% f
