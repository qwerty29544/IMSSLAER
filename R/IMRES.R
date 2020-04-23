# IMRES -------------------------------------------------------------------

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
#' @param A 
#' @param f 
#' @param u 
#' @param eps 
#'
#' @return
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