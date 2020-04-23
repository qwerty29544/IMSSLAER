# Jacobi ------------------------------------------------------------------

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
Jacobi <- function(A, f, u, eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A)  || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) stop("Operator must have dim >= 2")
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be quadratic")
    # Проверка диагонального преобладания
    if (!identical(diag(A), apply(A, 1, max))) stop("Operator must have a diagonal priority!")
    repeat {
        TempX <- f
        for (i in 1:dimA) {
            for (g in 1:dimA) {
                if (i != g) {
                    TempX[i] <-  TempX[i] - A[i, g] * u[g]
                }
            }
        }
        TempX <- TempX / diag(A)
        norm <- max(abs(u - TempX))
        u <- TempX
        if (norm < eps) break
    }
    return(u)
}


# Jacobi история ----------------------------------------------------------

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
Jacobi.history <- function(A, f, u, eps = 10e-4) {
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f), is.numeric(u), is.numeric(eps))
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
        
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
        
    }
    # Проверка диагонального преобладания
    if (!all.equal(diag(A), apply(A, 1, max))) {
        stop("Operator must have a diagonal priority!")
        
    }
    iter <- 0
    u.hist <- matrix(u ,nrow = length(u))
    t1 <- Sys.time()
    repeat {
        TempX <- f
        for (i in 1:dimA) {
            for (g in 1:dimA) {
                if (i != g) {
                    TempX[i] <-  TempX[i] - A[i, g] * u[g]
                }
            }
        }
        TempX <- TempX / diag(A)
        norm <- max(abs(u - TempX))
        u <- TempX
        iter <- iter + 1
        u.hist <- cbind(u.hist, u)
        if (norm < eps) break
        
    }
    t2 <- Sys.time()
    return(list(num.iter = iter, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}
