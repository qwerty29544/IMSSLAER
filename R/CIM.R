# Итерации Чебышев ---------------------------------------------------------------

#' Title
#'
#' @param A 
#' @param f 
#' @param u 
#' @param lambs 
#' @param layers 
#' @param eps 
#'
#' @return
#' @export
#'
#' @examples
Chebishev <- function(A, f, u, lambs, layers, eps = 10e-4) {
    stopifnot(is.numeric(A), is.matrix(A), is.numeric(f), is.numeric(u),
              is.numeric(lambs), is.integer(layers), is.numeric(eps))
    if (dim(A)[1] < 2) stop("Linear operator has dim less than 2x2")
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be a quadratic")
    mu <- (max(lambs) + min(lambs))/2 + ((max(lambs) - min(lambs))/2) * cos(((2*(1:layers) - 1)*pi)/(2*layers))
    repeat {
        for (iter in (1:layers)) {
            u <- u - (1/mu[iter]) * (A %*% u - f)
        }
        if (((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f)))[1,1] < eps) break
    }
    return(u)
}

# История Чебышев ---------------------------------------------------------------

#' Title
#' 
#' @details This method is necessary to preserve the history of sequential calculation of an unknown vector in order to visualize the convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A 
#' @param f 
#' @param u 
#' @param lambs 
#' @param layers 
#' @param eps 
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' num.layers - number of layers used (число использованных слоёв); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples
Chebishev.history <- function(A, f, u, lambs, layers, eps = 10e-4) {
    stopifnot(is.numeric(A), is.matrix(A), is.numeric(f), is.numeric(u),
              is.numeric(lambs), is.integer(layers), is.numeric(eps))
    if (dim(A)[1] < 2) stop("Linear operator has dim less than 2x2")
    if (dim(A)[1] != dim(A)[2]) stop("Operator must be a quadratic")
    mu <- (max(lambs) + min(lambs))/2 + ((max(lambs) - min(lambs))/2) * cos(((2*(1:layers) - 1)*pi)/(2*layers))
    i <- 0
    num.layers <- 0
    u.hist <- matrix(u, nrow = length(u))
    t1 <- Sys.time()
    repeat {
        num.layers <- num.layers + 1
        for (iter in (1:layers)) {
            u <- u - (1/mu[iter]) * (A %*% u - f)
            u.hist <- cbind(u.hist, u)
            i <- i + 1
        }
        if (((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f)))[1,1] < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, num.layers = num.layers, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}
