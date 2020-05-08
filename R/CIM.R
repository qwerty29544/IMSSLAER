
# CHebishev iteration method ----------------------------------------------


Chebishev <- function(A, f, u, lambs, layers = 5, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps),
              nrow(A) == ncol(A), ncol(A) == length(f), 
              length(f) == length(u),
              ncol(A) >= 2,
              is.numeric(iterations), length(iterations) == 1, 
              is.atomic(iterations),
              is.atomic(layers), layers >= 1, layers <= 1000)
    
    mu <- (max(lambs) + min(lambs))/2 + ((max(lambs) - min(lambs))/2) * cos(((2*(1:layers) - 1)*pi)/(2*layers))
    repeat {
        for (iter in (1:layers)) {
            u <- u - (1/mu[iter]) * (A %*% u - f)
        }
        if (((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f)))[1,1] < eps) break
    }
    return(u)
}


# Chebishev history -------------------------------------------------------


Chebishev.history <- function(A, f, u, lambs, layers = 5, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps),
              nrow(A) == ncol(A), ncol(A) == length(f), 
              length(f) == length(u),
              ncol(A) >= 2,
              is.numeric(iterations), length(iterations) == 1, 
              is.atomic(iterations),
              is.atomic(layers), layers >= 1, layers <= 1000)
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
