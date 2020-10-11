# MariaRemark -------------------------------------------------------------


# Out_of_trend ------------------------------------------------------------

MariaRemark.out_of_trend <- function(ts, delta, 
                                     mode = c("AP", "GP_log", "GP_lin", "W", "HP")) {
    lents <- length(ts)
    if (mode == "GP_log") {
        return(
            c(rep(0, times = delta), 
              log((ts[1:(lents - 2*delta)] * ts[(2*delta + 1):lents]) / 
                      (ts[(delta + 1):(lents - delta)]^2)), 
              rep(0, times = delta))
        )
    } else if (mode == "GP_lin") {
        return(
            c(rep(0, times = delta), 
              (ts[1:(lents - 2*delta)] * ts[(2*delta + 1):lents]) - 
                  (ts[(delta + 1):(lents - delta)]^2), 
              rep(0, times = delta))
        )
    } else if (mode == "W") {
        yt_ <- ts[1:(lents - 3*delta)]
        yt <- ts[(1 + delta):(lents - 2*delta)]
        yt_1 <- ts[(1 + 2*delta):(lents - delta)]
        yt_2 <- ts[(1 + 3*delta):(lents)]
        return(
            c(rep(0, times = delta),
              log((yt_ * yt_2 + yt * yt_1)/(yt_ * yt_1 + yt * yt_2)),
              rep(0, times = 2 * delta))
        )
    } else if (mode == "HP") {
        return(
            c(rep(0, times = delta), 
              log((ts[1:(lents - 2*delta)]^2 + ts[(2*delta + 1):lents]^2) / 
                      ((ts[1:(lents - 2*delta)] + ts[(2*delta + 1):lents]) * 
                           ts[(delta + 1):(lents - delta)])), 
              rep(0, times = delta))
        )
    } else {
        return(
            c(rep(0, times = delta), 
              log((ts[1:(lents - 2*delta)] + ts[(2*delta + 1):lents]) / 
                      (2 * ts[(delta + 1):(lents - delta)])), 
              rep(0, times = delta))
        )
    }
}


# Alter_Johns -------------------------------------------------------------



MariaRemark.Alter_Johns <- function(ts, p = 1) {
    lents <- length(ts)
    a <- numeric(length = lents)
    a[1] <- 0
    for (i in 1:(lents-1)) {
        a[i + 1] <- 1/(lents - i) * sum(abs(ts[1:(lents - i)] - ts[(1 + i):lents])^(p))^(1/p)
    }
    return(a)
}


# Shift_function ----------------------------------------------------------

MariaRemark.Shift_function <- function(ts, delta_bounds = c(1, 100), 
                                       mode = "AP", 
                                       p = 1) {
    lents <- length(ts)
    Z <- matrix(ncol = lents, nrow = (delta_bounds[2] - delta_bounds[1] + 1))
    for (row in ((delta_bounds[1]:delta_bounds[2]) - delta_bounds[1] + 1)) {
        Z[row, ] <- c(
            MariaRemark.Alter_Johns(
                MariaRemark.out_of_trend(ts = ts, 
                                         delta = (row - 1 + delta_bounds[1]), 
                                         mode = mode), 
                p = p
            )
        )
    }
    return(Z)
}





# A <- MariaRemark.Shift_function(sin(seq(-10, 10, 0.01)) + 2)
# plotly::plot_ly(z = ~A, type = "surface")



# (sin(seq(-10, 10, 0.01)) + 3) %>% 
#     MariaRemark.out_of_trend(delta = 100, mode = "AP") %>% 
#     MariaRemark.Alter_Johns(p = 1) %>% 
#     plot(type = "l")



     
# t <- (1:100000)[is_prime(1:100000)]
# x <- t * cos(t)
# y <- t * sin(t)
# plotly::plot_ly(x = x, y = y, type = "scatter", size = I(0.7)) %>% plotly::layout(scene = list(aspectration=list(x=1,y=1)))
