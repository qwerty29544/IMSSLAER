# ПРоверка на пуш из линукса
# Вспомогательная функция для ортисовки графиков --------------------------

.getYmult <- function() {
    if (grDevices::dev.cur() == 1) {
        base::warning("No graphics device open.")
        ymult <- 1
    }
    else {
        xyasp <- graphics::par("pin")
        xycr <- base::diff(graphics::par("usr"))[c(1, 3)]
        ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
    }
    return(ymult)
}


# Функция отрисовки окружности --------------------------------------------

#' Title
#'
#' @param x 
#' @param y 
#' @param radius 
#' @param nv 
#' @param border 
#' @param col 
#' @param lty 
#' @param density 
#' @param angle 
#' @param lwd 
#'
#' @return
#' @export
#'
#' @examples
.circle <- function(x, y, radius, nv = 100, border = NULL,
                    col = NA, lty = 1, density = NULL, angle = 45, lwd = 0.5) {
    ymult <- .getYmult()
    angle.inc <- 2 * pi/nv
    angles <- base::seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (base::length(col) < base::length(radius))
        col <- base::rep(col, length.out = base::length(radius))
    for (circle in 1:base::length(radius)) {
        xv <- cos(angles) * radius[circle] + x
        yv <- sin(angles) * radius[circle] * ymult + y
        graphics::polygon(xv, yv, border = border, col = col[circle], lty = lty,
                density = density, angle = angle, lwd = lwd)
    }
    invisible(list(x = xv, y = yv))
}


# Отрисовка комплексной плоскости -----------------------------------------

.drawComplex <- function(lambs, mu, R, fileName = "complexPlot") {
    stopifnot(is.complex(lambs) || is.numeric(lambs))
    stopifnot(is.complex(mu) || is.numeric(mu))
    stopifnot(is.complex(R) || is.numeric(R))
    dimA <- length(lambs)
    grDevices::jpeg(paste0(fileName,".jpg"), width = 900, height = 900)
    plot(x = c(Re(lambs), Re(mu)), asp = 1,
         y = c(Im(lambs), Im(mu)),
         xlim = c(Re(mu) - Re(R), Re(mu) + Re(R)),
         ylim = c(Im(mu) - Re(R), Im(mu) + Re(R)),
         main = "Спектр на комплексной плоскости линейного оператора",
         pch = 19,
         col = "black",
         cex = 1.5,
         cex.lab = 1.5,
         cex.main = 1.5)
    lines(x = Re(lambs),
          y = Im(lambs),
          type = "l")
    lines(x = c(Re(lambs[1]), Re(lambs[dimA])),
          y = c(Im(lambs[1]), Im(lambs[dimA])))
    abline(h = 0, lty = 1, lwd = 1.2)
    abline(v = 0, lty = 1, lwd = 1.2)
    # Функция из библиотеки 'plotrix' для рисования окружности с центром и радиусом
    .circle(x = Re(mu), y = Im(mu), radius = Re(R), border = "red", lwd = 1.2)
    grid()
    grDevices::dev.off()
}

# Нахождение центра окружности спектра - отрезка на комплексной пл --------
# Нахождение такого мю по двум крайним точкам отрезка спектра линейного
# оператора, что он лежит на серединном перпендикуляре к отрезку z1z2
# и лежит на луче, исходящем из начала координат

#' Title
#'
#' @param z1 - комплексное число, обозначающее начало отрезка спектра оператора на комплексной плоскости
#' @param z2 - комплексное число, обозначающее конец отрезка спектра оператора на комплексной плоскости
#'
#' @return mu - центр окружности на комплексной плоскости, выпуклая оболочка спектра
#' @export
#'
#' @examples
.mu2 <- function(z1, z2) {
    stopifnot(is.numeric(z1) || is.complex(z1))
    stopifnot(is.numeric(z2) || is.complex(z2))
    return(((z1 + z2) / 2) + ((1i * Im(z1 * Conj(z2)) * (z2 - z1)) / (2 * (Mod(z1 * Conj(z2)) + Re(z1 * Conj(z2))))))
}


# Нахождение радиуса такого круга -----------------------------------------

.R2 <- function(z1, z2) {
    stopifnot(is.numeric(z1) || is.complex(z1))
    stopifnot(is.numeric(z2) || is.complex(z2))
    return( sqrt((Mod(z1 - z2) * Mod(z1 - z2) * Mod(Conj(z1) * z2)) / (2 * (Mod(Conj(z1) * z2) + Re(Conj(z1) * z2)))))
}


# Нахождение центра окружности по 3 точкам --------------------------------

.mu3 <- function(z1, z2, z3) {
    stopifnot(is.numeric(z1) || is.complex(z1))
    stopifnot(is.numeric(z2) || is.complex(z2))
    stopifnot(is.numeric(z3) || is.complex(z3))
    return(1i * (Mod(z1)^2 * (z2 - z3) + Mod(z2)^2 * (z3 - z1) + Mod(z3)^2 * (z1 - z2)) / (2 * Im(z1 * Conj(z2) + z2 * Conj(z3) + z3 * Conj(z1))))
}


# Нахождение радиуса окружности по 3 точкам и центру ----------------------

.R3 <- function(mu, z) {
    stopifnot(is.numeric(mu) || is.complex(mu))
    stopifnot(is.numeric(z) || is.complex(z))
    return(sqrt(Mod(z - mu)^2))
}


# Функция перебора точек спектра ------------------------------------------

.FlagsCalc <- function(mu, R, n, lambs) {
    Flags <- matrix(FALSE, nrow = length(R), ncol = n + 1)
    for (j in 1:length(R)) {
        for (i in 1:n) {
            Flags[j, i] <- ((abs(R[j]) > abs(mu[j] - lambs[i])) || (abs(R[j]) - abs(mu[j] - lambs[i]) == 0))
        }
        Flags[j, n + 1] <- abs(R[j]) < abs(mu[j])
    }
    return(Flags)
}

# Нахождение мю по входным значениям краёв спектра оператора --------------

muFind <- function(lambs, draw = TRUE, path = "complexPlot") {
    # Проверка типа переменной - вектора
    stopifnot(is.numeric(lambs) || is.complex(lambs))
    stopifnot(is.logical(draw))
    n <- length(lambs) # Переменная количества считанных данных
    if (n <= 1) stop("dim of linear operator is emprty or not valid")
    # Если спектр оператора всё-таки действительная ось
    Rflash <- numeric(0)
    muflash <- complex(0)
    if (is.numeric(lambs) || all(Im(lambs) == 0)) {
        mu <- (max(Re(lambs)) + min(Re(lambs))) / 2 # считаем центр по двум крайним точкам
        R <- (max(Re(lambs)) - min(Re(lambs))) / 2 # считаем радиум круга по двум крайним точкам
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
    }
    # оказывается верным является выражение 5 < 5/0
    if (n == 2) {
        mu <- .mu2(lambs[1], lambs[2])
        R <- .R2(lambs[1], lambs[2])
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
    }
    if (n >= 3) {

        for (i in (1:(n - 1))) {
            for (j in ((i + 1):n)) {
                mu <- .mu2(lambs[i], lambs[j])
                R <- .R2(lambs[i], lambs[j])
                if (abs(mu) <= abs(R)) next
                if (any(round(abs(mu - lambs), digits = 8) > round(abs(R), digits = 8))) next
                muflash <- c(muflash, mu)
                Rflash <- c(Rflash, R)
            }
        }

        for (i in (1:(n - 2))) {
            for (j in ((i + 1):(n - 1))) {
                for (k in ((j + 1):n)) {
                    mu <- .mu3(lambs[i], lambs[j], lambs[k])
                    R <- .R3(mu, lambs[i])
                    if (abs(mu) <= abs(R)) next
                    if (any(round(abs(mu - lambs), digits = 8) > round(abs(R), digits = 8))) next
                    muflash <- c(muflash, mu)
                    Rflash <- c(Rflash, R)
                }
            }
        }
    }
    mu <- muflash[which.min(Rflash)]
    R <- min(Rflash)
    if (draw) .drawComplex(lambs = lambs, mu = mu, R = R, fileName = path)
    return(mu)
}


# muFind.File -------------------------------------------------------------

muFind.File <- function(input.file = "complexNumbeRS", output.file = "results", out = T,
                        draw = TRUE, plot.file = "complexPlot") {
    stopifnot(is.logical(draw), is.character(input.file), is.character(output.file), is.character(plot.file),
              is.logical(out))
    # Считывание данных о точках многоугольника из текстового документа
    lambs <- as.complex(read.table(file = input.file, header = F)[[1]])
    stopifnot((is.numeric(lambs) == TRUE) || (is.complex(lambs) == TRUE))
    cat(lambs, "\n") # Вывод считанных данных в консоль как есть
    n <- length(lambs) # Переменная количества считанных данных
    if (n <= 1) stop("dim of linear operator is emprty or not valid")
    # Если спектр оператора всё-таки действительная ось
    Rflash <- numeric(0)
    muflash <- complex(0)
    if (is.numeric(lambs) || all(Im(lambs) == 0)) {
        mu <- (max(Re(lambs)) + min(Re(lambs))) / 2 # считаем центр по двум крайним точкам
        R <- (max(Re(lambs)) - min(Re(lambs))) / 2 # считаем радиум круга по двум крайним точкам
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
    }
    # оказывается верным является выражение 5 < 5/0
    if (n == 2) {
        mu <- .mu2(lambs[1], lambs[2])
        R <- .R2(lambs[1], lambs[2])
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
    }
    if (n >= 3) {

        for (i in (1:(n - 1))) {
            for (j in ((i + 1):n)) {
                mu <- .mu2(lambs[i], lambs[j])
                R <- .R2(lambs[i], lambs[j])
                if (abs(mu) <= abs(R)) next
                if (any(round(abs(mu - lambs), digits = 8) > round(abs(R), digits = 8))) next
                muflash <- c(muflash, mu)
                Rflash <- c(Rflash, R)
            }
        }

        for (i in (1:(n - 2))) {
            for (j in ((i + 1):(n - 1))) {
                for (k in ((j + 1):n)) {
                    mu <- .mu3(lambs[i], lambs[j], lambs[k])
                    R <- .R3(mu, lambs[i])
                    if (abs(mu) <= abs(R)) next
                    if (any(round(abs(mu - lambs), digits = 8) > round(abs(R), digits = 8))) next
                    muflash <- c(muflash, mu)
                    Rflash <- c(Rflash, R)
                }
            }
        }
    }
    mu <- muflash[which.min(Rflash)]
    if (length(mu) == 0) stop("Нет подходящего центра")
    R <- min(Rflash)
    results <- c(Center = mu, Radius = R, Ro = sqrt(abs(R)^2/abs(mu)^2))
    if (draw) .drawComplex(lambs = lambs, mu = mu, R = R, fileName = plot.file)
    if (out) write.table(results, file = output.file, quote = FALSE)
    return(mu)
}


# muFindIter --------------------------------------------------------------


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

# GMSI --------------------------------------------------------------------

GMSI <- function(A, f, u, mu, eps = 10e-4) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps),
              is.complex(mu) || is.numeric(mu))
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    # Итерации
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        if ((sqrt(abs(t(A %*% u - f) %*% (A %*% u - f)))) / (sqrt(abs(t(f) %*% f))) < eps) break
    }
    return(u)
}


# GMSI.mu -----------------------------------------------------------------

GMSI.mu <- function(A, f, u, lambs, eps = 10e-4) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps),
              is.complex(lambs) || is.numeric(lambs))
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    # Итерации
    mu <- muFind(lambs = lambs, draw = F)[1]
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        if ((sqrt(abs(t(A %*% u - f) %*% (A %*% u - f)))) / (sqrt(abs(t(f) %*% f))) < eps) break
    }
    return(u)
}

# GMSI.mu -----------------------------------------------------------------

GMSI.mu.file <- function(A, f, u, input.file = "complexNumbeRS", out = F,
                         output.file = "results", draw = F, plot.file = "complexPlot", eps = 10e-4) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps))
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    # Итерации
    mu <- muFind.File(input.file = input.file, output.file = output.file, draw = draw,
                      plot.file = plot.file, out = out)[1]
    repeat{
        u <- u - (1/mu) * (A %*% u - f)
        if ((sqrt(abs(t(A %*% u - f) %*% (A %*% u - f)))) / (sqrt(abs(t(f) %*% f))) < eps) break
    }
    return(u)
}


# История по GMSI ---------------------------------------------------------

GMSI.history <- function(A, f, u, mu, eps = 10e-4) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps),
              is.complex(mu) || is.numeric(mu))
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    # Итерации
    i <- 0 # начало итераций
    u.hist <- matrix(u, nrow = dimA)

    t1 <- Sys.time()
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        i <- i + 1 # обновляем счетчик итераций
        u.hist <- cbind(u.hist, u) # обновляем историю
        if ((sqrt(abs(t(A %*% u - f) %*% (A %*% u - f)))) / (sqrt(abs(t(f) %*% f))) < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}

GMSI.mu.history <- function(A, f, u, lambs, eps = 10e-4) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps),
              is.complex(lambs) || is.numeric(lambs))
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    # Итерации
    i <- 0 # начало итераций
    u.hist <- matrix(u, nrow = dimA)
    t1 <- Sys.time()
    mu <- muFind(lambs = lambs, draw = F)[1]
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        i <- i + 1 # обновляем счетчик итераций
        u.hist <- cbind(u.hist, u) # обновляем историю
        if ((sqrt(abs(t(A %*% u - f) %*% (A %*% u - f)))) / (sqrt(abs(t(f) %*% f))) < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}

GMSI.mu.file.history <- function(A, f, u, input.file = "complexNumbeRS", out = F,
                                 output.file = "results", draw = F, plot.file = "complexPlot", eps = 10e-4) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps))
    # Получение размерности матрицы А
    dimA <- dim(A)[1]
    # Проверка на n >= 2
    if (dimA[1] < 2) {
        stop("Operator must have dim >= 2")
    }
    # Проверка на размерность матрицы оператора
    if (dim(A)[1] != dim(A)[2]) {
        stop("Operator must be quadratic")
    }
    # Итерации
    i <- 0 # начало итераций
    u.hist <- matrix(u, nrow = dimA)
    t1 <- Sys.time()
    mu <- muFind.File(input.file = input.file, output.file = output.file, draw = draw,
                      plot.file = plot.file, out = out)[1]
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        i <- i + 1 # обновляем счетчик итераций
        u.hist <- cbind(u.hist, u) # обновляем историю
        if ((sqrt(abs(t(A %*% u - f) %*% (A %*% u - f)))) / (sqrt(abs(t(f) %*% f))) < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}


# SIM ---------------------------------------------------------------------

SIM <- function(A, f, u, eps = 10e-4) {
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
    B <- diag(1, nrow = dimA, ncol = dimA) - A
    repeat {
        u <- B %*% u + f
        if ((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f)) < eps) break
    }
    return(u)
}


# SIM история -------------------------------------------------------------

SIM.history <- function(A, f, u, eps = 10e-4) {
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
    i <- 0
    u.hist <- matrix(u, nrow = dimA)
    t1 <- Sys.time()
    B <- diag(1, nrow = dimA, ncol = dimA) - A
    repeat {
        u <- B %*% u + f
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        if ((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f)) < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}


# Jacobi ------------------------------------------------------------------

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


# IMRES -------------------------------------------------------------------

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



# Итерации Чебышев ---------------------------------------------------------------

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


