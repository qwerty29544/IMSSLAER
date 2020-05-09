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

.circle <- function(x, y, radius, nv = 100, border = NULL, col = NA, 
                    lty = 1, density = NULL, angle = 45, lwd = 0.5) {
    ymult <- .getYmult()
    angle.inc <- 2 * pi/nv
    angles <- base::seq(0, 2 * pi - angle.inc, by = angle.inc)
    if (base::length(col) < base::length(radius))
        col <- base::rep(col, length.out = base::length(radius))
    for (circle in 1:base::length(radius)) {
        xv <- cos(angles) * radius[circle] + x
        yv <- sin(angles) * radius[circle] * ymult + y
        graphics::polygon(xv, yv, border = border, col = col[circle], 
                          lty = lty,
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
    # Функция из библиотеки 'plotrix' для рисования 
    # окружности с центром и радиусом
    .circle(x = Re(mu), y = Im(mu), radius = Re(R), 
            border = "red", lwd = 1.2)
    grid()
    grDevices::dev.off()
}

# Нахождение центра окружности спектра - отрезка на комплексной пл --------
# Нахождение такого мю по двум крайним точкам отрезка спектра линейного
# оператора, что он лежит на серединном перпендикуляре к отрезку z1z2
# и лежит на луче, исходящем из начала координат

.mu2 <- function(z1, z2) {
    stopifnot(is.numeric(z1) || is.complex(z1))
    stopifnot(is.numeric(z2) || is.complex(z2))
    return(
        ((z1 + z2) / 2) + 
               ((1i * Im(z1 * Conj(z2)) * (z2 - z1)) 
                / (2 * (Mod(z1 * Conj(z2)) + Re(z1 * Conj(z2)))))
    )
}


# Нахождение радиуса такого круга -----------------------------------------


.R2 <- function(z1, z2) {
    stopifnot(is.numeric(z1) || is.complex(z1))
    stopifnot(is.numeric(z2) || is.complex(z2))
    return( 
        sqrt((Mod(z1 - z2) * Mod(z1 - z2) 
              * Mod(Conj(z1) * z2)) 
             / (2 * (Mod(Conj(z1) * z2) + Re(Conj(z1) * z2))))
        )
}


# Нахождение центра окружности по 3 точкам --------------------------------


.mu3 <- function(z1, z2, z3) {
    stopifnot(is.numeric(z1) || is.complex(z1))
    stopifnot(is.numeric(z2) || is.complex(z2))
    stopifnot(is.numeric(z3) || is.complex(z3))
    return(
        1i * (Mod(z1)^2 * (z2 - z3) 
              + Mod(z2)^2 * (z3 - z1) 
              + Mod(z3)^2 * (z1 - z2)) / 
            (2 * Im(z1 * Conj(z2) + z2 * Conj(z3) + z3 * Conj(z1)))
        )
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
            Flags[j, i] <- ((abs(R[j]) > abs(mu[j] - lambs[i])) || 
                                (abs(R[j]) - abs(mu[j] - lambs[i]) == 0))
        }
        Flags[j, n + 1] <- abs(R[j]) < abs(mu[j])
    }
    return(Flags)
}

# Нахождение мю по входным значениям краёв спектра оператора --------------

#' Iterative parameter for GMSI iterations [muFind]
#' (Нахождение итерационного параметра для обобщенного метода простой итерации)
#' @description An algorithm for finding an iterative parameter for 
#' the generalized method of simple iteration based on information 
#' about the values of the spectrum point of a linear operator or 
#' on the basis of a priori knowledge about the location of the spectrum 
#' of operator A. For geometric reasons, a circle of smallest radius is 
#' constructed that does not contain the origin, which will contain all 
#' points of the spectrum of the operator.
#' (Алгоритм нахождения итерационного параметра для обобщенного метода 
#' простой итерации на основе информации о значениях точке спектра линейного 
#' оператора или на основе априорных знаний о расположении спектра 
#' оператора A. Из геометрических соображений строится круг наименьшего 
#' радиуса, не содержащий в себе начало координат, который будет содержать 
#' в себе все точки спектра оператора.)
#' @param lambs - input data for spectrum points 
#' (входные данные для точек спектра)
#' @param draw - boolean type variable that controls the drawing of 
#' the complex plane of the spectrum of the operator (пременная 
#' логического типа, управляет выводом рисунка комплексной плоскости 
#' спектра оператора)
#' @param path - operator plane image file name
#' (имя файла изображения комплексной плоскости оператора)
#'
#' @return mu - the complex value of the operator center 
#' (комплексное значение центра оператора);
#' complex.plot - ".jpg" image file of the operator's spectrum point on 
#' the complex plane and circle (файл изображения точке спектра оператора 
#' на комплексной плоскости и окружности)
#' @export
#'
#' @examples A <- diag(seq(0.1, 99.1, 1))
#' print(muFind(lambs = diag(A), draw = T))
#' 
#' A <- diag(seq(-0.1, -99.1, -1))
#' print(muFind(lambs = diag(A), draw = F))
#' 
#' print(muFind(lambs = c(5 + 5i, 5 - 5i, 4)))
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
    } else if (n == 2) {
        # оказывается верным является выражение 5 < 5/0
        mu <- .mu2(lambs[1], lambs[2])
        R <- .R2(lambs[1], lambs[2])
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
    } else if (n >= 3) {
        
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

#' Iterative parameter for GMSI iterations from file [muFind.File]
#' (Нахождение итерационного параметра для обобщенного метода простой итерации)
#' @description An algorithm for finding an iterative parameter for 
#' the generalized method of simple iteration based on information 
#' about the values of the spectrum point of a linear operator or 
#' on the basis of a priori knowledge about the location of the spectrum 
#' of operator A. For geometric reasons, a circle of smallest radius is 
#' constructed that does not contain the origin, which will contain all 
#' points of the spectrum of the operator.
#' (Алгоритм нахождения итерационного параметра для обобщенного метода 
#' простой итерации на основе информации о значениях точке спектра линейного 
#' оператора или на основе априорных знаний о расположении спектра 
#' оператора A. Из геометрических соображений строится круг наименьшего 
#' радиуса, не содержащий в себе начало координат, который будет содержать 
#' в себе все точки спектра оператора.)
#' @param input.file - file with input data for spectrum points 
#' (файл с входными данными для точек спектра)
#' @param output.file - file with input data for the center of the circle 
#' (файл с входными данными для центра окружности)
#' @param out - boolean variable, controls file output 
#' (переменная логического типа, управляет выводом в файл)
#' @param draw - boolean type variable that controls the drawing of 
#' the complex plane of the spectrum of the operator (пременная 
#' логического типа, управляет выводом рисунка комплексной плоскости 
#' спектра оператора)
#' @param plot.file - operator plane image file name
#' (имя файла изображения комплексной плоскости оператора)
#'
#' @return mu - the complex value of the operator center 
#' (комплексное значение центра оператора);
#' ouput.file - file with data of the center of the circle on the 
#' complex plane (файл с данными центра окружности на комплексной плоскости);
#' complex.plot - ".jpg" image file of the operator's spectrum point on 
#' the complex plane and circle (файл изображения точке спектра оператора 
#' на комплексной плоскости и окружности)
#' @export
#'
#' @examples
muFind.File <- function(input.file = "complexNumbeRS", output.file = "results", out = T,
                        draw = TRUE, plot.file = "complexPlot") {
    stopifnot(is.logical(draw), is.character(input.file), 
              is.character(output.file), is.character(plot.file),
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
    
    } else if (n == 2) {
        
        # оказывается верным является выражение 5 < 5/0
        mu <- .mu2(lambs[1], lambs[2])
        R <- .R2(lambs[1], lambs[2])
        if (abs(mu) <= abs(R)) stop("Окружность оператора лежит на начале координат")
        muflash <- c(muflash, mu)
        Rflash <- c(Rflash, R)
        
    } else if (n >= 3) {

        for (i in (1:(n - 1))) {
            for (j in ((i + 1):n)) {
                mu <- .mu2(lambs[i], lambs[j])
                R <- .R2(lambs[i], lambs[j])
                if (abs(mu) <= abs(R)) next
                if (any(round(abs(mu - lambs), digits = 8) > 
                        round(abs(R), digits = 8))) next
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
                    if (any(round(abs(mu - lambs), digits = 8) > 
                            round(abs(R), digits = 8))) next
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
# GMSI --------------------------------------------------------------------

#' Generalized method of simple iterations [GMSI]
#' (Обощённый метод простой итерации)
#' @description - Stationary iterative method for solving operator 
#' equations or systems of linear algebraic equations. The limitation 
#' is the absence among the values of the spectrum of the operator, 
#' the values of the mirror relative to the origin.
#' (Стационарный итерационный метод для решения операторных уравнений 
#' или систем линейных алгебраических уравнений. Ограничением ялвяется 
#' отсутствие среди значений спектра оператора, значений зеркальных 
#' относительно начала координат.)
#' @param A - the original matrix of the operator equation - 
#' numeric or complex matrix (исходная матрица операторного уравнения 
#' - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param mu - center of a circle on a complex plane containing points 
#' of the spectrum of the operator(центр окружности на комплексной 
#' плоскости, содержащей точки спектра оператора)
#' @param eps - accuracy of calculation of the desired vector - numeric 
#' (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return u - unknown vector in some approximation 
#' (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples A <- diag(rnorm(25, 50, 0.1) + 1i * rnorm(25, 50, 0.1))
#' f <- rnorm(25) + 1i * rnorm(25)
#' u <- rnorm(25)
#' print(IMSSLAER::GMSI(A, f, u, mu = IMSSLAER::muFind(lambs = diag(A), draw = F)))
GMSI <- function(A, f, u, mu, eps = 10e-4, iterations = 10000) {
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
              is.complex(mu) || is.numeric(mu))
    i <- 0
    # Итерации
    repeat {
        i <- i + 1
        u <- u - (1/mu) * (A %*% u - f)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    return(u)
}

# GMSI.mu -----------------------------------------------------------------

#' Generalized method of simple iterations [GMSI]
#' (Обощённый метод простой итерации)
#' @description - Stationary iterative method for solving operator 
#' equations or systems of linear algebraic equations. The limitation 
#' is the absence among the values of the spectrum of the operator, 
#' the values of the mirror relative to the origin.
#' (Стационарный итерационный метод для решения операторных уравнений 
#' или систем линейных алгебраических уравнений. Ограничением ялвяется 
#' отсутствие среди значений спектра оператора, значений зеркальных 
#' относительно начала координат.)
#' @param A - the original matrix of the operator equation - 
#' numeric or complex matrix (исходная матрица операторного уравнения 
#' - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param lambs - input data for spectrum points 
#' (входные данные для точек спектра)
#' @param eps - accuracy of calculation of the desired vector - numeric 
#' (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return u - unknown vector in some approximation 
#' (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples A <- diag(rnorm(25, 50, 0.1) + 1i * rnorm(25, 50, 0.1))
#' f <- rnorm(25) + 1i * rnorm(25)
#' u <- rnorm(25)
#' print(IMSSLAER::GMSI.mu(A, f, u, lambs = diag(A)))
GMSI.mu <- function(A, f, u, lambs, eps = 10e-4, iterations = 10000) {
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
              is.complex(lambs) || is.numeric(lambs))
    # Итерации
    i <- 0
    mu <- muFind(lambs = lambs, draw = F)[1]
    repeat {
        i <- i + 1
        u <- u - (1/mu) * (A %*% u - f)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    return(u)
}

# GMSI.mu.file -----------------------------------------------------------------

#' Generalized method of simple iterations [GMSI]
#' (Обощённый метод простой итерации)
#' @description - Stationary iterative method for solving operator 
#' equations or systems of linear algebraic equations. The limitation 
#' is the absence among the values of the spectrum of the operator, 
#' the values of the mirror relative to the origin.
#' (Стационарный итерационный метод для решения операторных уравнений 
#' или систем линейных алгебраических уравнений. Ограничением ялвяется 
#' отсутствие среди значений спектра оператора, значений зеркальных 
#' относительно начала координат.)
#' @param A - the original matrix of the operator equation - 
#' numeric or complex matrix (исходная матрица операторного уравнения 
#' - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param input.file - file with input data for spectrum points 
#' (файл с входными данными для точек спектра)
#' @param out - boolean variable, controls file output 
#' (переменная логического типа, управляет выводом в файл)
#' @param output.file - file with input data for the center of the circle 
#' (файл с входными данными для центра окружности)
#' @param draw - boolean type variable that controls the drawing of 
#' the complex plane of the spectrum of the operator (пременная 
#' логического типа, управляет выводом рисунка комплексной плоскости 
#' спектра оператора)
#' @param plot.file - operator plane image file name
#' (имя файла изображения комплексной плоскости оператора) 
#' @param eps - accuracy of calculation of the desired vector - numeric 
#' (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return u - unknown vector in some approximation 
#' (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples
GMSI.mu.file <- function(A, f, u, input.file = "complexNumbeRS", out = F,
                         output.file = "results", draw = F, 
                         plot.file = "complexPlot", eps = 10e-4,
                         iterations = 10000) {
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
              is.character(input.file), is.character(output.file), 
              is.character(plot.file), is.logical(draw), is.logical(out))
    # Итерации
    i <- 0
    mu <- muFind.File(input.file = input.file, output.file = output.file, draw = draw,
                      plot.file = plot.file, out = out)[1]
    repeat {
        i <- i + 1
        u <- u - (1/mu) * (A %*% u - f)
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    return(u)
}


# История по GMSI ---------------------------------------------------------

#' Generalized method of simple iterations [GMSI]
#' (Обощённый метод простой итерации)
#' @description - Stationary iterative method for solving operator 
#' equations or systems of linear algebraic equations. The limitation 
#' is the absence among the values of the spectrum of the operator, 
#' the values of the mirror relative to the origin.
#' (Стационарный итерационный метод для решения операторных уравнений 
#' или систем линейных алгебраических уравнений. Ограничением ялвяется 
#' отсутствие среди значений спектра оператора, значений зеркальных 
#' относительно начала координат.)
#' @details This method is necessary to preserve the history of 
#' sequential calculation of an unknown vector in order to visualize 
#' the convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного 
#' вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - 
#' numeric or complex matrix (исходная матрица операторного уравнения 
#' - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param mu - center of a circle on a complex plane containing points 
#' of the spectrum of the operator(центр окружности на комплексной плоскости, 
#' содержащей точки спектра оператора)
#' @param eps - accuracy of calculation of the desired vector - numeric 
#' (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история вычисления 
#' неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples A <- diag(rnorm(25, 50, 0.1) + 1i * rnorm(25, 50, 0.1))
#' f <- rnorm(25) + 1i * rnorm(25)
#' u <- rnorm(25)
#' print(IMSSLAER::GMSI.history(A, f, u, mu = IMSSLAER::muFind(lambs = diag(A), draw = F)))
GMSI.history <- function(A, f, u, mu, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps), length(eps) == 1, 
              is.atomic(eps),
              nrow(A) == ncol(A), ncol(A) == length(f), 
              length(f) == length(u),
              ncol(A) >= 2,
              is.numeric(iterations), 
              length(iterations) == 1, 
              is.atomic(iterations),
              is.complex(mu) || is.numeric(mu))
    # Итерации
    i <- 0 # начало итераций
    u.hist <- matrix(u, nrow = ncol(A))
    
    t1 <- Sys.time()
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        i <- i + 1 # обновляем счетчик итераций
        u.hist <- cbind(u.hist, u) # обновляем историю
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, 
                systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}


# GMSI.mu.history ---------------------------------------------------------

#' Generalized method of simple iterations [GMSI]
#' (Обощённый метод простой итерации)
#' @description - Stationary iterative method for solving operator 
#' equations or systems of linear algebraic equations. The limitation 
#' is the absence among the values of the spectrum of the operator, 
#' the values of the mirror relative to the origin.
#' (Стационарный итерационный метод для решения операторных уравнений 
#' или систем линейных алгебраических уравнений. Ограничением ялвяется 
#' отсутствие среди значений спектра оператора, значений зеркальных 
#' относительно начала координат.)
#' @details This method is necessary to preserve the history of 
#' sequential calculation of an unknown vector in order to visualize 
#' the convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного 
#' вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - 
#' numeric or complex matrix (исходная матрица операторного уравнения 
#' - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param lambs - input data for spectrum points 
#' (входные данные для точек спектра)
#' @param eps - accuracy of calculation of the desired vector - numeric 
#' (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история 
#' вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples A <- diag(rnorm(25, 50, 0.1) + 1i * rnorm(25, 50, 0.1))
#' f <- rnorm(25) + 1i * rnorm(25)
#' u <- rnorm(25)
#' print(IMSSLAER::GMSI.mu.history(A, f, u, lambs = diag(A)))
GMSI.mu.history <- function(A, f, u, lambs, eps = 10e-4, iterations = 10000) {
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
              is.complex(lambs) || is.numeric(lambs))
    # Итерации
    i <- 0 # начало итераций
    u.hist <- matrix(u, nrow = ncol(A))
    t1 <- Sys.time()
    mu <- muFind(lambs = lambs, draw = F)[1]
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        i <- i + 1 # обновляем счетчик итераций
        u.hist <- cbind(u.hist, u) # обновляем историю
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, 
                systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}

# GMSI.mu.file.history ----------------------------------------------------

#' Generalized method of simple iterations [GMSI]
#' (Обощённый метод простой итерации)
#' @description - Stationary iterative method for solving operator 
#' equations or systems of linear algebraic equations. The limitation 
#' is the absence among the values of the spectrum of the operator, 
#' the values of the mirror relative to the origin.
#' (Стационарный итерационный метод для решения операторных уравнений 
#' или систем линейных алгебраических уравнений. Ограничением ялвяется 
#' отсутствие среди значений спектра оператора, значений зеркальных 
#' относительно начала координат.)
#' @details This method is necessary to preserve the history of 
#' sequential calculation of an unknown vector in order to visualize 
#' the convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного 
#' вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - 
#' numeric or complex matrix (исходная матрица операторного уравнения 
#' - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор 
#' свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - 
#' numeric or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param input.file - file with input data for spectrum points 
#' (файл с входными данными для точек спектра)
#' @param out - boolean variable, controls file output 
#' (переменная логического типа, управляет выводом в файл)
#' @param output.file - file with input data for the center of the circle 
#' (файл с входными данными для центра окружности)
#' @param draw - boolean type variable that controls the drawing of 
#' the complex plane of the spectrum of the operator (пременная 
#' логического типа, управляет выводом рисунка комплексной плоскости 
#' спектра оператора)
#' @param plot.file - operator plane image file name
#' (имя файла изображения комплексной плоскости оператора) 
#' @param eps - accuracy of calculation of the desired vector - numeric 
#' (точность вычисления искомого вектора - вещественная)
#' @param iterations - the upper limit on the number of iterations when 
#' the method diverges (ограничение сверху на число итераций при 
#' расхождении метода)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история 
#' вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples
GMSI.mu.file.history <- function(A, f, u, input.file = "complexNumbeRS", 
                                 out = F, output.file = "results", 
                                 draw = F, plot.file = "complexPlot", 
                                 eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A) || is.complex(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps), length(eps) == 1, is.atomic(eps),
              nrow(A) == ncol(A), ncol(A) == length(f), 
              length(f) == length(u),ncol(A) >= 2,
              is.numeric(iterations), length(iterations) == 1, 
              is.atomic(iterations),
              is.character(input.file), is.character(output.file), 
              is.character(plot.file), is.logical(draw), is.logical(out))
    i <- 0
    u.hist <- matrix(u, nrow = ncol(A))
    t1 <- Sys.time()
    mu <- muFind.File(input.file = input.file, 
                      output.file = output.file, draw = draw,
                      plot.file = plot.file, out = out)[1]
    # Итерации
    repeat {
        u <- u - (1/mu) * (A %*% u - f)
        i <- i + 1 # обновляем счетчик итераций
        u.hist <- cbind(u.hist, u) # обновляем историю
        if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
        if (i > iterations) {
            message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
            break
        }
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, 
                systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}
