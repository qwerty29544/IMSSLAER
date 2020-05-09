# CHebishev iterative method ----------------------------------------------

#' Chebishev iterative method [CIM] 
#' (Чебышевский итерационный метод)
#' @description The Chebyshev iterative method is a stationary 
#' iterative method, a variation of the simple iteration method, 
#' in which, before starting iterations, a vector of parameters 
#' is determined based on a given number of iterative layers. 
#' The larger the number of layers in the algorithm, the better 
#' will be the convergence of the method.
#' (Чебышевский итерационный метод представляет собой стационарный 
#' итерационный метод, разновидность метода простой итерации, в 
#' котором перед началом итераций определяется вектор параметров 
#' исходя из заданного числа итерационных слоёв. 
#' Чем больше число слоёв алгоритма, тем лучше будет сходимость метода.)
#' @param A - the original matrix of the operator equation - numeric 
#' only matrix (исходная матрица операторного уравнения - 
#' вещественная)
#' @param f - bias - numeric or complex vector (вектор свободных членов 
#' вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric 
#' or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param lambs - vector of numeric Operator's spectre: sigma(A) 
#' is Numeric [Real], it's can be a vector of all or several lambs, 
#' but there must be a min and max of eigen(A) (Вещественный вектор 
#' спектра оператора. Ограничением метода является вещественность 
#' спектра оператора, может быть вектором из всех, двух или одной 
#' точки спектра, однако данный вектор должен содержать его минимум 
#' и максимум для оценки спектрального диаметра) 
#' @param layers - the number of layers of the iterative method 
#' determines the degree of accuracy of calculations (число слоёв 
#' итерационного метода, определяет степень точности вычислений)
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
#' @examples A <- diag(seq(0.1, 1, 0.1), ncol = 10)
#' f <- rnorm(10) + 1i * rnorm(10)
#' u <- rnorm(10)
#' print(Chebishev(A, f, u, lambs = diag(A), layers = 5, eps = 10e-08))
#' all.equal(solve(A) %*% f, Chebishev(A, f, u, lambs = diag(A), layers = 5, eps = 10e-08))
Chebishev <- function(A, f, u, lambs, layers = 5, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A),
              is.numeric(f) || is.complex(f),
              is.numeric(u) || is.complex(u),
              is.numeric(eps), length(eps) == 1, 
              is.vector(lambs), is.numeric(lambs),
              is.atomic(eps),
              nrow(A) == ncol(A), ncol(A) == length(f), 
              length(f) == length(u),
              ncol(A) >= 2,
              is.numeric(iterations), length(iterations) == 1, 
              is.atomic(iterations),
              is.atomic(layers), layers >= 1, layers <= 1000)
    
    mu <- (max(lambs) + min(lambs))/2 + 
        ((max(lambs) - min(lambs))/2) * 
        cos(((2*(1:layers) - 1)*pi)/(2*layers))
    for (iter in 1:layers) {
        i <- 0
        repeat {
            u <- u - (1/mu[iter]) * (A %*% u - f)
            i <- i + 1
            if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                    / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
            if (i > iterations) {
                message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
                break
            }
        }
    }
    return(u)
}

# Chebishev history -------------------------------------------------------

#' Chebishev iterative method [CIM] 
#' (Чебышевский итерационный метод)
#' @description The Chebyshev iterative method is a stationary 
#' iterative method, a variation of the simple iteration method, 
#' in which, before starting iterations, a vector of parameters 
#' is determined based on a given number of iterative layers. 
#' The larger the number of layers in the algorithm, the better 
#' will be the convergence of the method.
#' (Чебышевский итерационный метод представляет собой стационарный 
#' итерационный метод, разновидность метода простой итерации, в 
#' котором перед началом итераций определяется вектор параметров 
#' исходя из заданного числа итерационных слоёв. 
#' Чем больше число слоёв алгоритма, тем лучше будет сходимость метода.)
#' @details This method is necessary to preserve the history of sequential 
#' calculation of an unknown vector in order to visualize the 
#' convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного 
#' вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - numeric 
#' only matrix (исходная матрица операторного уравнения - 
#' вещественная)
#' @param f - bias - numeric or complex vector (вектор свободных членов 
#' вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric 
#' or complex vector (начальное приближение неизвестного вектора 
#' - вещественный или комплексный вектор)
#' @param lambs - vector of numeric Operator's spectre: sigma(A) 
#' is Numeric [Real], it's can be a vector of all or several lambs, 
#' but there must be a min and max of eigen(A) (Вещественный вектор 
#' спектра оператора. Ограничением метода является вещественность 
#' спектра оператора, может быть вектором из всех, двух или одной 
#' точки спектра, однако данный вектор должен содержать его минимум 
#' и максимум для оценки спектрального диаметра) 
#' @param layers - the number of layers of the iterative method 
#' determines the degree of accuracy of calculations (число слоёв 
#' итерационного метода, определяет степень точности вычислений)
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
#' @examples A <- diag(seq(0.1, 1, 0.1), ncol = 10)
#' f <- rnorm(10) + 1i * rnorm(10)
#' u <- rnorm(10)
#' print(Chebishev.history(A, f, u, lambs = diag(A), layers = 5, eps = 10e-08))
Chebishev.history <- function(A, f, u, lambs, layers = 5, eps = 10e-4, iterations = 10000) {
    stopifnot(is.matrix(A),
              is.numeric(A),
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
    t1 <- Sys.time()
    mu <- (max(lambs) + min(lambs))/2 + 
        ((max(lambs) - min(lambs))/2) * 
        cos(((2*(1:layers) - 1)*pi)/(2*layers))
    iterate <- 0
    iterate2 <- 0
    u.hist <- matrix(u, ncol = 1, nrow = length(u))
    for (iter in 1:layers) {
        i <- 0
        repeat {
            # одно умножение матрицы на вектор
            iterate <- iterate + 1
            iterate2 <- c(iterate2, iterate)
            u <- u - (1/mu[iter]) * (A %*% u - f)
            u.hist <- cbind(u.hist, u)
            i <- i + 1
            if (abs((sqrt(t(A %*% u - f) %*% Conj(A %*% u - f))) 
                    / (sqrt(t(f) %*% Conj(f))))[1, 1] < eps) break
            if (i > iterations) {
                message("Iterations of the method may not come close to 
                    the final result / allowed number of iterations 
                    is exceeded")
                break
            }
        }
    }
    t2 <- Sys.time()
    return(list(num.iter = iterate2, var = u, var.hist = u.hist, 
                systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}
