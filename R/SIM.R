# simple iterations method ---------------------------------------------------------------------

#' Simple iteration method 
#' (Метод простой итерации)
#' @description A stationary iterative method for solving systems of linear algebraic equations. The method is based on the operation of reducing the operator equation to an iterative form, in which, when the matrix is multiplied by a vector, the unknown vector u approaches the real desired solution, in form: Au = f. A significant limitation of this method is the need for strict inequality for the spectrum of the operator in order to converge the method: sigma(A) < 1
#' (Стационарный итерационный метод решения систем линейных алгебраических уравнений. В основе метода лежит операция приведения операторного уравнения к итерационной форме, в которой при умножении матрицы на вектор происходит приближение неизвестного вектора u к реальному искомому решению, в форме: Au = f. Существенным ограничением данного метода является необходимость строгого неравенства для спектра оператора в целях сходимости метода: sigma(A) < 1)
#' @param A - the original matrix of the operator equation - numeric or complex matrix (исходная матрица операторного уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric or complex vector (начальное приближение неизвестного вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - numeric (точность вычисления искомого вектора - вещественная)
#'
#' @return u - unknown vector in some approximation (неизвестный вектор в некотором приближении)
#' @export
#'
#' @examples A <- diag(c(0.3, 0.4, 0.5), nrow = 3, ncol = 3)
#' f <- rnorm(3)
#' u <- rnorm(3)
#' result <- SIM(A = A, u = u, f = f, eps = 10e-4)
#' print(result)
#' 
#' A <- diag(c(0.5, 0.6 + 0.3i, 0.8, 0.2), nrow = 4, ncol = 4)
#' f <- rnorm(4) + 1i * rnorm(4)
#' u <- rnorm(4) + 1i * rnorm(4)
#' result <- SIM(A = A, u = u, f = f, eps = 10e-4)
#' print(result)
SIM <- function(A, f, u, eps = 10e-4) {
    # Проверка на N >= 2 - размерность
    # Все размерности совпадают
    # Все типы данных соответствуют ограничениям
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f) || is.complex(f), is.numeric(u) || is.complex(u), is.numeric(eps), nrow(A) == ncol(A), ncol(A) == length(f), length(f) == length(u), ncol(A) >= 2)
    dimA <- dim(A)[1]
    B <- diag(1, nrow = dimA, ncol = dimA) - A
    repeat {
        u <- B %*% u + f
        if (abs((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f))) < eps) break
    }
    return(u)
}

# simple iterations method history -------------------------------------------------------------

#' Simple iteration method history
#' (Метод простой итерации)
#' @description A stationary iterative method for solving systems of linear algebraic equations. The method is based on the operation of reducing the operator equation to an iterative form, in which, when the matrix is multiplied by a vector, the unknown vector u approaches the real desired solution, in form: Au = f. A significant limitation of this method is the need for strict inequality for the spectrum of the operator in order to converge the method: sigma(A) < 1
#' (Стационарный итерационный метод решения систем линейных алгебраических уравнений. В основе метода лежит операция приведения операторного уравнения к итерационной форме, в которой при умножении матрицы на вектор происходит приближение неизвестного вектора u к реальному искомому решению, в форме: Au = f. Существенным ограничением данного метода является необходимость строгого неравенства для спектра оператора в целях сходимости метода: sigma(A) < 1)
#' @details This method is necessary to preserve the history of sequential calculation of an unknown vector in order to visualize the convergence of the method 
#' (Данный метод необходим для сохранения истории последовательного вычисления неизвестного вектора с целью визуализации сходимости метода)
#' @param A - the original matrix of the operator equation - numeric or complex matrix (исходная матрица операторного уравнения - вещественная или комплексная)
#' @param f - bias - numeric or complex vector (вектор свободных членов вещественный или комплексный)
#' @param u - initial approximation of an unknown vector - numeric or complex vector (начальное приближение неизвестного вектора - вещественный или комплексный вектор)
#' @param eps - accuracy of calculation of the desired vector - numeric (точность вычисления искомого вектора - вещественная)
#'
#' @return result - list: 
#' num.iter - number of iterations (число итераций); 
#' var - unknown vector result (результат вычисления неизвестного вектора); 
#' var.hist - history of computing an unknown vector (история вычисления неизвестного вектора); 
#' systime.iter - system time calculation (системное время вычисления); 
#' @export
#'
#' @examples A <- diag(c(0.3, 0.4, 0.5), nrow = 3, ncol = 3)
#' f <- rnorm(3)
#' u <- rnorm(3)
#' result <- SIM.history(A = A, u = u, f = f, eps = 10e-4)
#' print(result)
#' 
#' A <- diag(c(0.5, 0.6 + 0.3i, 0.8, 0.2), nrow = 4, ncol = 4)
#' f <- rnorm(4) + 1i * rnorm(4)
#' u <- rnorm(4) + 1i * rnorm(4)
#' result <- SIM.history(A = A, u = u, f = f, eps = 10e-4)
#' print(result)
SIM.history <- function(A, f, u, eps = 10e-4) {
    # Проверка на N >= 2 - размерность
    # Все размерности совпадают
    # Все типы данных соответствуют ограничениям
    stopifnot(is.matrix(A), is.numeric(A) || is.complex(A), is.numeric(f) || is.complex(f), is.numeric(u) || is.complex(u), is.numeric(eps), nrow(A) == ncol(A), ncol(A) == length(f), length(f) == length(u), ncol(A) >= 2)
    dimA <- dim(A)[1]
    i <- 0
    u.hist <- matrix(u, nrow = dimA)
    t1 <- Sys.time()
    B <- diag(1, nrow = dimA, ncol = dimA) - A
    repeat {
        u <- B %*% u + f
        i <- i + 1
        u.hist <- cbind(u.hist, u)
        if (abs((sqrt(t(A %*% u - f) %*% (A %*% u - f))) / (sqrt(t(f) %*% f))) < eps) break
    }
    t2 <- Sys.time()
    return(list(num.iter = i, var = u, var.hist = u.hist, systime.iter = difftime(t2, t1, units = "secs")[[1]]))
}