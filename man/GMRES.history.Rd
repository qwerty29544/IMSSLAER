% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GMRES.R
\name{GMRES.history}
\alias{GMRES.history}
\title{Generalized minimal residual algorithm history [GMRES.history]
(Обобщенный метод минимальных невязок)}
\usage{
GMRES.history(A, f, u, layers = 2, eps = 0.001, iterations = 10000)
}
\arguments{
\item{A}{- the original matrix of the operator equation - numeric 
or complex matrix (исходная матрица операторного уравнения 
- вещественная или комплексная)}

\item{f}{- bias - numeric or complex vector (вектор свободных 
членов вещественный или комплексный)}

\item{u}{- initial approximation of an unknown vector - 
numeric or complex vector (начальное приближение неизвестного вектора 
- вещественный или комплексный вектор)}

\item{layers}{- Krylov subspace order (порядок подпространства 
Крылова)}

\item{eps}{- accuracy of calculation of the desired vector 
- numeric (точность вычисления искомого вектора - вещественная)}

\item{iterations}{- the upper limit on the number of iterations 
when the method diverges (ограничение сверху на число итераций 
при расхождении метода)}
}
\value{
result - list: 
num.iter - number of iterations (число итераций); 
var - unknown vector result (результат вычисления неизвестного вектора); 
var.hist - history of computing an unknown vector (история вычисления 
неизвестного вектора); 
systime.iter - system time calculation (системное время вычисления);
}
\description{
Non-stationary iterative numerical method for solving 
systems of linear algebraic equations. The projection of the vector 
onto the Krylov subspace of arbitrary order is used as a residual.
(Нестационарный итерационный численный метод решения систем 
линейных алгебраических уравнений. В качестве невязки используется 
проекция вектора на подпространство Крылова произвольного порядка.)
}
\details{
This method is necessary to preserve the history of sequential 
calculation of an unknown vector in order to visualize the 
convergence of the method 
(Данный метод необходим для сохранения истории последовательного 
вычисления неизвестного вектора с целью визуализации сходимости метода)
}
\examples{
A <- diag(rnorm(25, 5) + 1i * rnorm(25, 1), ncol = 25, nrow = 25)
f <- rnorm(25, 2)
u <- rnorm(25) 
print(GMRES.history(A, f, u, layers = 5, eps = 10e-07))
print(IMRES.history(A, f, u, eps = 10e-7))
}
