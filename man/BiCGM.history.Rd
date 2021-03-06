% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BiCGM.R
\name{BiCGM.history}
\alias{BiCGM.history}
\title{Biconjugate gradient method [BiCGM]
(Метод бисопряженных градиентов)}
\usage{
BiCGM.history(A, f, u, eps = 0.001, iterations = 10000)
}
\arguments{
\item{A}{- the original matrix of the operator equation - 
numeric matrix only (исходная матрица операторного уравнения 
- вещественная только)}

\item{f}{- bias - numeric or complex vector (вектор свободных 
членов вещественный или комплексный)}

\item{u}{- initial approximation of an unknown vector - 
numeric or complex vector (начальное приближение неизвестного 
вектора - вещественный или комплексный вектор)}

\item{eps}{- accuracy of calculation of the desired vector 
- numeric (точность вычисления искомого вектора - вещественная)}

\item{iterations}{- the upper limit on the number of 
iterations when the method diverges (ограничение сверху 
на число итераций при расхождении метода)}
}
\value{
result - list: 
num.iter - number of iterations (число итераций); 
var - unknown vector result (результат вычисления неизвестного 
вектора); 
var.hist - history of computing an unknown vector 
(история вычисления неизвестного вектора); 
systime.iter - system time calculation (системное время 
вычисления);
}
\description{
Non-stationary iterative numerical method 
for solving SLAEs of the Krylov type. It is a 
generalization of the conjugate gradient method.
(Нестационарный итерационный численный метод решения 
СЛАУ крыловского типа. Является обобщением метода 
сопряжённых градиентов.)
}
\details{
This method is necessary to preserve the history 
of sequential calculation of an unknown vector in order 
to visualize the convergence of the method 
(Данный метод необходим для сохранения истории 
последовательного вычисления неизвестного вектора с целью 
визуализации сходимости метода)
}
\examples{
A <- diag(rnorm(5, 2), nrow = 5, ncol = 5)
u <- rnorm(5, 12)
f <- rnorm(5, 17)
print(BiCGM.history(A, f, u, iterations = 10000))
}
