# Setup -------------------------------------------------------------------

set.seed(123)                   
Reductor <- 10
AN <- diag(seq(0.1, 0.9, 0.8/Reductor)) 
AC <- diag(seq(0.1, 0.9, 0.8/Reductor) + 
               1i * seq(0.1, 0.9, 0.8/Reductor))
fN <- rnorm(Reductor + 1)
fC <- rnorm(Reductor + 1) + 1i * rnorm(Reductor + 1)
uN <- rnorm(Reductor + 1)
uC <- rnorm(Reductor + 1) + 1i * rnorm(Reductor + 1)
eps <- c(10e-3, 10e-4, 10e-5, 10e-6, 10e-7, 10e-8, 10e-09)

# testing and checking ----------------------------------------------------

testthat::test_that("Accuracy and numeric/complex type check", {
    for (e in eps) {
        # N/N/N
        result <- IMSSLAER::GDM(A = AN, f = fN, u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fN) %*% 
                          Conj(AN %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/N/C
        result <- IMSSLAER::GDM(A = AN, f = fN, u = uC, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fN) %*% 
                          Conj(AN %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/C/N
        result <- IMSSLAER::GDM(A = AN, f = fC, u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fC) %*% 
                          Conj(AN %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/N/N
        result <- IMSSLAER::GDM(A = AC, f = fN, u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fN) %*% 
                          Conj(AC %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # C/C/N
        result <- IMSSLAER::GDM(A = AC, f = fC, u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fC) %*% 
                          Conj(AC %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/N/C
        result <- IMSSLAER::GDM(A = AC, f = fN, u = uC, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fN) %*% 
                          Conj(AC %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/C/C
        result <- IMSSLAER::GDM(A = AN, f = fC, u = uC, eps = e, 
                                iterations = 40000)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fC) %*% 
                          Conj(AN %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/C/C
        result <- IMSSLAER::GDM(A = AC, f = fC, u = uC, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fC) %*% 
                          Conj(AC %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
    }
})

testthat::test_that("Lengths and DIMs check", {
    testthat::expect_length(
        IMSSLAER::GDM(A = AN, u = uN, f = fN, eps = eps[1]), 
        n = length(uN)
        )
    for (e in eps) {
        for (i in 1:(Reductor + 1)) {
            testthat::expect_error(
                IMSSLAER::GDM(A = AN, u = uN, f = fN[-i], 
                              eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN, u = uN[-i], 
                              f = fN[-i], eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN, u = uN[-i], 
                              f = fN, eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN[-i, ], u = uN[-i], 
                              f = fN, eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN[-i, ], u = uN, 
                              f = fN, eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN[-i, ], u = uN, 
                              f = fN[-i], eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN[-i, ][, -i], u = uN, 
                              f = fN[-i], eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN[-i, ][, -i], u = uN[-i], 
                              f = fN, eps = e)
                )
            testthat::expect_error(
                IMSSLAER::GDM(A = AN[-i, ], u = uN[-i], 
                              f = fN[-i], eps = e)
                )
            testthat::expect(
                IMSSLAER::GDM(A = AN[-i, ][, -i], u = uN[-i], 
                              f = fN[-i], eps = e)
                )
        }
    }
})

testthat::test_that("Error stuctures check", {
    set.seed(123)
    for (e in eps) {    
        M <- rnorm(Reductor + 1)
        testthat::expect_error(IMSSLAER::GDM(M, fN, uN, e))
        M <- list(vec = rnorm(Reductor + 1))
        testthat::expect_error(IMSSLAER::GDM(AN, M, uN, e))
        testthat::expect_error(IMSSLAER::GDM(AN, fN, M, e))
        testthat::expect_error(IMSSLAER::GDM(M, fN,  uN, e))
        M <- diag(letters[1:(Reductor + 1)])
        testthat::expect_error(IMSSLAER::GDM(M, fN, uN, e))
        M <- diag((1:(Reductor + 1)))
        testthat::expect_error(IMSSLAER::GDM(AN, M, uN, e))
        testthat::expect_error(IMSSLAER::GDM(AN, fN, M, e))
    }
    rm(M)
})

testthat::test_that("GDM works", {
    M <- matrix(rnorm((Reductor + 1) ^ 2), 
                nrow = Reductor + 1, ncol = Reductor + 1)
    testthat::expect_length(
        IMSSLAER::GDM(M, fC, uN, eps = 10e-6, 
                      iterations = 10000), 
        n = Reductor + 1
        )
})

# Clear -------------------------------------------------------------------

rm(list = ls())
