# Setup -------------------------------------------------------------------

set.seed(123) 
Reductor <- 10
AN <- diag(seq(15.1, 199.9, 184.8/Reductor)) 
AC <- diag(seq(15.1, 199.9, 184.8/Reductor) + 
               1i * seq(1.1, 5.9, 4.8/Reductor))
fN <- rnorm(Reductor + 1, 6)
fC <- rnorm(Reductor + 1, 5) + 1i * rnorm(Reductor + 1, 5)
uN <- rnorm(Reductor + 1, 7)
uC <- rnorm(Reductor + 1, 6) + 1i * rnorm(Reductor + 1, 8)
eps <- c(10e-3, 10e-4, 10e-5, 10e-6, 10e-7, 10e-8, 10e-9)

# testing and checking ----------------------------------------------------

testthat::test_that("Accuracy and numeric/complex type check", {
    for (e in eps) {
        # N/N/N
        result <- IMSSLAER::BiCGMStab(A = AN, f = fN, 
                                      u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fN) %*% 
                          Conj(AN %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/N/C
        result <- IMSSLAER::BiCGMStab(A = AN, f = fN, 
                                      u = uC, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fN) %*% 
                          Conj(AN %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        #N/C/N
        result <- IMSSLAER::BiCGMStab(A = AN, f = fC, 
                                      u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fC) %*%
                          Conj(AN %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/N/N
        result <- IMSSLAER::BiCGMStab(A = AC, f = fN, 
                                      u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fN) %*% 
                          Conj(AC %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # C/C/N
        result <- IMSSLAER::BiCGMStab(A = AC, f = fC, 
                                      u = uN, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fC) %*% 
                          Conj(AC %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/N/C
        result <- IMSSLAER::BiCGMStab(A = AC, f = fN, 
                                      u = uC, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fN) %*% 
                          Conj(AC %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/C/C
        result <- IMSSLAER::BiCGMStab(A = AN, f = fC, 
                                      u = uC, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fC) %*% 
                          Conj(AN %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # # C/C/C
        result <- IMSSLAER::BiCGMStab(A = AC, f = fC, 
                                      u = uC, eps = e)
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fC) %*% 
                          Conj(AC %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
    }
})

testthat::test_that("Lengths and DIMs check", {
    testthat::expect_length(
        IMSSLAER::BiCGMStab(A = AN, u = uN, f = fN, 
                            eps = eps[1]),
        n = length(uN)
        )
    for (e in eps) {
        for (i in 1:(Reductor + 1)) {
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN, u = uN, 
                                    f = fN[-i], eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN, u = uN[-i], 
                                    f = fN[-i], eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN, u = uN[-i], 
                                    f = fN, eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN[-i, ],
                                    u = uN[-i], f = fN, 
                                    eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN[-i, ], 
                                    u = uN, f = fN, 
                                    eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN[-i, ], 
                                    u = uN, f = fN[-i], 
                                    eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN[-i, ][, -i], 
                                    u = uN, f = fN[-i], 
                                    eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN[-i, ][, -i], 
                                    u = uN[-i], f = fN, 
                                    eps = e)
                )
            testthat::expect_error(
                IMSSLAER::BiCGMStab(A = AN[-i, ], 
                                    u = uN[-i], f = fN[-i], 
                                    eps = e)
                )
            testthat::expect_length(
                IMSSLAER::BiCGMStab(A = AN[-i, ][, -i],
                                    u = uN[-i], f = fN[-i], 
                                    eps = e), 
                n = length(uN[-i])
                )
        }
    }
})

testthat::test_that("Error stuctures check", {
    set.seed(123)
    for (e in eps) {    
        M <- rnorm(Reductor + 1)
        testthat::expect_error(
            IMSSLAER::BiCGMStab(M, fN, uN, e))
        M <- list(vec = rnorm(Reductor + 1))
        testthat::expect_error(
            IMSSLAER::BiCGMStab(AN, M, uN, e))
        testthat::expect_error(
            IMSSLAER::BiCGMStab(AN, fN, M, e))
        testthat::expect_error(
            IMSSLAER::BiCGMStab(M, fN,  uN, e))
        M <- diag(letters[1:(Reductor + 1)])
        testthat::expect_error(
            IMSSLAER::BiCGMStab(M, fN, uN, e))
        M <- diag((1:(Reductor + 1)))
        testthat::expect_error(
            IMSSLAER::BiCGMStab(AN, M, uN, e))
        testthat::expect_error(
            IMSSLAER::BiCGMStab(AN, fN, M, e))
    }
    rm(M)
})

testthat::test_that("BiCGMStab works", {
    M <- matrix(rnorm((Reductor + 1) ^ 2), 
                nrow = Reductor + 1, ncol = Reductor + 1)
    testthat::expect_length(
        IMSSLAER::BiCGMStab(M, fC, uN, eps = 10e-6, 
                            iterations = 10000), 
        n = Reductor + 1
        )
})

# Clear -------------------------------------------------------------------

rm(list = ls())
