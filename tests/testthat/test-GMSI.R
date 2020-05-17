# Setup -------------------------------------------------------------------

set.seed(123)
Reductor <- 10
AN <- diag(seq(0.1, 14.9, 14.8 / Reductor))
AC <- diag(seq(0.1, 14.9, 14.8 / Reductor) + 
               1i * seq(0.1, 14.9, 14.8 / Reductor))
fN <- rnorm(Reductor + 1, 6)
fC <- rnorm(Reductor + 1, 5) + 1i * rnorm(Reductor + 1, 5)
uN <- rnorm(Reductor + 1, 7)
uC <- rnorm(Reductor + 1, 6) + 1i * rnorm(Reductor + 1, 8)
eps <- c(10e-3, 10e-4, 10e-5, 10e-6, 10e-7, 10e-8, 10e-9)
# Testing -----------------------------------------------------------------

testthat::test_that("muFind Works", {
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(0.1, 0.9, 0.1) + 
                             1i * seq(0.1, 0.9, 0.1), 
                         draw = FALSE), 0.5 + 0.5i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(0.1, 1.5, 0.1) + 
                             1i * seq(0.1, 1.5, 0.1), 
                         draw = FALSE), 0.8 + 0.8i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(1, 1.5, 0.1) + 
                             1i * seq(1, 1.5, 0.1), 
                         draw = FALSE), 1.25 + 1.25i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(1, 2, 0.1), 
                         draw = FALSE), 1.5 + 0i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(-1, -2, -0.1), 
                         draw = FALSE), -1.5 + 0i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(-1, -1.5, -0.1) + 
                             1i * seq(1, 1.5, 0.1), 
                         draw = FALSE), -1.25 + 1.25i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(-1, -1.5, -0.1) + 
                             1i * seq(-1, -1.5, -0.1), 
                         draw = FALSE), -1.25 + -1.25i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = seq(1, 1.5, 0.1) + 
                             1i * seq(-1, -1.5, -0.1), 
                         draw = FALSE), 1.25 + -1.25i
    )
    testthat::expect_equal(
        IMSSLAER::muFind(lambs = c(5 + 5i, 5 - 5i, 4 + 0i), 
                         draw = FALSE), 10 + 0i
    )
})

testthat::test_that("GMSI.mu works", {
    testthat::expect_length(GMSI.mu(AC, fC, uC, 
                                    lambs = diag(AC)), 
                            n = length(uC))
})

testthat::test_that("Accuracy and numeric/complex type check", {
    for (e in eps) {
        # N/N/N
        result <- IMSSLAER::GMSI.mu(A = AN, f = fN, u = uN, 
                                    eps = e, lambs = diag(AN))
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fN) %*% 
                          Conj(AN %*% result - fN))) /
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/N/C
        result <- IMSSLAER::GMSI.mu(A = AN, f = fN, u = uC, 
                                    eps = e, lambs = diag(AN))
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fN) %*% 
                          Conj(AN %*% result - fN))) /
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/C/N
        result <- IMSSLAER::GMSI.mu(A = AN, f = fC, u = uN, 
                                    eps = e, lambs = diag(AN))
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fC) %*% 
                          Conj(AN %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/N/N
        result <- IMSSLAER::GMSI.mu(A = AC, f = fN, u = uN, 
                                    eps = e, lambs = diag(AC))
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fN) %*% 
                          Conj(AC %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # C/C/N
        result <- IMSSLAER::GMSI.mu(A = AC, f = fC, u = uN, 
                                    eps = e, lambs = diag(AC))
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fC) %*% 
                          Conj(AC %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/N/C
        result <- IMSSLAER::GMSI.mu(A = AC, f = fN, u = uC, 
                                    eps = e, lambs = diag(AC))
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fN) %*% 
                          Conj(AC %*% result - fN))) / 
                    (sqrt(t(fN) %*% Conj(fN)))) < e
            )
        # N/C/C
        result <- IMSSLAER::GMSI.mu(A = AN, f = fC, u = uC, 
                                    eps = e, lambs = diag(AN))
        testthat::expect_true(
            abs((sqrt(t(AN %*% result - fC) %*% 
                          Conj(AN %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
        # C/C/C
        result <- IMSSLAER::GMSI.mu(A = AC, f = fC, u = uC, 
                                    eps = e, lambs = diag(AC))
        testthat::expect_true(
            abs((sqrt(t(AC %*% result - fC) %*% 
                          Conj(AC %*% result - fC))) / 
                    (sqrt(t(fC) %*% Conj(fC)))) < e
            )
    }
})

testthat::test_that("Lengths and DIMs check", {
    testthat::expect_length(
        IMSSLAER::GMSI.mu(A = AN, u = uN, f = fN, eps = eps[1], 
                          lambs = diag(AN)),
        n = length(uN)
        )
    for (e in eps) {
        for (i in 1:(Reductor + 1)) {
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN, u = uN, f = fN[-i], 
                                  eps = e, lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN, u = uN[-i], 
                                  f = fN[-i], eps = e, 
                                  lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN, u = uN[-i], f = fN, 
                                  eps = e, lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN[-i, ], u = uN[-i], 
                                  f = fN, eps = e, 
                                  lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN[-i, ], u = uN, f = fN, 
                                  eps = e, lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN[-i, ], u = uN, 
                                  f = fN[-i], eps = e, 
                                  lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN[-i, ][, -i], u = uN, 
                                  f = fN[-i], eps = e, 
                                  lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN[-i, ][, -i], u = uN[-i], 
                                  f = fN, eps = e, 
                                  lambs = diag(AN))
                )
            testthat::expect_error(
                IMSSLAER::GMSI.mu(A = AN[-i, ], u = uN[-i], 
                                  f = fN[-i], eps = e, 
                                  lambs = diag(AN))
                )
            testthat::expect_length(
                IMSSLAER::GMSI.mu(A = AN[-i, ][, -i], u = uN[-i], 
                                  f = fN[-i], eps = e, 
                                  lambs = diag(AN[-i, ][, -i])), 
                n = length(as.vector(uN[-i]))
                )
        }
    }
})
# Clearing ----------------------------------------------------------------

rm(list = ls())
