
# prime numbers funcs -----------------------------------------------------

Primes.is_prime <- Vectorize(function(number) {
    if (abs(number) == 1 || abs(number) == 2) return(T)
    is_prime <- T
    for (i in 2:ceiling(sqrt(number))) {
        if (number %% i == 0) {
            is_prime <- F
            break
        }
    }
    return(is_prime)
})


Primes.euler_fun <- Vectorize(function(number) {
    result <- number
    i <- 2
    while (i * i <= number) {
        if (number %% i == 0) {
            while (number %% i == 0) {
                number = number / i
            }
            result = result - result / i
        }
        i = i + 1
    }
    if (number > 1) {
        result = result - result / number
    }
    return(result)
    
})

# Tests -------------------------------------------------------------------
# 
# 
# testthat::test_that("is_prime fucntion unit tests", {
#     testthat::expect_equal(Primes.is_prime(1), T)
#     testthat::expect_equal(Primes.is_prime(-1), T)
#     testthat::expect_equal(Primes.is_prime(2), T)
#     testthat::expect_equal(Primes.is_prime(-2), T)
#     testthat::expect_equal(Primes.is_prime(3), T)
#     testthat::expect_equal(Primes.is_prime(5), T)
#     for (i in 1:4) {
#         testthat::expect_equal(Primes.is_prime(4 * i), F)
#     }
#     
# })
# 
# testthat::test_that("euler_fun function unit tests", {
#     eulers <- c(1, 1, 2, 2, 4, 2, 6, 4, 6,
#                 4, 10, 4, 12, 6, 8, 8, 16, 6, 18,
#                 8, 12, 10, 22, 8, 20, 12, 18, 12, 28, 8)
#     for (i in 1:30) {
#         testthat::expect_equal(Primes.euler_fun(i), eulers[i])
#     }
# })
