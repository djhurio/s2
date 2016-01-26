#### S^2 estimators


# Packages ####
require(vardpoor)


# Definition of S^2 estimators ####

# s2a -- sample variance -- equal to var()

s2a <- function(y) {
  n <- length(y)
  s2 <- (sum(y ^ 2) - sum(y) ^ 2 / n) / (n - 1)
  return(s2)
}


# s2b

s2b <- function(y, w = rep(1, length(y))) {
  N <- sum(w)
  s2 <- (sum(y ^ 2 * w) - sum(y * w) ^ 2 / N) / (N - 1)
  return(s2)
}


# s2c

s2c <- function(y, w = rep(1, length(y))) {
  n <- length(y)
  N <- sum(w)
  s2 <- (sum(y ^ 2 * w) - sum(y * w) ^ 2 / N) / N * n / (n - 1)
  return(s2)
}


# s2d

s2d <- function(y, w = rep(1, length(y)),
                H = rep(1, length(y)), PSU = seq_along(y)) {
  N <- sum(w)
  n <- length(y)
  VarY <- vardom(Y = y, H = H, PSU = PSU, w_final = w)$all_result$var
  s2 <- (sum(y ^ 2 * w) - (sum(y * w) ^ 2 - VarY) / N) / (N - 1)
  return(s2)
}



# Test ####

# Test 1

y <- runif(10)

# var(y)
# s2a(y)
# s2b(y)
# s2c(y)
# s2d(y)

all.equal(rep(var(y), 4), c(s2a(y), s2b(y), s2c(y), s2d(y)))


# Test 2

y <- runif(10)
w <- rep(5, 10)

# var(y)
# s2a(y)
# s2b(y, w)
# s2c(y, w)
# s2d(y, w)

all.equal(rep(var(y), 3), c(s2a(y), s2c(y, w), s2d(y, w)))
