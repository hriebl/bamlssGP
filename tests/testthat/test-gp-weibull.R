set.seed(1337)

q_fun <- function(s) {
  ifelse(s < 30, 0.1 + 0.9 / 30 * s, 1)
}

N <- 30 # 170 in the application
n <- 30 # ~180 in the application

limit <- 2000
shape <- 3.5
scale <- 100
linear1 <- 1
stddev <- 50
range <- 2

s <- seq.int(1, 180, length.out = n)
D <- as.matrix(dist(s))

q <- q_fun(s)
qq <- tcrossprod(q)

mean <- m_fun(s, limit, shape, scale)
covariance <- c_fun_matern(D, stddev, range, qq)
y <- rmvnorm_c(N, mean, covariance)

Z <- rmvnorm_c(N, 0, c_fun_gauss(D, 80, 28, qq))

y <- array_branch(y + linear1 * Z, 2)
Z <- map(array_branch(Z, 2), matrix, ncol = 1)
s <- rep.int(list(s), N)
D <- rep.int(list(D), N)
q <- rep.int(list(q), N)

initial_par <- data.frame(
  limit = rep.int(limit, N),
  shape = rep.int(shape, N),
  scale = rep.int(scale, N),
  linear1 = rep.int(linear1, N),
  stddev = rep.int(stddev, N),
  range = rep.int(range, N)
)

family <- gp_weibull_bamlss(y, s, Z, D, q, initial_par)

test_that("bamlss runs", {
  expect_error(
    object = bamlss(
      formula = list(rep.int(0, N) ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1),
      family = family,
      n.iter = 30,
      burnin = 0
    ),
    regexp = NA
  )
})

test_that("score of limit parameter is correct", {
  test_score("limit", family, initial_par)
})

test_that("score of shape parameter is correct", {
  test_score("shape", family, initial_par)
})

test_that("score of scale parameter is correct", {
  test_score("scale", family, initial_par)
})

test_that("score of linear1 parameter is correct", {
  test_score("linear1", family, initial_par)
})

test_that("score of stddev parameter is correct", {
  test_score("stddev", family, initial_par)
})

test_that("score of range parameter is correct", {
  test_score("range", family, initial_par)
})

test_that("hessian of limit parameter is correct", {
  test_hess("limit", family, initial_par, tolerance = 30)
})

test_that("hessian of shape parameter is correct", {
  test_hess("shape", family, initial_par)
})

test_that("hessian of scale parameter is correct", {
  test_hess("scale", family, initial_par, tolerance = 20)
})

test_that("hessian of linear1 parameter is correct", {
  test_hess("linear1", family, initial_par)
})

test_that("hessian of stddev parameter is correct", {
  test_hess("stddev", family, initial_par, tolerance = 5)
})

test_that("hessian of range parameter is correct", {
  test_hess("range", family, initial_par, tolerance = 5)
})
