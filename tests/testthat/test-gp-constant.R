set.seed(1337)

N <- 30
n <- 30

mean <- 0
stddev <- exp(0)
range <- exp(0)

s <- seq.int(0, 1, length.out = n)
D <- as.matrix(dist(s))

y <- rmvnorm_c(N, mean, c_fun_matern(D, stddev, range))

y <- array_branch(y, 2)
D <- rep.int(list(D), N)

initial_par <- data.frame(
  mean = rep.int(mean, N),
  stddev = rep.int(stddev, N),
  range = rep.int(range, N)
)

family <- gp_constant_bamlss(y, D, initial_par)

test_that("bamlss runs", {
  expect_error(
    object = bamlss(
      formula = list(rep.int(0, N) ~ 1, ~ 1, ~ 1),
      family = family,
      n.iter = 30,
      burnin = 0
    ),
    regexp = NA
  )
})

test_that("score of mean parameter is correct", {
  test_score("mean", family, initial_par)
})

test_that("score of stddev parameter is correct", {
  test_score("stddev", family, initial_par)
})

test_that("score of range parameter is correct", {
  test_score("range", family, initial_par)
})

test_that("hessian of mean parameter is correct", {
  test_hess("mean", family, initial_par)
})

test_that("hessian of stddev parameter is correct", {
  # increase N to decrease tolerance
  test_hess("stddev", family, initial_par, tolerance = 5)
})

test_that("hessian of range parameter is correct", {
  # increase N to decrease tolerance
  test_hess("range", family, initial_par, tolerance = 10)
})
