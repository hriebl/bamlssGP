library(bamlss)
library(numDeriv)
library(purrr)

test_score <- function(wrt, family, par, tolerance = 0.001) {
  link <- make.link(family$links[[wrt]])

  fun <- function(x) {
    par[[wrt]] <- link$linkinv(x)
    sum(family$d(par = par))
  }

  expect_equal(
    object = family$score[[wrt]](par = par),
    expected = grad(fun, link$linkfun(par[[wrt]])),
    tolerance = tolerance
  )
}

test_hess <- function(wrt, family, par, tolerance = 0.001) {
  link <- make.link(family$links[[wrt]])

  fun <- function(x) {
    par[[wrt]] <- link$linkinv(x)
    sum(family$score[[wrt]](par = par))
  }

  expect_equal(
    object = mean(family$hess[[wrt]](par = par)),
    expected = -mean(grad(fun, link$linkfun(par[[wrt]]))),
    tolerance = tolerance
  )
}
